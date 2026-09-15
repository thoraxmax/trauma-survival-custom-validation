# 02_RQ2_time_to_onset.R
#
# =============================================================================
# ANALYSES INCLUDED
#
# Statistical model = SM
# Statistical package version/name = SP
# Independent / dependent variables = IV/DV
# Model formula = MF
# Motivation = M
# =============================================================================

# 1. Latency descriptives from study entry, most recent trauma, and worst trauma.
# SM: Descriptive statistics; SP: base R; IV/DV: time origin / time to first detected OCD; MF: none; M: summarize observed onset timing.

# 2. Midpoint-onset sensitivity using the last negative/first positive interval.
# SM: Alternative event-time definition; SP: base R + survival; IV/DV: onset definition / midpoint-estimated OCD onset; MF: midpoint = (last negative + first positive) / 2; M: assess interval-censoring uncertainty.

# 3. Restricted mean OCD-free time (RMST) with participant-bootstrap CIs.
# SM: Delayed-entry Kaplan-Meier RMST with participant bootstrap; SP: survival::survfit; IV/DV: time origin/group / OCD-free time; MF: Surv(entry, duration, event) ~ 1; M: summarize OCD-free time without a PH assumption.

# 4. Delayed-entry Cox model for female versus male sex.
# SM: Delayed-entry attained-age Cox PH; SP: survival::coxph; IV/DV: sex / first detected OCD; MF: Surv(age_entry, age_exit, event) ~ sex; M: estimate the sex association while accounting for delayed entry.

# 5. One-time sex PH assessment before subsequent Cox models.
# SM: Initial attained-age Cox model plus PH diagnostic; SP: survival::coxph + survival::cox.zph; IV/DV: sex / OCD; MF: Surv(age_entry, age_exit, event) ~ sex, followed by consistent sex adjustment or strata(sex) in all subsequent Cox models; M: avoid repeated data-dependent switching of sex handling across models.

# 6. Rank-transformed Schoenfeld PH tests using survival::cox.zph.
# SM: Scaled Schoenfeld residual PH diagnostic; SP: survival::cox.zph; IV/DV: fitted Cox terms / residual-time association; MF: cox.zph(fit, transform = "rank"); M: assess the proportional-hazards assumption.

# 7. Age-varying coefficient models when a substantive fixed predictor is PH-flagged.
# SM: Extended Cox model with tt(); SP: survival::coxph; IV/DV: PH-flagged exposure x log attained age / OCD; MF: Surv(...) ~ exposure + tt(exposure) [+ strata(sex)]; M: estimate age-dependent effects when PH is violated.

# 8. Start-stop proximal-risk models for 3-12 months after most recent/worst trauma.
# SM: Start-stop Cox with participant-clustered robust SE; SP: survival::coxph; IV/DV: within-cutoff + age at trauma / OCD; MF: Surv(start, stop, event) ~ within_cutoff + age_at_trauma + strata(sex) + cluster(record_id); M: test whether OCD hazard is elevated in the early post-trauma period.

# 9. Midpoint-onset versions of the 3-12-month proximal-risk models.
# SM: Start-stop Cox using midpoint event times; SP: survival::coxph; IV/DV: within-cutoff + age at trauma / midpoint OCD onset; MF: same as Analysis 8 using midpoint exit time; M: assess robustness to uncertainty in OCD onset timing.

# 10. Exploratory pooled sex interactions for age at most recent/worst trauma.
# SM: Pooled delayed-entry Cox interaction; SP: survival::coxph; IV/DV: centered age at trauma x female / OCD; MF: Surv(...) ~ age_trauma_c + age_trauma_c_x_female + strata(sex); M: formally test sex effect modification without splitting the sample.

# 11. Prespecified 6-month proximal-risk-by-sex interaction sensitivity models.
# SM: Start-stop pooled Cox interaction with clustered SE; SP: survival::coxph; IV/DV: within-cutoff x female / OCD; MF: Surv(...) ~ within_cutoff + within_cutoff_x_female + age_at_trauma + strata(sex) + cluster(record_id); M: test whether 6-month proximal risk differs by sex.

# 12. Time-updated new lifetime-worst-trauma model.
# SM: Start-stop time-varying Cox with participant-clustered robust SE; SP: survival::coxph; IV/DV: new-worst-trauma time-varying indicator / OCD; MF: Surv(start, stop, event) ~ new_worst_trauma_tv [+ trauma_type_count] + sex handling + cluster(record_id); M: preserve exposure timing and avoid immortal-time bias.
# Standard PH tests are not applied to the explicitly time-varying cutoff models.

resolve_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  frames <- sys.frames()
  source_files <- vapply(
    frames,
    function(frame) {
      value <- frame$ofile
      if (is.null(value) || !length(value)) "" else as.character(value[[1L]])
    },
    character(1)
  )
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files)) {
    return(dirname(normalizePath(tail(source_files, 1L), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

script_dir <- resolve_script_dir()
Sys.setenv(SURVIVAL_PACKAGE_DIR = script_dir)
helper_path <- file.path(script_dir, "custom_functions.R")
if (!file.exists(helper_path)) {
  stop("Cannot find custom_functions.R beside this script: ", helper_path, call. = FALSE)
}
source(helper_path, local = FALSE)
assert_helper_api(21L, c("progress_log", "format_elapsed", "read_survival_data", "save_cox_zph_plots", "bootstrap_km", "save_apa_bundle"))
SCRIPT_STARTED <- Sys.time()

# =============================================================================
# LOCAL COX ENGINE -- model fitting is defined in this RQ script, not hidden in
# the shared helper. The helper is used only for I/O, extraction, diagnostics,
# bootstrap utilities, and reporting.
# =============================================================================
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("This analysis requires the R package 'survival'.", call. = FALSE)
}

fit_age_scale_cox <- function(data, covariates, strata = character(), cluster = NULL, model_name = "Cox model") {
  required <- unique(c("age_entry", "age_exit", "event", covariates, strata, cluster))
  required <- required[!is.na(required) & nzchar(required)]
  require_columns(data, required, model_name)
  fit_data <- complete_model_data(data, required)
  fit_data <- fit_data[
    fit_data$event %in% c(0, 1) & is.finite(fit_data$age_entry) &
      is.finite(fit_data$age_exit) & fit_data$age_exit > fit_data$age_entry,
    , drop = FALSE
  ]
  if (!nrow(fit_data) || sum(fit_data$event) == 0L) {
    stop(model_name, " has no estimable rows/events after complete-case filtering.", call. = FALSE)
  }
  rhs <- c(covariates, if (length(strata)) sprintf("strata(%s)", strata) else character())
  if (!is.null(cluster) && length(cluster) && nzchar(cluster)) rhs <- c(rhs, sprintf("cluster(%s)", cluster))
  if (!length(rhs)) rhs <- "1"
  form <- stats::as.formula(paste("survival::Surv(age_entry, age_exit, event) ~", paste(rhs, collapse = " + ")))
  progress_log("Cox model", paste0(model_name, "; survival::coxph; formula: ", paste(deparse(form), collapse = " ")))
  fit <- survival::coxph(
    form, data = fit_data, ties = "efron", x = TRUE, y = TRUE,
    model = TRUE, singular.ok = TRUE
  )
  attr(fit, "analysis_data") <- fit_data
  attr(fit, "model_name") <- model_name
  fit
}

progress_log("Statistical engine", paste0("R package survival ", as.character(utils::packageVersion("survival")), "; Cox models use survival::coxph; PH diagnostics use survival::cox.zph"))

fit_extended_age_effect <- function(
  data, exposure, adjustment_terms=character(), strata=character(),
  model_name=paste("Extended Cox:", exposure), report_ages=REPORT_AGES
) {
  required <- unique(c("age_entry","age_exit","event",exposure,adjustment_terms,strata))
  required <- required[!is.na(required) & nzchar(required)]
  require_columns(data, required, model_name)
  d <- complete_model_data(data, required)
  d <- d[d$event %in% c(0,1) & d$age_exit > d$age_entry,,drop=FALSE]
  event_ages <- d$age_exit[d$event==1]
  if (length(event_ages) < 3L) stop(model_name, " has too few events.", call.=FALSE)
  center_log_age <- mean(log(event_ages), na.rm=TRUE)
  rhs <- c(exposure, adjustment_terms, sprintf("tt(%s)", exposure))
  if (length(strata) && !is.na(strata[[1L]]) && nzchar(strata[[1L]])) rhs <- c(rhs, sprintf("strata(%s)", strata[[1L]]))
  form <- stats::as.formula(paste("survival::Surv(age_entry, age_exit, event) ~", paste(rhs, collapse=" + ")))
  progress_log("Cox model", paste0(model_name, "; survival::coxph time-varying coefficient; formula: ", paste(deparse(form),collapse=" ")))
  fit <- survival::coxph(form, data=d, ties="efron", x=TRUE, y=TRUE, singular.ok=TRUE,
    tt=function(x,t,...) x * (log(pmax(t,1e-8)) - center_log_age))
  attr(fit,"model_name") <- model_name; attr(fit,"analysis_data") <- d
  b <- stats::coef(fit); main_name <- exposure; tt_name <- paste0("tt(",exposure,")")
  if (!(main_name %in% names(b)) || !(tt_name %in% names(b))) stop("Could not identify extended-model coefficients for ", exposure, call.=FALSE)
  V <- stats::vcov(fit)[c(main_name,tt_name),c(main_name,tt_name),drop=FALSE]
  ages <- report_ages[report_ages >= min(d$age_entry,na.rm=TRUE) & report_ages <= max(d$age_exit,na.rm=TRUE)]
  if (!length(ages)) ages <- unique(round(stats::quantile(event_ages,c(.25,.5,.75),na.rm=TRUE),1))
  age_rows <- lapply(ages,function(age){
    contrast <- c(1,log(age)-center_log_age); log_hr <- sum(contrast*b[c(main_name,tt_name)]); se <- sqrt(as.numeric(t(contrast)%*%V%*%contrast))
    data.frame(Model=model_name,Exposure=exposure,Age=age,log_HR=log_hr,SE=se,HR=exp(log_hr),
      CI_lower=exp(log_hr-stats::qnorm(.975)*se),CI_upper=exp(log_hr+stats::qnorm(.975)*se),
      z=log_hr/se,p=2*stats::pnorm(abs(log_hr/se),lower.tail=FALSE),Center_log_age=center_log_age,N=fit$n,Events=fit$nevent,stringsAsFactors=FALSE)
  })
  list(fit=fit,age_specific=do.call(rbind,age_rows),interaction=extract_cox_terms(fit,model_name),center_log_age=center_log_age)
}


# RQ2 RMST uses the same delayed-entry Kaplan-Meier estimator directly.
make_origin_data <- function(data, origin) {
  if (origin == "Study entry") {
    out <- data.frame(record_id=data$record_id,sex=data$sex,event=data$event,entry=0,duration=data$time_since_start/12)
  } else if (origin == "Most recent trauma") {
    out <- data.frame(record_id=data$record_id,sex=data$sex,event=data$event,entry=data$months_since_most_recent_at_baseline/12,duration=data$time_since_most_recent_trauma/12)
  } else if (origin == "Worst trauma") {
    out <- data.frame(record_id=data$record_id,sex=data$sex,event=data$event,entry=data$months_since_worst_at_baseline/12,duration=data$time_since_worst_trauma/12)
  } else stop("Unknown time origin: ",origin,call.=FALSE)
  out <- out[stats::complete.cases(out[,c("event","entry","duration")]),,drop=FALSE]
  out[out$event %in% c(0,1) & out$entry>=0 & out$duration>out$entry,,drop=FALSE]
}

km_fit <- function(data) {
  survival::survfit(survival::Surv(entry,duration,event) ~ 1, data=data, conf.type="log-log")
}

# RQ2-specific start-stop and interaction Cox models. These are deliberately
# implemented here so the exact survival::coxph formulas are visible in RQ2.
fit_age_trauma_sex_interaction <- function(data, exposure, model_name, ph_plot_dir = NULL) {
  required <- c("age_entry", "age_exit", "event", "sex", exposure)
  require_columns(data, required, model_name)
  d <- complete_model_data(data, required)
  d <- d[d$event %in% c(0,1) & is.finite(d$age_entry) & is.finite(d$age_exit) & d$age_exit > d$age_entry & is.finite(d[[exposure]]), , drop=FALSE]
  d$sex <- droplevels(factor(as.character(d$sex), levels=c("Male","Female")))
  if (nlevels(d$sex) < 2L) stop(model_name, " requires both sex strata.", call.=FALSE)
  center_value <- mean(d[[exposure]], na.rm=TRUE)
  d$age_trauma_c <- d[[exposure]] - center_value
  d$age_trauma_c_x_female <- d$age_trauma_c * as.numeric(d$sex == "Female")
  form <- survival::Surv(age_entry, age_exit, event) ~ age_trauma_c + age_trauma_c_x_female + strata(sex)
  progress_log("Cox model", paste0(model_name, "; survival::coxph; formula: ", paste(deparse(form), collapse=" ")))
  fit <- survival::coxph(form, data=d, ties="efron", x=TRUE, y=TRUE, model=TRUE, singular.ok=TRUE)
  attr(fit,"model_name") <- model_name; attr(fit,"analysis_data") <- d; attr(fit,"center_value") <- center_value
  ph <- cox_ph_table(fit, model_name, main_exposures=c("age_trauma_c","age_trauma_c_x_female")); ph$Model_phase <- "Exploratory sex interaction"
  ph_plot_audit <- if (!is.null(ph_plot_dir) && nzchar(ph_plot_dir)) save_cox_zph_plots(fit, ph_plot_dir, model_name, "Exploratory sex interaction", ph_table=ph) else data.frame()
  effects <- extract_pooled_sex_interaction(fit, "age_trauma_c", "age_trauma_c_x_female", exposure, "HR per 1-year older age at trauma", center_value)
  audit <- data.frame(Model=model_name, Exposure=exposure, Center_value=center_value, N=fit$n, Events=fit$nevent,
    Male_N=sum(d$sex=="Male"), Male_events=sum(d$event[d$sex=="Male"]), Female_N=sum(d$sex=="Female"), Female_events=sum(d$event[d$sex=="Female"]),
    Exposure_min=min(d[[exposure]],na.rm=TRUE), Exposure_max=max(d[[exposure]],na.rm=TRUE), Formula=paste(deparse(form),collapse=" "), Exploratory=TRUE, stringsAsFactors=FALSE)
  list(fit=fit, terms=extract_cox_terms(fit,model_name), stats=extract_model_stats(fit,model_name), ph=ph,
       ph_plot_audit=ph_plot_audit, sex_specific=effects$sex_specific, interaction_test=effects$interaction_test, audit=audit)
}

build_cutoff_long <- function(data, trauma_age_col, exit_col="age_exit", cutoff_months=6L) {
  required <- c("record_id","sex","event","age_entry",exit_col,trauma_age_col); require_columns(data, required, "RQ2 cutoff model")
  d <- complete_model_data(data, required); d <- d[d[[exit_col]] > d$age_entry & d$event %in% c(0,1), , drop=FALSE]
  rows <- vector("list", nrow(d)*2L); k <- 0L
  for (i in seq_len(nrow(d))) {
    start <- d$age_entry[[i]]; stop <- d[[exit_col]][[i]]; cut_age <- d[[trauma_age_col]][[i]] + cutoff_months/12
    if (start < cut_age) { k <- k+1L; first_stop <- min(stop,cut_age); rows[[k]] <- data.frame(record_id=d$record_id[[i]],sex=as.character(d$sex[[i]]),start=start,stop=first_stop,event=as.integer(d$event[[i]]==1 && stop<=cut_age),within_cutoff=1,age_at_trauma=d[[trauma_age_col]][[i]]) }
    if (stop > cut_age) { k <- k+1L; rows[[k]] <- data.frame(record_id=d$record_id[[i]],sex=as.character(d$sex[[i]]),start=max(start,cut_age),stop=stop,event=as.integer(d$event[[i]]==1),within_cutoff=0,age_at_trauma=d[[trauma_age_col]][[i]]) }
  }
  if (!k) return(data.frame())
  out <- do.call(rbind, rows[seq_len(k)]); out <- out[out$stop > out$start,,drop=FALSE]; out$sex <- factor(out$sex,levels=c("Male","Female")); out
}

fit_cutoff_cox <- function(long_data, model_name, sex_mode=c("strata","covariate")) {
  sex_mode <- match.arg(sex_mode)
  form <- if (sex_mode == "strata") {
    survival::Surv(start, stop, event) ~ within_cutoff + age_at_trauma + strata(sex) + cluster(record_id)
  } else {
    survival::Surv(start, stop, event) ~ within_cutoff + age_at_trauma + sex + cluster(record_id)
  }
  progress_log("Cox model", paste0(model_name, "; survival::coxph start-stop; formula: ", paste(deparse(form),collapse=" ")))
  fit <- survival::coxph(form, data=long_data, ties="efron", x=TRUE, y=TRUE, model=TRUE, singular.ok=TRUE)
  attr(fit,"model_name") <- model_name; attr(fit,"analysis_data") <- long_data; attr(fit,"sex_handling") <- sex_mode; fit
}

fit_cutoff_sex_interaction <- function(long_data, model_name) {
  d <- complete_model_data(long_data, c("record_id","sex","start","stop","event","within_cutoff","age_at_trauma"))
  d <- d[d$event %in% c(0,1) & is.finite(d$start) & is.finite(d$stop) & d$stop>d$start & d$within_cutoff %in% c(0,1),,drop=FALSE]
  d$sex <- droplevels(factor(as.character(d$sex),levels=c("Male","Female"))); d$within_cutoff_x_female <- d$within_cutoff * as.numeric(d$sex=="Female")
  form <- survival::Surv(start, stop, event) ~ within_cutoff + within_cutoff_x_female + age_at_trauma + strata(sex) + cluster(record_id)
  progress_log("Cox model", paste0(model_name, "; survival::coxph interaction; formula: ", paste(deparse(form),collapse=" ")))
  fit <- survival::coxph(form,data=d,ties="efron",x=TRUE,y=TRUE,model=TRUE,singular.ok=TRUE)
  attr(fit,"model_name") <- model_name; attr(fit,"analysis_data") <- d
  effects <- extract_pooled_sex_interaction(fit,"within_cutoff","within_cutoff_x_female","within versus after the prespecified cutoff","Within-cutoff versus after-cutoff hazard ratio",NA_real_)
  list(fit=fit,terms=extract_cox_terms(fit,model_name),stats=extract_model_stats(fit,model_name),sex_specific=effects$sex_specific,interaction_test=effects$interaction_test)
}

progress_log("RQ2", sprintf("helper API %d loaded from %s", SURVIVAL_HELPER_API, normalizePath(helper_path, mustWork = FALSE)))

INVENTORY <- c(
  "Latency descriptives from three time origins",
  "Midpoint-onset interval-censoring sensitivity",
  "RMST with participant-bootstrap confidence intervals",
  "Delayed-entry Cox model for female versus male sex",
  "One-time sex PH assessment with consistent subsequent sex handling",
  "cox.zph proportional-hazards diagnostics with scaled Schoenfeld residual plots, including sex",
  "Age-varying coefficient models for PH-flagged substantive predictors",
  "Start-stop proximal-risk models for 3-12-month cutoffs",
  "Midpoint-onset proximal-risk sensitivity models",
  "Exploratory pooled age-at-trauma-by-sex interaction models for most recent and worst trauma",
  "Prespecified 6-month proximal-risk-by-sex interaction sensitivity models",
  "Time-updated new lifetime-worst-trauma model",
  "Figure-ready cutoff-HR table for Python/Matplotlib rendering"
)
analysis_inventory("RQ2: Time to onset", INVENTORY)

out_dir <- file.path(OUTPUT_ROOT, "RQ2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ph_plot_dir <- file.path(out_dir, "PH_residual_plots")
dir.create(ph_plot_dir, recursive = TRUE, showWarnings = FALSE)
RMST_TAU_OVERRIDE <- trimws(Sys.getenv("RQ2_RMST_TAU", ""))
RMST_TAU_QUANTILE <- as.numeric(Sys.getenv("RQ2_RMST_TAU_QUANTILE", "0.95"))
SEX_INTERACTION_CUTOFF_MONTHS <- as.integer(Sys.getenv("RQ2_SEX_INTERACTION_CUTOFF_MONTHS", "6"))
if (!is.finite(SEX_INTERACTION_CUTOFF_MONTHS) || SEX_INTERACTION_CUTOFF_MONTHS < 1L || SEX_INTERACTION_CUTOFF_MONTHS > 24L) {
  stop("RQ2_SEX_INTERACTION_CUTOFF_MONTHS must be an integer from 1 through 24.", call. = FALSE)
}
if (!is.finite(RMST_TAU_QUANTILE) || RMST_TAU_QUANTILE <= 0 || RMST_TAU_QUANTILE > 1) {
  stop("RQ2_RMST_TAU_QUANTILE must be in (0, 1].", call. = FALSE)
}

progress_log("RQ2", "loading and validating the prepared dataset")
data <- read_survival_data()
progress_log("RQ2", sprintf("loaded N=%d participants; events=%d", nrow(data), sum(data$event == 1, na.rm = TRUE)))
require_columns(
  data,
  c(
    "record_id", "sex", "sex_female", "event", "age_entry", "age_exit",
    "time_since_start", "time_since_most_recent_trauma", "time_since_worst_trauma",
    "months_since_most_recent_at_baseline", "months_since_worst_at_baseline",
    "age_at_most_recent_trauma", "age_at_worst_trauma",
    "first_positive_months", "last_negative_before_first_positive_months"
  ),
  "RQ2"
)
midpoint_data <- build_midpoint_dataset(data)
midpoint_audit <- attr(midpoint_data, "midpoint_audit")
# -----------------------------------------------------------------------------
# Latency descriptives
# -----------------------------------------------------------------------------
progress_log("RQ2", "calculating latency descriptives")
origin_duration_cols <- c(
  "Study entry" = "time_since_start",
  "Most recent trauma" = "time_since_most_recent_trauma",
  "Worst trauma" = "time_since_worst_trauma"
)
latency_rows <- list()
for (origin in names(origin_duration_cols)) {
  col <- origin_duration_cols[[origin]]
  primary <- data[[col]][data$event == 1] / 12
  midpoint <- midpoint_data[[paste0(col, "_midpoint")]][midpoint_data$event == 1] / 12
  latency_rows[[length(latency_rows) + 1L]] <- latency_descriptives(primary, "First detected", origin)
  latency_rows[[length(latency_rows) + 1L]] <- latency_descriptives(midpoint, "Midpoint sensitivity", origin)
}
latency <- rbind_fill(latency_rows)

# -----------------------------------------------------------------------------
# RMST from primary and midpoint survival times
# -----------------------------------------------------------------------------
make_scenario_origin <- function(d, origin, scenario) {
  if (scenario == "First detected") return(make_origin_data(d, origin))
  if (origin == "Study entry") {
    out <- data.frame(record_id = d$record_id, sex = d$sex, event = d$event,
                      entry = 0, duration = d$time_since_start_midpoint / 12)
  } else if (origin == "Most recent trauma") {
    out <- data.frame(record_id = d$record_id, sex = d$sex, event = d$event,
                      entry = d$months_since_most_recent_at_baseline / 12,
                      duration = d$time_since_most_recent_trauma_midpoint / 12)
  } else {
    out <- data.frame(record_id = d$record_id, sex = d$sex, event = d$event,
                      entry = d$months_since_worst_at_baseline / 12,
                      duration = d$time_since_worst_trauma_midpoint / 12)
  }
  out <- out[complete.cases(out[, c("event", "entry", "duration")]), , drop = FALSE]
  out[out$event %in% c(0, 1) & out$entry >= 0 & out$duration > out$entry, , drop = FALSE]
}

# Select one common tau per time origin. By default, tau is the 95th percentile
# of the primary first-detected duration distribution, capped at common observed
# support across the primary/midpoint and sex-specific analysis cells. A numeric
# RQ2_RMST_TAU environment variable overrides the percentile candidate.
rmst_tau_rows <- list()
rmst_tau_by_origin <- setNames(rep(NA_real_, length(origin_duration_cols)), names(origin_duration_cols))
progress_log("RQ2", "selecting RMST tau for each time origin")
for (origin in names(origin_duration_cols)) {
  primary_overall <- make_scenario_origin(data, origin, "First detected")
  primary_duration <- primary_overall$duration[is.finite(primary_overall$duration)]
  override <- suppressWarnings(as.numeric(RMST_TAU_OVERRIDE))
  candidate <- if (nzchar(RMST_TAU_OVERRIDE) && is.finite(override) && override > 0) {
    override
  } else {
    unname(quantile(primary_duration, probs = RMST_TAU_QUANTILE, na.rm = TRUE, type = 7))
  }

  support_values <- numeric()
  for (scenario in c("First detected", "Midpoint sensitivity")) {
    source_data <- if (scenario == "First detected") data else midpoint_data
    scenario_data <- make_scenario_origin(source_data, origin, scenario)
    for (group_name in c("Overall", "Female", "Male")) {
      g <- if (group_name == "Overall") scenario_data else scenario_data[scenario_data$sex == group_name, , drop = FALSE]
      if (nrow(g) && any(is.finite(g$duration))) support_values <- c(support_values, max(g$duration, na.rm = TRUE))
    }
  }
  common_support_cap <- if (length(support_values)) min(support_values) else NA_real_
  tau <- min(candidate, common_support_cap, na.rm = TRUE)
  if (!is.finite(tau) || tau <= 0) tau <- NA_real_
  rmst_tau_by_origin[[origin]] <- tau
  rmst_tau_rows[[length(rmst_tau_rows) + 1L]] <- data.frame(
    Time_origin = origin,
    Tau_method = if (nzchar(RMST_TAU_OVERRIDE) && is.finite(override) && override > 0) "User-specified override" else paste0("Primary duration quantile: ", RMST_TAU_QUANTILE),
    Quantile_probability = if (nzchar(RMST_TAU_OVERRIDE) && is.finite(override) && override > 0) NA_real_ else RMST_TAU_QUANTILE,
    Candidate_tau_years = candidate,
    Common_support_cap_years = common_support_cap,
    Final_tau_years = tau,
    Primary_analysis_N = nrow(primary_overall),
    Primary_events = sum(primary_overall$event),
    stringsAsFactors = FALSE
  )
}
rmst_tau_selection <- rbind_fill(rmst_tau_rows)

rmst_rows <- list()
rmst_apa <- character()
for (scenario_index in seq_along(c("First detected", "Midpoint sensitivity"))) {
  scenario <- c("First detected", "Midpoint sensitivity")[[scenario_index]]
  d_source <- if (scenario == "First detected") data else midpoint_data
  for (origin_index in seq_along(names(origin_duration_cols))) {
    origin <- names(origin_duration_cols)[[origin_index]]
    origin_data <- make_scenario_origin(d_source, origin, scenario)
    tau <- rmst_tau_by_origin[[origin]]
    for (group_name in c("Overall", "Female", "Male")) {
      d <- if (group_name == "Overall") origin_data else origin_data[origin_data$sex == group_name, , drop = FALSE]
      if (!nrow(d) || !sum(d$event) || !is.finite(tau)) next
      fit <- km_fit(d)
      rmst <- rmst_from_fit(fit, tau)
      boot <- bootstrap_km(
        d,
        times = numeric(0),
        tau = tau,
        B = BOOT_B,
        seed = RNG_SEED + scenario_index * 10000L + origin_index * 1000L + match(group_name, c("Overall", "Female", "Male")),
        label = paste("RQ2 RMST", scenario, origin, group_name, sep = " | ")
      )
      ci <- percentile_ci(boot$rmst)
      rmst_rows[[length(rmst_rows) + 1L]] <- data.frame(
        Scenario = scenario,
        Time_origin = origin,
        Group = group_name,
        Tau_years = tau,
        Tau_reestimated_within_bootstrap = FALSE,
        RMST_years = rmst,
        Bootstrap_CI_lower = ci[[1L]],
        Bootstrap_CI_upper = ci[[2L]],
        Analysis_N = nrow(d),
        Events = sum(d$event),
        Bootstrap_replicates_requested = BOOT_B,
        Bootstrap_replicates_valid = sum(is.finite(boot$rmst)),
        stringsAsFactors = FALSE
      )
      rmst_apa <- c(rmst_apa, sprintf(
        "%s, %s, %s: RMST through tau = %s years was %s years, bootstrap 95%% CI [%s, %s] (N = %s; events = %s).",
        scenario, origin, group_name, fmt_num(tau), fmt_num(rmst), fmt_num(ci[[1L]]), fmt_num(ci[[2L]]), nrow(d), sum(d$event)
      ))
    }
  }
}
rmst_results <- rbind_fill(rmst_rows)

# -----------------------------------------------------------------------------
# Fixed-coefficient Cox models and cox.zph diagnostics
# -----------------------------------------------------------------------------
progress_log("RQ2", "fitting fixed-coefficient Cox models and PH diagnostics")
term_tables <- list()
model_tables <- list()
ph_tables <- list()
extended_rows <- list()
extended_terms <- list()
model_errors <- list()
sex_decision_tables <- list()
ph_plot_audits <- list()

run_fixed_model <- function(
  model_name,
  covariates,
  strata = character(),
  main_exposures = covariates,
  nuisance = character()
) {
  tryCatch({
    fit <- fit_age_scale_cox(data, covariates, strata = strata, model_name = model_name)
    ph <- cox_ph_table(fit, model_name, main_exposures = main_exposures, nuisance_terms = nuisance)
    ph_plot_audit <- save_cox_zph_plots(
      fit, ph_plot_dir, model_name = model_name, model_phase = "Final model", ph_table = ph
    )
    ph_plot_audits[[length(ph_plot_audits) + 1L]] <<- ph_plot_audit
    obj <- list(
      fit = fit,
      ph = ph,
      final_ph = ph,
      sex_stratified = identical(strata, "sex"),
      sex_as_covariate = "sex" %in% covariates
    )
    terms <- extract_cox_terms(fit, model_name)
    stats <- extract_model_stats(fit, model_name)
    term_tables[[length(term_tables) + 1L]] <<- terms
    model_tables[[length(model_tables) + 1L]] <<- stats
    ph_tables[[length(ph_tables) + 1L]] <<- ph
    obj$terms <- terms
    obj$stats <- stats
    obj
  }, error = function(e) {
    model_errors[[length(model_errors) + 1L]] <<- data.frame(Model = model_name, Error = conditionMessage(e), stringsAsFactors = FALSE)
    NULL
  })
}

# First Cox model: assess sex PH once on the full analytic cohort.
sex_model <- run_fixed_model(
  "RQ2 sex model",
  covariates = "sex",
  main_exposures = "sex",
  nuisance = "sex"
)
if (is.null(sex_model)) stop("RQ2 initial sex PH assessment could not be fitted.", call. = FALSE)
SEX_PH_P <- ph_term_p(sex_model$ph, "sex")
SEX_STRATIFY <- is.finite(SEX_PH_P) && SEX_PH_P < PH_ALPHA
GLOBAL_SEX_MODE <- if (SEX_STRATIFY) "strata" else "covariate"
GLOBAL_SEX_STRATA <- if (SEX_STRATIFY) "sex" else character()
GLOBAL_SEX_COVARIATES <- if (SEX_STRATIFY) character() else "sex"
sex_row <- sex_model$ph[sex_model$ph$PH_term == "sex", , drop = FALSE]
if (!nrow(sex_row)) sex_row <- sex_model$ph[grep("^sex", sex_model$ph$PH_term), , drop = FALSE]
sex_effect <- sex_effect_from_fit(sex_model$fit)
sex_decision_tables[[1L]] <- data.frame(
  Model = "RQ2 initial sex PH assessment",
  Sex_PH_tested = is.finite(SEX_PH_P),
  Sex_PH_chisq = if (nrow(sex_row)) suppressWarnings(as.numeric(sex_row$PH_chisq[[1L]])) else NA_real_,
  Sex_PH_df = if (nrow(sex_row)) suppressWarnings(as.numeric(sex_row$PH_df[[1L]])) else NA_real_,
  Sex_PH_p = SEX_PH_P,
  PH_alpha = PH_ALPHA,
  Sex_PH_concern = if (is.finite(SEX_PH_P)) SEX_STRATIFY else NA,
  Final_sex_handling = if (SEX_STRATIFY) "strata(sex) used consistently in all subsequent Cox models" else "sex retained as a covariate consistently in all subsequent Cox models",
  Sex_stratified = SEX_STRATIFY,
  Sex_retained_as_covariate = !SEX_STRATIFY,
  Screening_sex_HR = scalar_numeric_or_na(sex_effect$HR),
  Screening_sex_CI_lower = scalar_numeric_or_na(sex_effect$CI_lower),
  Screening_sex_CI_upper = scalar_numeric_or_na(sex_effect$CI_upper),
  Screening_sex_coefficient_p = scalar_numeric_or_na(sex_effect$coefficient_p),
  Screening_N = sex_model$fit$n,
  Screening_events = sex_model$fit$nevent,
  Final_formula = if (SEX_STRATIFY) "Subsequent models: ... + strata(sex)" else "Subsequent models: ... + sex",
  stringsAsFactors = FALSE
)
global_sex_decision <- sex_decision_tables[[1L]]
global_sex_decision$Decision_owner <- "RQ2"
global_sex_decision$Decision_scope <- "RQ2, RQ3, and prespecified RQ3 sensitivity Cox models"
global_sex_decision$Sex_mode <- GLOBAL_SEX_MODE
global_sex_decision$Decision_file_purpose <- "Read by RQ3 and sensitivity scripts so sex is not re-tested outside RQ2"
global_sex_decision_path <- file.path(out_dir, "RQ2_global_sex_handling_decision.csv")
write_csv_safe(global_sex_decision, global_sex_decision_path)
progress_log("RQ2 sex handling", paste0("one-time sex PH p ", fmt_p(SEX_PH_P), "; subsequent models use ", if (SEX_STRATIFY) "strata(sex)" else "sex as a covariate", "; decision saved to ", normalizePath(global_sex_decision_path, mustWork = FALSE)))

age_recent_model <- run_fixed_model(
  "RQ2 age at most recent trauma",
  covariates = c("age_at_most_recent_trauma", GLOBAL_SEX_COVARIATES),
  strata = GLOBAL_SEX_STRATA,
  main_exposures = "age_at_most_recent_trauma",
  nuisance = GLOBAL_SEX_COVARIATES
)

age_worst_model <- run_fixed_model(
  "RQ2 age at worst trauma",
  covariates = c("age_at_worst_trauma", GLOBAL_SEX_COVARIATES),
  strata = GLOBAL_SEX_STRATA,
  main_exposures = "age_at_worst_trauma",
  nuisance = GLOBAL_SEX_COVARIATES
)

# Retain primary Cox estimates; add extended models as supplementary results when flagged.
if (!is.null(sex_model) && any(sex_model$ph$PH_term == "sex" & sex_model$ph$PH_concern %in% TRUE)) {
  ext <- tryCatch(
    fit_extended_age_effect(data, "sex_female", strata = character(), model_name = "RQ2 age-varying female versus male effect"),
    error = function(e) { model_errors[[length(model_errors) + 1L]] <<- data.frame(Model = "RQ2 age-varying female versus male effect", Error = conditionMessage(e)); NULL }
  )
  if (!is.null(ext)) {
    extended_rows[[length(extended_rows) + 1L]] <- ext$age_specific
    extended_terms[[length(extended_terms) + 1L]] <- ext$interaction
    model_tables[[length(model_tables) + 1L]] <- extract_model_stats(ext$fit, "RQ2 age-varying female versus male effect")
  }
}

for (spec in list(
  list(object = age_recent_model, exposure = "age_at_most_recent_trauma", name = "RQ2 age-varying effect: age at most recent trauma"),
  list(object = age_worst_model, exposure = "age_at_worst_trauma", name = "RQ2 age-varying effect: age at worst trauma")
)) {
  final_ph <- if (!is.null(spec$object)) spec$object$final_ph %||% spec$object$ph else data.frame()
  if (!is.null(spec$object) && any(final_ph$PH_term == spec$exposure & final_ph$PH_concern %in% TRUE)) {
    ext <- tryCatch(
      fit_extended_age_effect(
        data,
        spec$exposure,
        adjustment_terms = GLOBAL_SEX_COVARIATES,
        strata = GLOBAL_SEX_STRATA,
        model_name = spec$name
      ),
      error = function(e) { model_errors[[length(model_errors) + 1L]] <<- data.frame(Model = spec$name, Error = conditionMessage(e)); NULL }
    )
    if (!is.null(ext)) {
      extended_rows[[length(extended_rows) + 1L]] <- ext$age_specific
      extended_terms[[length(extended_terms) + 1L]] <- ext$interaction
      model_tables[[length(model_tables) + 1L]] <- extract_model_stats(ext$fit, spec$name)
    }
  }
}


# -----------------------------------------------------------------------------
# Exploratory pooled sex-effect-modification models for age at trauma
# -----------------------------------------------------------------------------
progress_log("RQ2", "fitting exploratory age-at-trauma-by-sex interaction models")
age_sex_interaction_terms <- list()
age_sex_interaction_stats <- list()
age_sex_interaction_ph <- list()
age_sex_specific_rows <- list()
age_sex_interaction_tests <- list()
age_sex_interaction_audits <- list()
age_sex_interaction_errors <- list()

for (spec in list(
  list(exposure = "age_at_most_recent_trauma", label = "most recent trauma"),
  list(exposure = "age_at_worst_trauma", label = "worst trauma")
)) {
  model_name <- paste0("RQ2 exploratory sex interaction: age at ", spec$label)
  obj <- tryCatch(
    fit_age_trauma_sex_interaction(
      data,
      exposure = spec$exposure,
      model_name = model_name,
      ph_plot_dir = ph_plot_dir
    ),
    error = function(e) {
      age_sex_interaction_errors[[length(age_sex_interaction_errors) + 1L]] <<-
        data.frame(Model = model_name, Error = conditionMessage(e), stringsAsFactors = FALSE)
      NULL
    }
  )
  if (!is.null(obj)) {
    age_sex_interaction_terms[[length(age_sex_interaction_terms) + 1L]] <- obj$terms
    age_sex_interaction_stats[[length(age_sex_interaction_stats) + 1L]] <- obj$stats
    age_sex_interaction_ph[[length(age_sex_interaction_ph) + 1L]] <- obj$ph
    obj$sex_specific$Exploratory <- TRUE
    obj$interaction_test$Exploratory <- TRUE
    age_sex_specific_rows[[length(age_sex_specific_rows) + 1L]] <- obj$sex_specific
    age_sex_interaction_tests[[length(age_sex_interaction_tests) + 1L]] <- obj$interaction_test
    age_sex_interaction_audits[[length(age_sex_interaction_audits) + 1L]] <- obj$audit
    ph_plot_audits[[length(ph_plot_audits) + 1L]] <- obj$ph_plot_audit
  }
}

age_sex_interaction_term_results <- rbind_fill(age_sex_interaction_terms)
age_sex_interaction_model_results <- rbind_fill(age_sex_interaction_stats)
age_sex_interaction_ph_results <- rbind_fill(age_sex_interaction_ph)
age_sex_specific_results <- rbind_fill(age_sex_specific_rows)
age_sex_interaction_test_results <- rbind_fill(age_sex_interaction_tests)
age_sex_interaction_audit_results <- rbind_fill(age_sex_interaction_audits)
if (nrow(age_sex_interaction_test_results)) {
  age_sex_interaction_test_results$p_FDR_BH <- stats::p.adjust(age_sex_interaction_test_results$p, method = "BH")
}

# -----------------------------------------------------------------------------
# Proximal-risk start-stop models for 3-12-month cutoffs
# -----------------------------------------------------------------------------
cutoff_terms <- list()
cutoff_stats <- list()
cutoff_audits <- list()
cutoff_errors <- list()

# The explicitly time-updated cutoff models do not receive their own standard
# cox.zph diagnostic. They inherit the single sex-handling decision above.
cutoff_sex_mode <- GLOBAL_SEX_MODE

cutoff_origins <- list(
  "Most recent trauma" = "age_at_most_recent_trauma",
  "Worst trauma" = "age_at_worst_trauma"
)

for (scenario in c("First detected", "Midpoint sensitivity")) {
  d_source <- if (scenario == "First detected") data else midpoint_data
  exit_col <- if (scenario == "First detected") "age_exit" else "age_exit_midpoint"
  for (origin in names(cutoff_origins)) {
    progress_log("RQ2 cutoff models", paste(scenario, origin, "cutoffs 3-12 months", sep = " | "))
    trauma_col <- cutoff_origins[[origin]]
    for (cutoff in 3:12) {
      model_name <- sprintf("RQ2 %s, %s, %d-month cutoff", scenario, origin, cutoff)
      tryCatch({
        long <- build_cutoff_long(d_source, trauma_col, exit_col = exit_col, cutoff_months = cutoff)
        cut_age <- d_source[[trauma_col]] + cutoff / 12
        
        eligible <- (
          is.finite(d_source$age_entry) &
            is.finite(d_source[[exit_col]]) &
            is.finite(cut_age)
        )
        
        n_at_risk_cutoff <- sum(
          eligible &
            d_source$age_entry < cut_age &
            d_source[[exit_col]] >= cut_age
        )
        
        fit <- fit_cutoff_cox(long, model_name, sex_mode = cutoff_sex_mode)
        terms <- extract_cox_terms(fit, model_name)
        terms$Scenario <- scenario
        terms$Time_origin <- origin
        terms$Cutoff_months <- cutoff
        stats <- extract_model_stats(fit, model_name)
        stats$Scenario <- scenario
        stats$Time_origin <- origin
        stats$Cutoff_months <- cutoff
        cutoff_terms[[length(cutoff_terms) + 1L]] <- terms
        cutoff_stats[[length(cutoff_stats) + 1L]] <- stats
        cutoff_audits[[length(cutoff_audits) + 1L]] <- data.frame(
          Model = model_name,
          Scenario = scenario,
          Time_origin = origin,
          Cutoff_months = cutoff,
          
          Participants = length(unique(long$record_id)),
          N_at_risk_cutoff = n_at_risk_cutoff,
          
          Interval_rows = nrow(long),
          Events = sum(long$event),
          Person_years_within = sum(
            (long$stop - long$start)[long$within_cutoff == 1]
          ),
          Person_years_after = sum(
            (long$stop - long$start)[long$within_cutoff == 0]
          ),
          Events_within = sum(
            long$event[long$within_cutoff == 1]
          ),
          Events_after = sum(
            long$event[long$within_cutoff == 0]
          ),
          
          PH_test_applied = FALSE,
          Sex_handling = if (cutoff_sex_mode == "strata")
            "strata(sex)" else "sex covariate",
          Sex_decision_source = "One-time RQ2 sex PH assessment",
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        cutoff_errors[[length(cutoff_errors) + 1L]] <- data.frame(Model = model_name, Error = conditionMessage(e), stringsAsFactors = FALSE)
      })
    }
  
  }
}
cutoff_audit_results <- rbind_fill(cutoff_audits)

figure2_counts <- cutoff_audit_results[
  cutoff_audit_results$Scenario == "First detected",
  c(
    "Time_origin",
    "Cutoff_months",
    "N_at_risk_cutoff",
    "Events_within"
  ),
  drop = FALSE
]

write_csv_safe(
  figure2_counts,
  file.path(out_dir, "Figure2_post_trauma_counts.csv")
)

# -----------------------------------------------------------------------------
# Exploratory sex interaction at one prespecified proximal-risk cutoff
# -----------------------------------------------------------------------------
progress_log(
  "RQ2",
  paste0("fitting exploratory proximal-risk-by-sex interaction models at ", SEX_INTERACTION_CUTOFF_MONTHS, " months")
)
proximal_sex_interaction_terms <- list()
proximal_sex_interaction_stats <- list()
proximal_sex_specific_rows <- list()
proximal_sex_interaction_tests <- list()
proximal_sex_interaction_audits <- list()
proximal_sex_interaction_errors <- list()

for (scenario in c("First detected", "Midpoint sensitivity")) {
  d_source <- if (scenario == "First detected") data else midpoint_data
  exit_col <- if (scenario == "First detected") "age_exit" else "age_exit_midpoint"
  for (origin in names(cutoff_origins)) {
    trauma_col <- cutoff_origins[[origin]]
    model_name <- sprintf(
      "RQ2 exploratory sex interaction: %s, %s, %d-month cutoff",
      scenario, origin, SEX_INTERACTION_CUTOFF_MONTHS
    )
    obj <- tryCatch({
      long <- build_cutoff_long(
        d_source,
        trauma_col,
        exit_col = exit_col,
        cutoff_months = SEX_INTERACTION_CUTOFF_MONTHS
      )
      fit_cutoff_sex_interaction(long, model_name)
    }, error = function(e) {
      proximal_sex_interaction_errors[[length(proximal_sex_interaction_errors) + 1L]] <<-
        data.frame(Model = model_name, Error = conditionMessage(e), stringsAsFactors = FALSE)
      NULL
    })
    if (!is.null(obj)) {
      terms <- obj$terms
      terms$Scenario <- scenario
      terms$Time_origin <- origin
      terms$Cutoff_months <- SEX_INTERACTION_CUTOFF_MONTHS
      stats <- obj$stats
      stats$Scenario <- scenario
      stats$Time_origin <- origin
      stats$Cutoff_months <- SEX_INTERACTION_CUTOFF_MONTHS
      sex_specific <- obj$sex_specific
      sex_specific$Scenario <- scenario
      sex_specific$Time_origin <- origin
      sex_specific$Cutoff_months <- SEX_INTERACTION_CUTOFF_MONTHS
      sex_specific$Exploratory <- TRUE
      test <- obj$interaction_test
      test$Scenario <- scenario
      test$Time_origin <- origin
      test$Cutoff_months <- SEX_INTERACTION_CUTOFF_MONTHS
      test$Exploratory <- TRUE
      fitted_data <- attr(obj$fit, "analysis_data")
      audit_rows <- lapply(c("Male", "Female"), function(sex_label) {
        x <- fitted_data[fitted_data$sex == sex_label, , drop = FALSE]
        data.frame(
          Model = model_name,
          Scenario = scenario,
          Time_origin = origin,
          Cutoff_months = SEX_INTERACTION_CUTOFF_MONTHS,
          Sex = sex_label,
          Participants = length(unique(x$record_id)),
          Interval_rows = nrow(x),
          Events = sum(x$event),
          Person_years_within = sum((x$stop - x$start)[x$within_cutoff == 1]),
          Person_years_after = sum((x$stop - x$start)[x$within_cutoff == 0]),
          Events_within = sum(x$event[x$within_cutoff == 1]),
          Events_after = sum(x$event[x$within_cutoff == 0]),
          Sex_handling = "strata(sex) with within-cutoff-by-female interaction",
          PH_test_applied = FALSE,
          Exploratory = TRUE,
          stringsAsFactors = FALSE
        )
      })
      proximal_sex_interaction_terms[[length(proximal_sex_interaction_terms) + 1L]] <- terms
      proximal_sex_interaction_stats[[length(proximal_sex_interaction_stats) + 1L]] <- stats
      proximal_sex_specific_rows[[length(proximal_sex_specific_rows) + 1L]] <- sex_specific
      proximal_sex_interaction_tests[[length(proximal_sex_interaction_tests) + 1L]] <- test
      proximal_sex_interaction_audits[[length(proximal_sex_interaction_audits) + 1L]] <- do.call(rbind, audit_rows)
    }
  }
}

proximal_sex_interaction_term_results <- rbind_fill(proximal_sex_interaction_terms)
proximal_sex_interaction_model_results <- rbind_fill(proximal_sex_interaction_stats)
proximal_sex_specific_results <- rbind_fill(proximal_sex_specific_rows)
proximal_sex_interaction_test_results <- rbind_fill(proximal_sex_interaction_tests)
proximal_sex_interaction_audit_results <- rbind_fill(proximal_sex_interaction_audits)
if (nrow(proximal_sex_interaction_test_results)) {
  proximal_sex_interaction_test_results$p_FDR_BH <- stats::p.adjust(proximal_sex_interaction_test_results$p, method = "BH")
}



# -----------------------------------------------------------------------------
# New lifetime-worst-trauma analyses
# -----------------------------------------------------------------------------
# 1. First report the simple descriptive/binary new-worst-trauma indicator.
# 2. Then fit the time-updated model without trauma_type_count, to avoid
#    immortal-time bias while preserving the specific new-worst-trauma estimand.
# -----------------------------------------------------------------------------

new_worst_terms <- list()
new_worst_stats <- list()
new_worst_audit <- data.frame()
new_worst_binary_summary <- data.frame()
new_worst_binary_terms <- list()
new_worst_binary_stats <- list()
new_worst_errors <- list()

if (all(c("new_worst_trauma_flag", "new_worst_trauma_months_since_start") %in% names(data))) {
  
  # ---------------------------------------------------------------------------
  # A. Descriptive binary summary: ever had a new lifetime-worst trauma after study start
  # ---------------------------------------------------------------------------
  binary_flag <- data$new_worst_trauma_flag
  binary_flag <- ifelse(binary_flag %in% c(0, 1), binary_flag, NA_real_)
  
  new_worst_binary_summary <- data.frame(
    Measure = "Ever reported a new lifetime-worst trauma after study start",
    N_with_nonmissing_flag = sum(!is.na(binary_flag)),
    N_new_worst_trauma = sum(binary_flag == 1, na.rm = TRUE),
    Percent_new_worst_trauma = if (sum(!is.na(binary_flag)) > 0L) {
      100 * sum(binary_flag == 1, na.rm = TRUE) / sum(!is.na(binary_flag))
    } else {
      NA_real_
    },
    N_without_new_worst_trauma = sum(binary_flag == 0, na.rm = TRUE),
    Missing_flag = sum(is.na(binary_flag)),
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # B. Simple binary Cox model
  # This treats new lifetime-worst trauma as ever/never and does not preserve timing.
  # Therefore it is descriptive/supplementary. The time-updated model below is
  # the primary inferential model for this exposure.
  # ---------------------------------------------------------------------------
  binary_model_name <- "RQ2 binary ever-new lifetime-worst trauma"
  
  tryCatch({
    d_binary <- data
    d_binary$new_worst_trauma_binary <- binary_flag
    
    binary_covars <- "new_worst_trauma_binary"
    
    if (cutoff_sex_mode == "covariate") {
      binary_covars <- c(binary_covars, "sex")
    }
    
    binary_fit <- fit_age_scale_cox(
      d_binary,
      covariates = binary_covars,
      strata = if (cutoff_sex_mode == "strata") "sex" else character(),
      model_name = binary_model_name
    )
    
    binary_terms <- extract_cox_terms(binary_fit, binary_model_name)
    binary_terms$Analysis_note <- "Ever/never descriptive Cox model; timing not preserved"
    new_worst_binary_terms[[1L]] <- binary_terms
    
    binary_stats <- extract_model_stats(binary_fit, binary_model_name)
    binary_stats$Analysis_note <- "Ever/never descriptive Cox model; timing not preserved"
    new_worst_binary_stats[[1L]] <- binary_stats
    
  }, error = function(e) {
    new_worst_errors[[length(new_worst_errors) + 1L]] <<- data.frame(
      Model = binary_model_name,
      Error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  
  # ---------------------------------------------------------------------------
  # C. Primary timing-preserving model: time-updated new lifetime-worst trauma
  # IMPORTANT: trauma_type_count is intentionally skipped here.
  # ---------------------------------------------------------------------------
  model_name <- "RQ2 time-updated new lifetime-worst trauma"
  
  tryCatch({
    long <- build_new_worst_long(data)
    
    covars <- "new_worst_trauma_tv"
    
    # Do not add trauma_type_count here.
    # This model tests the timing-preserving effect of new lifetime-worst trauma itself.
    if (cutoff_sex_mode == "covariate") {
      covars <- c(covars, "sex")
    }
    
    form <- make_formula(
      "Surv(start, stop, event)",
      covariates = covars,
      strata = if (cutoff_sex_mode == "strata") "sex" else character(),
      cluster = "record_id"
    )
    
    fit <- survival::coxph(
      form,
      data = long,
      ties = "efron",
      x = TRUE,
      y = TRUE,
      model = TRUE,
      singular.ok = TRUE
    )
    
    attr(fit, "model_name") <- model_name
    attr(fit, "analysis_data") <- long
    
    new_worst_terms[[1L]] <- extract_cox_terms(fit, model_name)
    new_worst_stats[[1L]] <- extract_model_stats(fit, model_name)
    
    new_worst_audit <- attr(long, "audit")
    new_worst_audit$Model <- model_name
    new_worst_audit$PH_test_applied <- FALSE
    new_worst_audit$Sex_handling <- if (cutoff_sex_mode == "strata") {
      "strata(sex)"
    } else {
      "sex covariate"
    }
    new_worst_audit$Sex_decision_source <- "One-time RQ2 sex PH assessment"
    new_worst_audit$Trauma_type_count_adjusted <- FALSE
    new_worst_audit$Analysis_note <- "Time-updated model preserves exposure timing and avoids immortal-time bias"
    
  }, error = function(e) {
    new_worst_errors[[length(new_worst_errors) + 1L]] <<- data.frame(
      Model = model_name,
      Error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  
} else {
  new_worst_errors[[1L]] <- data.frame(
    Model = "RQ2 new lifetime-worst trauma analyses",
    Error = "Prepared new-worst-trauma date/flag columns were absent",
    stringsAsFactors = FALSE
  )
}

fixed_terms <- rbind_fill(c(
  term_tables,
  extended_terms,
  cutoff_terms,
  new_worst_binary_terms,
  new_worst_terms
))

fixed_stats <- rbind_fill(c(
  model_tables,
  cutoff_stats,
  new_worst_binary_stats,
  new_worst_stats
))
ph_results <- rbind_fill(ph_tables)
cutoff_term_results <- rbind_fill(cutoff_terms)
cutoff_model_results <- rbind_fill(cutoff_stats)
extended_results <- rbind_fill(extended_rows)
errors <- rbind_fill(c(model_errors, cutoff_errors, age_sex_interaction_errors, proximal_sex_interaction_errors, new_worst_errors))
sex_decisions <- rbind_fill(sex_decision_tables)
ph_plot_audit <- rbind_fill(ph_plot_audits)

progress_log("RQ2", "saving Figure 2 data for Python/Matplotlib rendering")
figure2_data <- cutoff_term_results[
  cutoff_term_results$Scenario == "First detected" & cutoff_term_results$Term == "within_cutoff" &
    is.finite(cutoff_term_results$HR) & is.finite(cutoff_term_results$CI_lower) & is.finite(cutoff_term_results$CI_upper),
  , drop = FALSE
]
if (nrow(figure2_data)) figure2_data$Significance <- p_stars(figure2_data$p)
write_csv_safe(figure2_data, file.path(out_dir, "Figure2_post_trauma_cutoff_HRs_data.csv"))

rq2_model_specifications <- data.frame(
  Analysis = c(
    "Female versus male fixed-coefficient model",
    "One-time sex PH assessment",
    "Age at trauma fixed-coefficient model with global sex handling",
    "Age-at-trauma by sex effect-modification model",
    "Proximal cutoff model",
    "Proximal cutoff by sex interaction model",
    "Time-updated new lifetime-worst-trauma model",
    "PH diagnostic for fixed-coefficient Cox models"
  ),
  R_package = "survival",
  R_function = c("survival::coxph","survival::coxph","survival::coxph","survival::coxph","survival::coxph","survival::coxph","survival::coxph","survival::cox.zph"),
  Model_type = c(
    "Delayed-entry attained-age Cox PH",
    "Delayed-entry attained-age Cox PH diagnostic for sex",
    "Delayed-entry attained-age Cox with globally selected sex handling",
    "Delayed-entry attained-age pooled Cox interaction",
    "Start-stop Cox with time-updated within-cutoff exposure and clustered robust SE",
    "Start-stop pooled Cox interaction with clustered robust SE",
    "Start-stop Cox with time-updated trauma exposure and clustered robust SE",
    "Scaled Schoenfeld residual PH test"
  ),
  Formula = c(
    "Surv(age_entry, age_exit, event) ~ sex",
    "Surv(age_entry, age_exit, event) ~ sex; cox.zph used once to choose subsequent sex handling",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ age_at_trauma + strata(sex)" else "Surv(age_entry, age_exit, event) ~ age_at_trauma + sex",
    "Surv(age_entry, age_exit, event) ~ age_trauma_c + age_trauma_c_x_female + strata(sex)",
    if (SEX_STRATIFY) "Surv(start, stop, event) ~ within_cutoff + age_at_trauma + strata(sex) + cluster(record_id)" else "Surv(start, stop, event) ~ within_cutoff + age_at_trauma + sex + cluster(record_id)",
    "Surv(start, stop, event) ~ within_cutoff + within_cutoff_x_female + age_at_trauma + strata(sex) + cluster(record_id)",
    if (SEX_STRATIFY) "Surv(start, stop, event) ~ new_worst_trauma_tv + strata(sex) + cluster(record_id)" else "Surv(start, stop, event) ~ new_worst_trauma_tv + sex + cluster(record_id)",
    "cox.zph(fit, transform='rank')"
  ),
  Time_scale = c(rep("Attained age (years)",4), rep("Attained age start-stop intervals (years)",3), "Rank-transformed event time for diagnostic"),
  Ties = c(rep("Efron",7), NA),
  Variance = c(rep("Model-based",4), rep("Participant-clustered robust sandwich",3), NA),
  Figure_renderer = "Python / Matplotlib from saved CSV; no manuscript Figure 2 is drawn in R",
  stringsAsFactors = FALSE
)
write_csv_safe(rq2_model_specifications, file.path(out_dir, "RQ2_model_specifications.csv"))

sex_decision_lines <- character()
if (nrow(sex_decisions)) {
  sex_decision_lines <- apply(sex_decisions, 1, function(x) sprintf(
    "%s: sex PH test chi-square(%s) = %s, p %s; final handling: %s.",
    x[["Model"]], fmt_num(as.numeric(x[["Sex_PH_df"]]), 0),
    fmt_num(as.numeric(x[["Sex_PH_chisq"]])), fmt_p(as.numeric(x[["Sex_PH_p"]])),
    x[["Final_sex_handling"]]
  ))
}

latency_apa <- apply(latency, 1, function(x) sprintf(
  "%s, %s: N events = %s, M = %s years, SD = %s, median = %s, mode = %s, range [%s, %s].",
  x[["Scenario"]], x[["Time_origin"]], x[["N_events"]],
  fmt_num(as.numeric(x[["Mean_years"]])), fmt_num(as.numeric(x[["SD_years"]])),
  fmt_num(as.numeric(x[["Median_years"]])), fmt_num(as.numeric(x[["Mode_years"]])),
  fmt_num(as.numeric(x[["Minimum_years"]])), fmt_num(as.numeric(x[["Maximum_years"]]))
))
new_worst_apa_lines <- character()

if (nrow(new_worst_binary_summary)) {
  x <- new_worst_binary_summary[1, ]
  
  new_worst_apa_lines <- c(new_worst_apa_lines, sprintf(
    "New lifetime-worst trauma after study start was reported by %s of %s participants (%s%%). This binary ever/never summary is descriptive and does not preserve exposure timing.",
    x$N_new_worst_trauma,
    x$N_with_nonmissing_flag,
    fmt_num(x$Percent_new_worst_trauma)
  ))
}

binary_terms_for_apa <- rbind_fill(new_worst_binary_terms)

if (
  nrow(binary_terms_for_apa) &&
  "Term" %in% names(binary_terms_for_apa) &&
  any(binary_terms_for_apa$Term == "new_worst_trauma_binary")
) {
  x <- binary_terms_for_apa[binary_terms_for_apa$Term == "new_worst_trauma_binary", , drop = FALSE][1, ]
  
  new_worst_apa_lines <- c(new_worst_apa_lines, sprintf(
    "In the descriptive ever/never Cox model, new lifetime-worst trauma was associated with HR = %s, 95%% CI [%s, %s], p %s. Because this model does not preserve exposure timing, the time-updated model is the primary inferential analysis.",
    fmt_num(x$HR),
    fmt_num(x$CI_lower),
    fmt_num(x$CI_upper),
    fmt_p(x$p)
  ))
}

time_updated_terms_for_apa <- rbind_fill(new_worst_terms)

if (
  nrow(time_updated_terms_for_apa) &&
  "Term" %in% names(time_updated_terms_for_apa) &&
  any(time_updated_terms_for_apa$Term == "new_worst_trauma_tv")
) {
  x <- time_updated_terms_for_apa[time_updated_terms_for_apa$Term == "new_worst_trauma_tv", , drop = FALSE][1, ]
  
  new_worst_apa_lines <- c(new_worst_apa_lines, sprintf(
    "In the timing-preserving start-stop model, new lifetime-worst trauma was modeled as a time-updated exposure and was associated with HR = %s, 95%% CI [%s, %s], p %s. This model did not adjust for trauma_type_count.",
    fmt_num(x$HR),
    fmt_num(x$CI_lower),
    fmt_num(x$CI_upper),
    fmt_p(x$p)
  ))
}
save_apa_bundle(
  output_dir = out_dir,
  prefix = "RQ2",
  title = "RQ2: Time to onset",
  inventory = INVENTORY,
  model_terms = fixed_terms,
  model_stats = fixed_stats,
  ph = ph_results,
  extra_tables = list(midpoint_audit = midpoint_audit,
    latency_descriptives = latency,
    RMST_results = rmst_results,
    RMST_tau_selection = rmst_tau_selection,
    cutoff_coefficient_results = cutoff_term_results,
    cutoff_model_statistics = cutoff_model_results,
    cutoff_person_time_audit = rbind_fill(cutoff_audits),
    age_at_trauma_sex_interaction_coefficients = age_sex_interaction_term_results,
    age_at_trauma_sex_specific_HRs = age_sex_specific_results,
    age_at_trauma_sex_interaction_tests = age_sex_interaction_test_results,
    age_at_trauma_sex_interaction_model_statistics = age_sex_interaction_model_results,
    age_at_trauma_sex_interaction_PH_diagnostics = age_sex_interaction_ph_results,
    age_at_trauma_sex_interaction_audit = age_sex_interaction_audit_results,
    proximal_sex_interaction_coefficients = proximal_sex_interaction_term_results,
    proximal_sex_specific_HRs = proximal_sex_specific_results,
    proximal_sex_interaction_tests = proximal_sex_interaction_test_results,
    proximal_sex_interaction_model_statistics = proximal_sex_interaction_model_results,
    proximal_sex_interaction_audit = proximal_sex_interaction_audit_results,
    new_worst_trauma_binary_summary = new_worst_binary_summary,
    new_worst_trauma_binary_model_terms = rbind_fill(new_worst_binary_terms),
    new_worst_trauma_binary_model_statistics = rbind_fill(new_worst_binary_stats),
    new_worst_trauma_audit = new_worst_audit,
    extended_age_specific_HRs = extended_results,
    sex_stratification_decisions = sex_decisions,
    PH_residual_plot_audit = ph_plot_audit,
    model_errors = errors,
    model_specifications = rq2_model_specifications
  ),
  additional_lines = c(
    "SCALED SCHOENFELD RESIDUAL PLOTS",
    paste0("Residual plots were saved under ", normalizePath(ph_plot_dir, mustWork = FALSE), ". The initial sex model provides the sex residual plot; subsequent strata(sex) models do not estimate a sex coefficient."),
    "",
    "ONE-TIME SEX PH ASSESSMENT AND GLOBAL HANDLING",
    if (length(sex_decision_lines)) sex_decision_lines else "The initial sex PH decision was not estimable.",
    sprintf("All subsequent Cox models use %s based on this single initial diagnostic; sex is not re-tested model by model.", if (cutoff_sex_mode == "strata") "strata(sex)" else "sex as a covariate"),
    "",
    "EXPLORATORY SEX-EFFECT-MODIFICATION SENSITIVITY MODELS",
    "Two pooled attained-age models tested age-at-trauma-by-sex interactions for most recent and worst trauma while using strata(sex); no sex-specific dataset split was used.",
    paste0("A separate proximal-risk-by-sex interaction was fitted at the prespecified ", SEX_INTERACTION_CUTOFF_MONTHS, "-month cutoff for each trauma origin and onset scenario. BH adjustment was applied separately within the age-at-trauma and proximal-risk interaction families."),
    "Sex-specific HRs and formal interaction tests are in dedicated CSV files and are exploratory rather than primary.",
    "",
    "LATENCY AND RMST RESULTS",
    sprintf("Tau was selected once per time origin and held fixed across bootstrap replicates. Default selection used the %.0fth percentile of the primary duration distribution, capped at common observed support.", 100 * RMST_TAU_QUANTILE),
    latency_apa,
    rmst_apa,
    "",
    "NEW LIFETIME-WORST TRAUMA",
    if (length(new_worst_apa_lines)) new_worst_apa_lines else "New lifetime-worst-trauma analyses were not estimable.",
    "The descriptive binary model is reported before the timing-preserving start-stop model. The time-updated new lifetime-worst-trauma model uses the same sex-handling decision as the cutoff models and intentionally excludes trauma_type_count.",
    "",
    sprintf(
      "The 3-12-month cutoff models use start-stop data with a time-updated within-cutoff indicator, %s, and participant-clustered robust standard errors. Standard cox.zph diagnostics are intentionally not applied to those explicitly time-varying exposure models.",
      if (cutoff_sex_mode == "strata") "sex-specific baseline hazards" else "sex retained as a fixed covariate"
    )
  )
)

progress_log("RQ2 complete", paste0("outputs: ", normalizePath(out_dir, mustWork = FALSE), "; elapsed=", format_elapsed(SCRIPT_STARTED)))
