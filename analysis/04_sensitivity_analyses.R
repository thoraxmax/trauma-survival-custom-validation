# 04_sensitivity_analyses.R
#
# 04_sensitivity_analyses.R
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

# 1. Exclusion of OCD cases first detected within one month.
# SM: Restricted-sample delayed-entry Cox sensitivity models; SP: survival::coxph; IV/DV: trauma predictors / first detected OCD; MF: primary RQ3 trauma-count and trauma-type formulas after excluding events <=1 month; M: assess whether associations are driven by very early detections or possible baseline/prodromal cases.

# 2. Combined psychiatric/first-degree-family confounder adjustment.
# SM: Covariate-adjusted delayed-entry Cox models; SP: survival::coxph; IV/DV: trauma predictors + baseline MINI psychiatric categories + first-degree psychiatric history / OCD; MF: Surv(age_entry, age_exit, event) ~ trauma predictors + baseline_psychiatric_categories + first_degree_psychiatric_history_any + inherited sex handling; M: address reviewers' psychiatric and familial confounding concern in one combined adjustment set using parent-or-sibling history only.

# 3. Site-clustered robust standard errors when a site field is available.
# SM: Delayed-entry Cox with site-clustered robust sandwich SE; SP: survival::coxph; IV/DV: trauma predictors / OCD; MF: Surv(...) ~ trauma predictors + inherited sex handling + cluster(site); M: address multisite dependence with one analytic strategy.

# 4. Inherited RQ2 sex PH decision and Cox-model PH diagnostics.
# SM: RQ2 initial sex Cox model plus PH diagnostic, followed by fixed global sex handling in all sensitivity Cox models; SP: base R import + survival::cox.zph for sensitivity-model predictors; IV/DV: sensitivity-model predictors / OCD; MF: sensitivity models use the RQ2 sex decision without refitting sex locally; M: keep sex handling fixed across sensitivity subsets rather than re-selecting it according to subset-specific p-values.

# 5. Sensitivity coefficient p-values are reported nominally; PH p-values remain raw.
# SM: No multiplicity adjustment is applied to sensitivity coefficients; SP: base R output only; IV/DV: substantive central trauma-test p-values / nominal sensitivity evidence; MF: none; M: reserve BH-FDR for the 10 primary RQ3 trauma-type coefficients and avoid repeated adjustment across sensitivity models.

# 6. Repeated PH-concern summary and targeted extended sensitivity models.
# SM: Extended Cox models with time-varying coefficients for repeatedly PH-flagged exposures; SP: survival::coxph with tt(); IV/DV: PH-flagged exposure x centered log attained age / OCD; MF: Surv(...) ~ exposure + tt(exposure) + covariates + inherited sex handling; M: determine whether sensitivity conclusions remain when non-proportional exposure effects are modeled explicitly.
#
# Reporting output: each run writes Supplementary_Material_Sensitivity_Analyses.txt and Supplementary_Material_Sensitivity_Analyses.docx in r_results/Sensitivity. If officer/flextable are unavailable, a dependency-free DOCX fallback is written.
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
helper_path <- file.path(script_dir, "00_survival_helpers.R")
if (!file.exists(helper_path)) {
  stop("Cannot find 00_survival_helpers.R beside this script: ", helper_path, call. = FALSE)
}
source(helper_path, local = FALSE)
assert_helper_api(21L, c("progress_log", "format_elapsed", "read_survival_data", "save_cox_zph_plots", "save_apa_bundle"))
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


INVENTORY <- c(
  "Exclude OCD cases first detected within one month",
  "Combined adjustment for baseline MINI psychiatric categories and first-degree psychiatric history",
  "Site-clustered robust standard errors",
  "Inherited RQ2 sex PH decision; no sensitivity-specific sex re-test",
  "Nominal sensitivity coefficient p-values only; BH-FDR reserved for primary RQ3 trauma-type coefficients",
  "Repeated PH-concern summary and targeted extended models",
  "APA-style supplementary material export as TXT and DOCX"
)
analysis_inventory("Sensitivity analyses", INVENTORY)

out_dir <- file.path(OUTPUT_ROOT, "Sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ph_plot_dir <- file.path(out_dir, "PH_residual_plots")
dir.create(ph_plot_dir, recursive = TRUE, showWarnings = FALSE)

progress_log("Sensitivity", "loading and validating the prepared dataset")
data <- read_survival_data()
progress_log("Sensitivity", sprintf("loaded N=%d participants; events=%d", nrow(data), sum(data$event == 1, na.rm = TRUE)))
require_columns(data, c("record_id", "sex", "event", "age_entry", "age_exit"), "Sensitivity analyses")

TRAUMA_VARS <- c(
  "sexual_abuse", "physical_abuse", "emotional_abuse", "bullying",
  "accident_experienced_witnessed", "natural_disaster", "illness_injury",
  "kidnapping", "family_violence", "community_violence"
)

screening <- screen_binary_exposures(
  data,
  TRAUMA_VARS,
  max_missing_prop = as.numeric(Sys.getenv("TRAUMA_MAX_MISSING_PROP", "0.20")),
  min_exposed_n = as.integer(Sys.getenv("TRAUMA_MIN_EXPOSED_N", "10")),
  min_unexposed_n = as.integer(Sys.getenv("TRAUMA_MIN_UNEXPOSED_N", "10")),
  min_events_per_level = as.integer(Sys.getenv("TRAUMA_MIN_EVENTS_PER_LEVEL", "3"))
)
eligible_trauma <- screening$Variable[screening$Eligible]

BASELINE_PSYCH <- c(
  "affective_disorder",
  "anxiety_disorders",
  "ptsd",
  "substance_use_disorders",
  "neurodevelopmental_disorders",
  "behavioral_disorders",
  "eating_disorders"
)

FIRST_DEGREE_HISTORY_CANDIDATES <- c(
  "first_degree_psychiatric_history_any",
  "first_relative_psychiatric_history_any",
  "parent_sibling_psychiatric_history_any",
  "parent_psychiatric_history_any"
)
FIRST_DEGREE_HISTORY_VAR <- first_existing(data, FIRST_DEGREE_HISTORY_CANDIDATES)
if (is.na(FIRST_DEGREE_HISTORY_VAR) || !nzchar(FIRST_DEGREE_HISTORY_VAR)) {
  progress_log("Sensitivity", "no first-degree psychiatric-history variable found; confounder-adjusted models use baseline psychiatric categories only")
  FIRST_DEGREE_HISTORY_VAR <- ""
} else if (!identical(FIRST_DEGREE_HISTORY_VAR, "first_degree_psychiatric_history_any")) {
  progress_log("Sensitivity", paste0("using legacy first-degree psychiatric-history fallback: ", FIRST_DEGREE_HISTORY_VAR))
}

term_tables <- list()
model_tables <- list()
ph_tables <- list()
omnibus_tables <- list()
audit_tables <- list()
errors <- list()
skipped <- list()
extended_rows <- list()
extended_terms <- list()
sex_decision_tables <- list()
ph_plot_audits <- list()

record_error <- function(model, e) {
  errors[[length(errors) + 1L]] <<- data.frame(Model = model, Error = conditionMessage(e), stringsAsFactors = FALSE)
}
record_skip <- function(analysis, reason) {
  skipped[[length(skipped) + 1L]] <<- data.frame(Analysis = analysis, Reason = reason, stringsAsFactors = FALSE)
}

usable_adjustments <- function(d, vars) {
  vars <- intersect(vars, names(d))
  vars[vapply(vars, function(v) {
    x <- d[[v]]
    x <- x[!is.na(x)]
    length(unique(x)) >= 2L
  }, logical(1))]
}

# Sex handling is inherited from the RQ2 one-time sex association/PH assessment.
# Sensitivity models do not refit or re-test sex in subsets.
load_rq2_sex_decision <- function(context = "Sensitivity") {
  decision_path <- trimws(Sys.getenv(
    "RQ2_SEX_DECISION_PATH",
    file.path(OUTPUT_ROOT, "RQ2", "RQ2_global_sex_handling_decision.csv")
  ))
  if (!file.exists(decision_path)) {
    stop(
      context,
      " requires the RQ2 global sex-handling decision. Run 02_RQ2_time_to_onset.R first, ",
      "or set RQ2_SEX_DECISION_PATH to the saved RQ2_global_sex_handling_decision.csv file.",
      call. = FALSE
    )
  }
  decision <- read.csv(decision_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(decision)) stop("The RQ2 sex-decision file is empty: ", decision_path, call. = FALSE)
  decision <- decision[1L, , drop = FALSE]
  decision$Decision_source_file <- normalizePath(decision_path, mustWork = FALSE)
  decision$Decision_used_by <- context
  decision
}

rq2_sex_decision <- load_rq2_sex_decision("Sensitivity")
SEX_PH_P <- suppressWarnings(as.numeric(rq2_sex_decision$Sex_PH_p[[1L]]))
SEX_STRATIFY <- isTRUE(rq2_sex_decision$Sex_stratified[[1L]]) ||
  tolower(as.character(rq2_sex_decision$Sex_stratified[[1L]])) %in% c("true", "t", "1", "yes")
GLOBAL_SEX_STRATA <- if (SEX_STRATIFY) "sex" else character()
GLOBAL_SEX_COVARIATES <- if (SEX_STRATIFY) character() else "sex"
rq2_sex_decision$Model <- "Sensitivity inherited RQ2 sex PH decision"
rq2_sex_decision$Final_formula <- if (SEX_STRATIFY) "Sensitivity Cox models: ... + strata(sex)" else "Sensitivity Cox models: ... + sex"
rq2_sex_decision$Final_sex_handling <- if (SEX_STRATIFY) "strata(sex) inherited from RQ2" else "sex covariate inherited from RQ2"
sex_decision_tables[[1L]] <- rq2_sex_decision
progress_log("Sensitivity sex handling", paste0("inherited RQ2 sex PH p ", fmt_p(SEX_PH_P), "; sensitivity models use ", if (SEX_STRATIFY) "strata(sex)" else "sex as a covariate"))

run_sensitivity_pair <- function(label, d, adjustments = character(), cluster = NULL) {
  started <- Sys.time()
  progress_log("Sensitivity", paste0("starting: ", label))
  adjustments <- usable_adjustments(d, adjustments)
  audit_tables[[length(audit_tables) + 1L]] <<- data.frame(
    Sensitivity = label,
    Input_N = nrow(d),
    Input_events = sum(d$event == 1, na.rm = TRUE),
    Eligible_trauma_variables = paste(eligible_trauma, collapse = "; "),
    Adjustment_variables = paste(adjustments, collapse = "; "),
    Cluster_variable = cluster %||% "",
    Sex_handling = if (SEX_STRATIFY) "strata(sex)" else "sex covariate",
    Sex_decision_source = "Inherited one-time RQ2 sex PH assessment",
    stringsAsFactors = FALSE
  )

  fit_sensitivity_model <- function(model_name, exposures) {
    covariates <- c(exposures, adjustments, GLOBAL_SEX_COVARIATES)
    fit <- fit_age_scale_cox(d, covariates, strata = GLOBAL_SEX_STRATA, cluster = cluster, model_name = model_name)
    ph <- cox_ph_table(
      fit, model_name,
      main_exposures = exposures,
      nuisance_terms = c(adjustments, GLOBAL_SEX_COVARIATES)
    )
    ph$Model_phase <- "Final model"
    ph_plot_audit <- save_cox_zph_plots(fit, ph_plot_dir, model_name, "Final model", ph_table = ph)
    ph_plot_audits[[length(ph_plot_audits) + 1L]] <<- ph_plot_audit
    term_tables[[length(term_tables) + 1L]] <<- extract_cox_terms(fit, model_name)
    model_tables[[length(model_tables) + 1L]] <<- extract_model_stats(fit, model_name)
    ph_tables[[length(ph_tables) + 1L]] <<- ph
    fit
  }

  if (length(eligible_trauma)) {
    model_name <- paste0(label, ": trauma-type model")
    tryCatch({
      fit <- fit_sensitivity_model(model_name, eligible_trauma)
      joint <- wald_joint_test(fit, eligible_trauma, paste0(label, ": trauma-type joint Wald test"))
      joint$Model <- model_name
      omnibus_tables[[length(omnibus_tables) + 1L]] <<- joint
    }, error = function(e) record_error(model_name, e))
  } else {
    record_skip(paste0(label, ": trauma-type model"), "No trauma variables passed screening")
  }

  if ("trauma_type_count" %in% names(d)) {
    model_name <- paste0(label, ": trauma-type-count model")
    tryCatch({
      fit <- fit_sensitivity_model(model_name, "trauma_type_count")
      joint <- wald_joint_test(fit, "trauma_type_count", paste0(label, ": trauma-count Wald test"))
      joint$Model <- model_name
      omnibus_tables[[length(omnibus_tables) + 1L]] <<- joint
    }, error = function(e) record_error(model_name, e))
  } else {
    record_skip(paste0(label, ": trauma-type-count model"), "trauma_type_count is absent")
  }
  progress_log("Sensitivity", sprintf("completed: %s; elapsed=%s", label, format_elapsed(started)))
}

# 1. Exclude one-month cases.
if ("first_positive_months" %in% names(data)) {
  exclude_one_month <- data[!(data$event == 1 & is.finite(data$first_positive_months) & data$first_positive_months <= 1), , drop = FALSE]
  run_sensitivity_pair("Exclude cases detected within 1 month", exclude_one_month)
} else {
  record_skip("Exclude cases detected within 1 month", "first_positive_months is absent")
}

# 2. Combined baseline psychiatric and first-degree-family-history adjustment.
family_history_adjustment <- if (!is.na(FIRST_DEGREE_HISTORY_VAR) && nzchar(FIRST_DEGREE_HISTORY_VAR)) FIRST_DEGREE_HISTORY_VAR else character()
combined_adjustments <- usable_adjustments(data, c(BASELINE_PSYCH, family_history_adjustment))
if (length(combined_adjustments)) {
  run_sensitivity_pair("Confounder-adjusted psychiatric and first-degree family history", data, adjustments = combined_adjustments)
} else {
  record_skip(
    "Confounder-adjusted psychiatric and first-degree family history",
    "No usable baseline psychiatric or first-degree family-history adjustment variables were available"
  )
}

# 3. Site-clustered robust standard errors only; no site-fixed-effects model.
site_var <- first_existing(data, c("site", "site_id", "record_dag_name", "dag", "data_access_group"))
if (!is.na(site_var)) {
  data[[site_var]] <- droplevels(factor(data[[site_var]]))
  if (nlevels(data[[site_var]]) >= 2L) {
    run_sensitivity_pair("Site-clustered robust SE", data, cluster = site_var)
  } else {
    record_skip("Site-clustered robust SE", paste(site_var, "has fewer than two observed levels"))
  }
} else {
  record_skip("Site-clustered robust SE", "No prepared site variable was found")
}


# -----------------------------------------------------------------------------
# Repeated central-exposure PH concerns and targeted extended sensitivity models
# -----------------------------------------------------------------------------
ph_results_pre <- rbind_fill(ph_tables)
repeated_ph <- data.frame()
if (nrow(ph_results_pre)) {
  central <- c(eligible_trauma, "trauma_type_count")
  final_phase <- if ("Model_phase" %in% names(ph_results_pre)) {
    is.na(ph_results_pre$Model_phase) | ph_results_pre$Model_phase == "Final model"
  } else {
    rep(TRUE, nrow(ph_results_pre))
  }
  concerns <- ph_results_pre[
    final_phase & ph_results_pre$PH_term %in% central & ph_results_pre$PH_concern %in% TRUE,
    c("Model", "PH_term", "PH_p"),
    drop = FALSE
  ]
  if (nrow(concerns)) {
    repeated_ph <- aggregate(
      Model ~ PH_term,
      data = concerns,
      FUN = function(x) length(unique(x))
    )
    names(repeated_ph) <- c("Exposure", "Sensitivity_models_with_PH_concern")
    repeated_ph$Targeted_extended_model <- repeated_ph$Sensitivity_models_with_PH_concern >= 2L

    repeated_exposures <- repeated_ph$Exposure[repeated_ph$Targeted_extended_model]
    for (exposure in repeated_exposures) {
      adjustment <- character()
      if (exposure %in% eligible_trauma) adjustment <- c(setdiff(eligible_trauma, exposure), psych_available)
      if (exposure == "trauma_type_count") adjustment <- psych_available
      model_name <- paste("Targeted extended sensitivity:", exposure)
      adjustment <- usable_adjustments(data, adjustment)
      ext <- tryCatch(
        fit_extended_age_effect(
          data,
          exposure = exposure,
          adjustment_terms = c(adjustment, GLOBAL_SEX_COVARIATES),
          strata = GLOBAL_SEX_STRATA,
          model_name = model_name
        ),
        error = function(e) { record_error(model_name, e); NULL }
      )
      if (!is.null(ext)) {
        extended_rows[[length(extended_rows) + 1L]] <- ext$age_specific
        extended_terms[[length(extended_terms) + 1L]] <- ext$interaction
        model_tables[[length(model_tables) + 1L]] <- extract_model_stats(ext$fit, model_name)
      }
    }
  }
}

# -----------------------------------------------------------------------------
# Final output tables
# -----------------------------------------------------------------------------
coefficient_results <- rbind_fill(c(term_tables, extended_terms))
model_results <- rbind_fill(model_tables)
ph_results <- rbind_fill(ph_tables)
omnibus_results <- rbind_fill(omnibus_tables)
audit_results <- rbind_fill(audit_tables)
error_results <- rbind_fill(errors)
skipped_results <- rbind_fill(skipped)
extended_results <- rbind_fill(extended_rows)
sex_decisions <- rbind_fill(sex_decision_tables)
ph_plot_audit <- rbind_fill(ph_plot_audits)

if (nrow(coefficient_results)) {
  coefficient_results$Sensitivity_p_values_adjusted <- FALSE
  coefficient_results$Multiplicity_note <- "Nominal sensitivity p-values; BH-FDR is reserved for the primary RQ3 trauma-type coefficient family."
} else {
  coefficient_results$Sensitivity_p_values_adjusted <- logical(0)
  coefficient_results$Multiplicity_note <- character(0)
}

sex_decision_lines <- character()
if (nrow(sex_decisions)) {
  sex_decision_lines <- apply(sex_decisions, 1, function(x) sprintf(
    "%s: sex PH test chi-square(%s) = %s, p %s; final handling: %s.",
    x[["Model"]], fmt_num(as.numeric(x[["Sex_PH_df"]]), 0),
    fmt_num(as.numeric(x[["Sex_PH_chisq"]])), fmt_p(as.numeric(x[["Sex_PH_p"]])),
    x[["Final_sex_handling"]]
  ))
}

robustness_lines <- character()
if (nrow(coefficient_results)) {
  key <- coefficient_results[coefficient_results$Term %in% c(eligible_trauma, "trauma_type_count"), , drop = FALSE]
  if (nrow(key)) {
    robustness_lines <- apply(key, 1, function(x) sprintf(
      "%s, %s: HR = %s, 95%% CI [%s, %s], z = %s, p %s.",
      x[["Model"]], x[["Term"]], fmt_num(as.numeric(x[["HR"]])),
      fmt_num(as.numeric(x[["CI_lower"]])), fmt_num(as.numeric(x[["CI_upper"]])),
      fmt_num(as.numeric(x[["z"]])), fmt_p(as.numeric(x[["p"]]))
    ))
  }
}

extended_lines <- character()
if (nrow(extended_results)) {
  extended_lines <- apply(extended_results, 1, function(x) sprintf(
    "%s at age %s: HR = %s, 95%% CI [%s, %s], z = %s, p %s.",
    x[["Model"]], fmt_num(as.numeric(x[["Age"]]), 1), fmt_num(as.numeric(x[["HR"]])),
    fmt_num(as.numeric(x[["CI_lower"]])), fmt_num(as.numeric(x[["CI_upper"]])),
    fmt_num(as.numeric(x[["z"]])), fmt_p(as.numeric(x[["p"]]))
  ))
}


adjusted_family_formula_term <- if (nzchar(FIRST_DEGREE_HISTORY_VAR)) {
  paste("baseline psychiatric categories +", FIRST_DEGREE_HISTORY_VAR)
} else {
  "baseline psychiatric categories"
}

sensitivity_model_specifications <- data.frame(
  Analysis = c(
    "Inherited RQ2 sex PH decision",
    "Early-case exclusion trauma-type sensitivity models",
    "Early-case exclusion trauma-count sensitivity models",
    "Combined psychiatric/first-degree-family-history adjustment",
    "Site-clustered robust sensitivity where available",
    "Age-varying coefficient extension for repeated PH concern"
  ),
  R_package = c("base R", rep("survival", 5)),
  R_function = c("read.csv", rep("survival::coxph", 5)),
  Model_type = c(
    "Imported RQ2 one-time sex association/PH decision; no sensitivity-specific sex model fitted",
    "Restricted-sample delayed-entry attained-age Cox using the inherited sex decision",
    "Restricted-sample delayed-entry attained-age Cox using the inherited sex decision",
    "Covariate-adjusted delayed-entry attained-age Cox using the inherited sex decision",
    "Delayed-entry attained-age Cox with clustered sandwich variance and inherited sex decision",
    "Delayed-entry attained-age Cox with exposure x centered log attained-age term and inherited sex decision"
  ),
  Formula = c(
    "RQ2: Surv(age_entry, age_exit, event) ~ sex; sensitivity scripts inherit the resulting sex handling",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ eligible trauma indicators + strata(sex)" else "Surv(age_entry, age_exit, event) ~ eligible trauma indicators + sex",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ trauma_type_count + strata(sex)" else "Surv(age_entry, age_exit, event) ~ trauma_type_count + sex",
    if (SEX_STRATIFY) paste0("Surv(age_entry, age_exit, event) ~ trauma predictors + ", adjusted_family_formula_term, " + strata(sex)") else paste0("Surv(age_entry, age_exit, event) ~ trauma predictors + ", adjusted_family_formula_term, " + sex"),
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ trauma predictors + strata(sex) + cluster(site)" else "Surv(age_entry, age_exit, event) ~ trauma predictors + sex + cluster(site)",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ exposure + tt(exposure) + strata(sex)" else "Surv(age_entry, age_exit, event) ~ exposure + tt(exposure) + sex"
  ),
  Ties = c(NA, rep("Efron", 5)),
  stringsAsFactors = FALSE
)
write_csv_safe(sensitivity_model_specifications, file.path(out_dir, "Sensitivity_model_specifications.csv"))

save_apa_bundle(
  output_dir = out_dir,
  prefix = "Sensitivity",
  title = "Sensitivity analyses",
  inventory = INVENTORY,
  model_terms = coefficient_results,
  model_stats = model_results,
  ph = ph_results,
  extra_tables = list(
    trauma_screening = screening,
    analysis_audits = audit_results,
    omnibus_Wald_tests = omnibus_results,
    repeated_PH_concerns = repeated_ph,
    targeted_extended_age_specific_HRs = extended_results,
    sex_stratification_decisions = sex_decisions,
    PH_residual_plot_audit = ph_plot_audit,
    skipped_analyses = skipped_results,
    model_errors = error_results,
    model_specifications = sensitivity_model_specifications
  ),
  additional_lines = c(
    "SCALED SCHOENFELD RESIDUAL PLOTS",
    paste0("Residual plots were saved under ", normalizePath(ph_plot_dir, mustWork = FALSE), ". Sex residual diagnostics are inherited from RQ2; sensitivity models save plots only for substantive fixed-coefficient predictors."),
    "",
    "INHERITED RQ2 SEX PH DECISION AND GLOBAL HANDLING",
    if (length(sex_decision_lines)) sex_decision_lines else "The RQ2 sex PH decision was not available.",
    "",
    "KEY ROBUSTNESS RESULTS",
    if (length(robustness_lines)) robustness_lines else "No key sensitivity coefficient was estimable.",
    "",
    "TARGETED EXTENDED SENSITIVITY RESULTS",
    if (length(extended_lines)) extended_lines else "No central exposure showed PH concern in at least two fixed-coefficient sensitivity models, or no targeted extended model was estimable.",
    "",
    "Sex is tested for PH once in RQ2. The resulting sex handling is inherited and held fixed across all sensitivity analyses, without age-based restricted subsets; sex is not re-tested within each sensitivity model. Sensitivity coefficient p-values are nominal; BH-FDR is reserved for the primary RQ3 trauma-type coefficient family."
  )
)

# =============================================================================
# MANUSCRIPT SUPPLEMENT EXPORT
# Always write a dependency-free APA-style TXT. If officer + flextable are
# installed, also write a formatted DOCX suitable for Supplementary Information.
# =============================================================================
supplement_txt_path <- file.path(out_dir, "Supplementary_Material_Sensitivity_Analyses.txt")
supplement_docx_path <- file.path(out_dir, "Supplementary_Material_Sensitivity_Analyses.docx")

supp_select <- function(x, columns) {
  if (!is.data.frame(x) || !nrow(x)) return(data.frame(Note = "No estimable results.", stringsAsFactors = FALSE))
  keep <- intersect(columns, names(x))
  if (!length(keep)) return(data.frame(Note = "No requested columns were available.", stringsAsFactors = FALSE))
  x[, keep, drop = FALSE]
}

supp_format_table <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) return(data.frame(Note = "No estimable results.", stringsAsFactors = FALSE))
  out <- x
  p_cols <- grep("(^p$|_p$|^PH_p$|p_FDR|_p_value$)", names(out), value = TRUE)
  num_2_cols <- intersect(
    c("B", "SE", "HR", "CI_lower", "CI_upper", "z", "PH_chisq", "Wald_chisq", "LR_chisq", "Score_chisq", "Concordance", "AIC", "Age"),
    names(out)
  )
  int_cols <- intersect(c("N", "Events", "Input_N", "Input_events", "Analysis_rows", "PH_df", "Wald_df", "LR_df", "Score_df"), names(out))
  for (nm in num_2_cols) {
    values <- suppressWarnings(as.numeric(out[[nm]]))
    out[[nm]] <- ifelse(is.finite(values), formatC(values, digits = 2L, format = "f"), "")
  }
  for (nm in int_cols) {
    values <- suppressWarnings(as.numeric(out[[nm]]))
    out[[nm]] <- ifelse(is.finite(values), formatC(values, digits = 0L, format = "f"), "")
  }
  for (nm in p_cols) {
    values <- suppressWarnings(as.numeric(out[[nm]]))
    out[[nm]] <- ifelse(is.finite(values), fmt_p(values), "")
  }
  for (nm in names(out)) {
    if (is.logical(out[[nm]])) out[[nm]] <- ifelse(is.na(out[[nm]]), "", ifelse(out[[nm]], "Yes", "No"))
  }
  out
}

supp_table_lines <- function(label, title, table) {
  formatted <- supp_format_table(table)
  body <- capture.output(utils::write.table(formatted, sep = "\t", row.names = FALSE, quote = FALSE, na = ""))
  c("", label, title, body)
}

supplement_audit <- supp_select(
  audit_results,
  c("Sensitivity", "Input_N", "Input_events", "Adjustment_variables", "Cluster_variable", "Sex_handling", "Sex_decision_source")
)
supplement_coefficients <- supp_select(
  coefficient_results,
  c("Model", "Term", "N", "Events", "HR", "CI_lower", "CI_upper", "p", "Sensitivity_p_values_adjusted", "Multiplicity_note")
)
supplement_omnibus <- supp_select(
  omnibus_results,
  c("Model", "Test", "Wald_chisq", "df", "p", "N", "Events")
)
supplement_ph <- supp_select(
  ph_results,
  c("Model", "Model_phase", "PH_term", "PH_chisq", "PH_df", "PH_p", "PH_concern", "Planned_action")
)
supplement_extended <- supp_select(
  extended_results,
  c("Model", "Exposure", "Age", "N", "Events", "HR", "CI_lower", "CI_upper", "p")
)
supplement_repeated_ph <- supp_select(
  repeated_ph,
  c("Exposure", "Models_flagged", "Number_models_flagged", "Targeted_extended_model")
)
supplement_skipped <- supp_select(skipped_results, c("Analysis", "Reason"))
supplement_errors <- supp_select(error_results, c("Model", "Error"))

sex_handling_sentence <- if (SEX_STRATIFY) {
  paste0(
    "Sex was evaluated once in RQ2 using a delayed-entry attained-age Cox model and rank-transformed Schoenfeld residual test (p ",
    fmt_p(SEX_PH_P), "). Because the proportional-hazards diagnostic indicated non-proportionality, all sensitivity Cox models inherited strata(sex), allowing sex-specific baseline hazards while estimating common exposure effects."
  )
} else {
  paste0(
    "Sex was evaluated once in RQ2 using a delayed-entry attained-age Cox model and rank-transformed Schoenfeld residual test (p ",
    fmt_p(SEX_PH_P), "). Because the proportional-hazards diagnostic did not indicate non-proportionality, sex was retained as an inherited fixed covariate in all sensitivity Cox models."
  )
}

supp_key <- if (nrow(coefficient_results)) {
  keep <- coefficient_results$Term %in% c(eligible_trauma, "trauma_type_count") & coefficient_results$p < .05
  coefficient_results[keep, , drop = FALSE]
} else data.frame()

supp_key_lines <- if (nrow(supp_key)) {
  apply(supp_key, 1, function(x) sprintf(
    "%s, %s: HR = %s, 95%% CI [%s, %s], p %s.",
    x[["Model"]], x[["Term"]], fmt_num(as.numeric(x[["HR"]])),
    fmt_num(as.numeric(x[["CI_lower"]])), fmt_num(as.numeric(x[["CI_upper"]])),
    fmt_p(as.numeric(x[["p"]]))
  ))
} else {
  "No central sensitivity coefficient met the nominal p < .05 threshold; complete estimates are reported in Table S2."
}

skipped_sentence <- if (nrow(skipped_results)) {
  paste0(
    "Analyses that could not be estimated because the required prepared variable or sufficient variation was unavailable are documented in Table S7 (n = ",
    nrow(skipped_results), ")."
  )
} else {
  "All prespecified sensitivity analyses with available variables were estimable."
}

supplement_lines <- c(
  "SUPPLEMENTARY INFORMATION",
  "Sensitivity Analyses",
  "",
  "Overview",
  "These supplementary analyses evaluated whether the primary trauma findings were sensitive to early-case exclusion, combined psychiatric/first-degree-family confounder adjustment, site clustering when available, multiplicity control, and non-proportional hazards.",
  "",
  "Methods",
  "All fixed-coefficient sensitivity analyses used delayed-entry Cox regression with attained age as the time scale and Efron handling of tied event times. Participants entered the risk set at age at study entry and were censored at their final completed OCD assessment.",
  sex_handling_sentence,
  "Sensitivity analyses excluded OCD cases first detected within one month; used one combined adjustment set containing available baseline MINI psychiatric categories and first-degree psychiatric history (parent or sibling); and evaluated site clustering when a prepared site field was available.",
  "The proportional-hazards assumption for fixed-coefficient sensitivity models was assessed using rank-transformed scaled Schoenfeld residuals. PH p-values were treated as diagnostics and were not FDR-adjusted. Central trauma exposures showing repeated PH concern were evaluated in targeted models containing an exposure-by-centered-log-attained-age term, with age-specific hazard ratios derived from the fitted coefficients.",
  "Sensitivity coefficient p-values were reported nominally; BH-FDR adjustment was reserved for the primary RQ3 trauma-type coefficient family and was not repeated across sensitivity models. All tests were two-sided.",
  "",
  "Results",
  paste0("The sensitivity-analysis dataset contained ", nrow(data), " participants and ", sum(data$event == 1, na.rm = TRUE), " first-detected OCD events before analysis-specific exclusions or complete-case restrictions."),
  supp_key_lines,
  skipped_sentence,
  "Complete estimates, diagnostics, and analysis audits are provided in Tables S1-S8 below."
)

supplement_lines <- c(
  supplement_lines,
  supp_table_lines("Table S1", "Sensitivity-analysis sample, adjustment, and sex-handling audit.", supplement_audit),
  supp_table_lines("Table S2", "Sensitivity-model coefficient estimates. P-values are nominal; BH-FDR is not repeated in sensitivity models.", supplement_coefficients),
  supp_table_lines("Table S3", "Omnibus Wald tests for trauma predictor sets.", supplement_omnibus),
  supp_table_lines("Table S4", "Proportional-hazards diagnostics for fixed-coefficient sensitivity models.", supplement_ph),
  supp_table_lines("Table S5", "Repeated PH concerns used to identify targeted extended sensitivity models.", supplement_repeated_ph),
  supp_table_lines("Table S6", "Age-specific hazard ratios from targeted exposure-by-log-attained-age models.", supplement_extended),
  supp_table_lines("Table S7", "Prespecified sensitivity analyses that were not estimable or were skipped.", supplement_skipped),
  supp_table_lines("Table S8", "Sensitivity-model errors, if any.", supplement_errors),
  "",
  "Note. HR = hazard ratio; CI = confidence interval; PH = proportional hazards. Sex stratification does not split the analytic sample; it permits separate baseline hazard functions by sex while retaining a common exposure coefficient."
)
writeLines(supplement_lines, con = supplement_txt_path, useBytes = TRUE)
progress_log("Supplementary material", paste0("APA-style TXT saved: ", normalizePath(supplement_txt_path, mustWork = FALSE)))

write_supplement_docx <- function(path) {
  if (!requireNamespace("officer", quietly = TRUE) || !requireNamespace("flextable", quietly = TRUE)) {
    progress_log("Supplementary material", "officer/flextable not available; writing dependency-free APA-style DOCX fallback")
    return(write_basic_docx_from_lines(path, supplement_lines))
  }

  add_para <- function(doc, text, size = 11, bold = FALSE, italic = FALSE, align = "left", after = 4) {
    props <- officer::fp_text(font.family = "Times New Roman", font.size = size, bold = bold, italic = italic)
    para <- officer::fpar(
      officer::ftext(text, props),
      fp_p = officer::fp_par(text.align = align, padding.bottom = after)
    )
    officer::body_add_fpar(doc, para)
  }

  add_table <- function(doc, number, title, x) {
    doc <- add_para(doc, number, size = 10, bold = TRUE, after = 0)
    doc <- add_para(doc, title, size = 10, italic = TRUE, after = 3)
    x <- supp_format_table(x)
    ft <- flextable::flextable(x)
    ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
    ft <- flextable::fontsize(ft, size = 8.5, part = "all")
    ft <- flextable::bold(ft, part = "header")
    if ("theme_apa" %in% getNamespaceExports("flextable")) {
      ft <- flextable::theme_apa(ft)
    } else {
      ft <- flextable::theme_booktabs(ft)
    }
    ft <- flextable::autofit(ft)
    doc <- flextable::body_add_flextable(doc, ft)
    add_para(doc, "", size = 5, after = 2)
  }

  doc <- officer::read_docx()
  doc <- add_para(doc, "Supplementary Information", size = 14, bold = TRUE, align = "center", after = 2)
  doc <- add_para(doc, "Sensitivity Analyses", size = 13, bold = TRUE, align = "center", after = 8)
  doc <- add_para(doc, "Overview", size = 11, bold = TRUE)
  doc <- add_para(doc, "These supplementary analyses evaluated whether the primary trauma findings were sensitive to early-case exclusion, combined psychiatric/first-degree-family confounder adjustment, site clustering when available, multiplicity control, and non-proportional hazards.")
  doc <- add_para(doc, "Methods", size = 11, bold = TRUE)
  for (line in supplement_lines[match("All fixed-coefficient sensitivity analyses used delayed-entry Cox regression with attained age as the time scale and Efron handling of tied event times. Participants entered the risk set at age at study entry and were censored at their final completed OCD assessment.", supplement_lines):match("Sensitivity coefficient p-values were reported nominally; BH-FDR adjustment was reserved for the primary RQ3 trauma-type coefficient family and was not repeated across sensitivity models. All tests were two-sided.", supplement_lines)]) {
    doc <- add_para(doc, line)
  }
  doc <- add_para(doc, "Results", size = 11, bold = TRUE)
  result_start <- match(paste0("The sensitivity-analysis dataset contained ", nrow(data), " participants and ", sum(data$event == 1, na.rm = TRUE), " first-detected OCD events before analysis-specific exclusions or complete-case restrictions."), supplement_lines)
  result_end <- match("Complete estimates, diagnostics, and analysis audits are provided in Tables S1-S8 below.", supplement_lines)
  if (is.finite(result_start) && is.finite(result_end)) {
    for (line in supplement_lines[result_start:result_end]) doc <- add_para(doc, line)
  }
  doc <- add_table(doc, "Table S1", "Sensitivity-analysis sample, adjustment, and sex-handling audit.", supplement_audit)
  doc <- add_table(doc, "Table S2", "Sensitivity-model coefficient estimates.", supplement_coefficients)
  doc <- add_table(doc, "Table S3", "Omnibus Wald tests for trauma predictor sets.", supplement_omnibus)
  doc <- add_table(doc, "Table S4", "Proportional-hazards diagnostics for fixed-coefficient sensitivity models.", supplement_ph)
  doc <- add_table(doc, "Table S5", "Repeated PH concerns used to identify targeted extended sensitivity models.", supplement_repeated_ph)
  doc <- add_table(doc, "Table S6", "Age-specific hazard ratios from targeted exposure-by-log-attained-age models.", supplement_extended)
  doc <- add_table(doc, "Table S7", "Prespecified sensitivity analyses that were not estimable or were skipped.", supplement_skipped)
  doc <- add_table(doc, "Table S8", "Sensitivity-model errors, if any.", supplement_errors)
  doc <- add_para(doc, "Note. HR = hazard ratio; CI = confidence interval; PH = proportional hazards. Sex stratification does not split the analytic sample; it permits separate baseline hazard functions by sex while retaining a common exposure coefficient.", size = 9, italic = TRUE)
  print(doc, target = path)
  TRUE
}


write_basic_docx_from_lines <- function(path, lines) {
  xml_escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  paragraph_xml <- function(text, style = "BodyText", bold = FALSE, italic = FALSE, size = 22L, align = "left") {
    text <- ifelse(is.na(text), "", as.character(text))
    tab_parts <- strsplit(text, "\t", fixed = TRUE)[[1L]]
    runs <- paste(vapply(seq_along(tab_parts), function(i) {
      prefix <- if (i == 1L) "" else "<w:tab/>"
      paste0(
        "<w:r><w:rPr><w:rFonts w:ascii=\"Times New Roman\" w:hAnsi=\"Times New Roman\"/>",
        if (bold) "<w:b/>" else "",
        if (italic) "<w:i/>" else "",
        "<w:sz w:val=\"", size, "\"/></w:rPr>", prefix,
        "<w:t xml:space=\"preserve\">", xml_escape(tab_parts[[i]]), "</w:t></w:r>"
      )
    }, character(1)), collapse = "")
    paste0(
      "<w:p><w:pPr><w:pStyle w:val=\"", style, "\"/><w:jc w:val=\"", align, "\"/>",
      "<w:spacing w:after=\"120\" w:line=\"480\" w:lineRule=\"auto\"/></w:pPr>",
      runs,
      "</w:p>"
    )
  }
  tmp <- tempfile("docx_build_")
  dir.create(file.path(tmp, "_rels"), recursive = TRUE)
  dir.create(file.path(tmp, "word"), recursive = TRUE)
  writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
</Types>', file.path(tmp, "[Content_Types].xml"), useBytes = TRUE)
  writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>', file.path(tmp, "_rels", ".rels"), useBytes = TRUE)
  styles <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:style w:type="paragraph" w:styleId="BodyText"><w:name w:val="Body Text"/><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="24"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:b/><w:sz w:val="32"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:b/><w:sz w:val="26"/></w:rPr></w:style>
</w:styles>'
  writeLines(styles, file.path(tmp, "word", "styles.xml"), useBytes = TRUE)
  body <- character()
  for (line in lines) {
    if (!nzchar(line)) {
      body <- c(body, paragraph_xml("", size = 22L))
    } else if (line %in% c("SUPPLEMENTARY INFORMATION", "Sensitivity Analyses")) {
      body <- c(body, paragraph_xml(line, style = "Title", bold = TRUE, size = 32L, align = "center"))
    } else if (line %in% c("Overview", "Methods", "Results") || grepl("^Table S[0-9]", line)) {
      body <- c(body, paragraph_xml(line, style = "Heading1", bold = TRUE, size = 24L))
    } else {
      body <- c(body, paragraph_xml(line, style = "BodyText", size = 22L))
    }
  }
  sect <- '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr>'
  doc_xml <- paste0('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>', paste(body, collapse = "\n"), sect, '</w:body></w:document>')
  writeLines(doc_xml, file.path(tmp, "word", "document.xml"), useBytes = TRUE)
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(tmp)
  if (file.exists(path)) unlink(path)
  ok <- tryCatch({ utils::zip(zipfile = path, files = c("[Content_Types].xml", "_rels", "word"), flags = "-r9Xq"); TRUE }, error = function(e) FALSE)
  ok && file.exists(path) && is.finite(file.info(path)$size) && file.info(path)$size > 0
}

if (isTRUE(write_supplement_docx(supplement_docx_path))) {
  progress_log("Supplementary material", paste0("DOCX saved: ", normalizePath(supplement_docx_path, mustWork = FALSE)))
}

progress_log("Sensitivity complete", paste0("outputs: ", normalizePath(out_dir, mustWork = FALSE), "; elapsed=", format_elapsed(SCRIPT_STARTED)))
