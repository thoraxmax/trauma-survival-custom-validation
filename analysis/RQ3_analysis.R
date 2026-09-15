#
# 03_RQ3_risk_factors.R
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

# 1. Prespecified trauma-variable eligibility and sparsity screening.
# SM: Descriptive variable-screening procedure; SP: base R; IV/DV: trauma indicators + event support / eligibility for modeling; MF: none; M: exclude predictors with excessive missingness, sparse exposure groups, or inadequate event support before multivariable modeling.

# 2. Primary sex-stratified delayed-entry multivariable Cox model for trauma types.
# SM: Delayed-entry attained-age multivariable Cox PH with sex-stratified baseline hazards; SP: survival::coxph; IV/DV: mutually adjusted trauma-type indicators / first detected OCD; MF: Surv(age_entry, age_exit, event) ~ trauma_type_1 + ... + trauma_type_k + strata(sex); M: estimate independent trauma-type associations while allowing sex-specific baseline hazards.

# 3. Primary sex-stratified delayed-entry Cox model for trauma-type count.
# SM: Delayed-entry attained-age Cox PH with sex-stratified baseline hazards; SP: survival::coxph; IV/DV: number of trauma types / first detected OCD; MF: Surv(age_entry, age_exit, event) ~ trauma_type_count + strata(sex); M: test whether cumulative trauma burden predicts OCD hazard.

# 4. Inherited global sex handling from RQ2.
# SM: No RQ3 sex model is fitted; SP: base R import of the RQ2 decision table; IV/DV: none; MF: RQ2 sex PH decision determines sex covariate versus strata(sex); M: avoid repeated data-dependent switching and keep the sex association/PH test in RQ2.

# 5. Raw rank-Schoenfeld PH diagnostics and scaled Schoenfeld residual plots for RQ3 substantive models.
# SM: Proportional-hazards diagnostics; SP: survival::cox.zph; IV/DV: fitted Cox coefficients / residual-time association; MF: cox.zph(fit, transform = "rank"); M: assess whether fitted log-HRs remain constant over attained age and identify terms requiring extended modeling.

# 6. Supplementary age-varying coefficient models for PH-flagged main exposures.
# SM: Extended delayed-entry Cox model with time-varying coefficient via tt(); SP: survival::coxph; IV/DV: PH-flagged trauma exposure x centered log attained age / OCD; MF: Surv(...) ~ exposure + tt(exposure) + other covariates + strata(sex); M: characterize non-proportional trauma effects rather than relying only on a time-averaged HR.

# 7. Age-specific HRs with delta-method 95% confidence intervals.
# SM: Linear contrasts from extended Cox coefficients with delta-method variance; SP: survival::coxph + base R matrix calculations; IV/DV: attained age / age-specific trauma HR; MF: logHR(age) = beta_exposure + beta_tt × [log(age) - center]; M: translate time-varying coefficients into interpretable HRs at specific ages.

# 8. No RQ3 sex-association model and no sex-split trauma models.
# SM: Not applicable; SP: not applicable; IV/DV: sex is not a new RQ3 exposure; MF: none; M: RQ2 owns the sex association and sex PH diagnostic, while RQ3 uses the inherited global sex handling.

# 9. Age-at-most-recent/worst-trauma Cox models using the inherited global sex handling.
# SM: Delayed-entry attained-age Cox PH with the one-time sex decision applied; SP: survival::coxph; IV/DV: age at most recent or worst trauma / OCD; MF: Surv(...) ~ age_at_trauma + strata(sex) for the current data; M: test whether age at trauma predicts OCD hazard while consistently accounting for sex non-proportionality.

# 10. Race and birth-country Cox models using the same global sex handling when prepared variables exist.
# SM: Delayed-entry attained-age Cox models for categorical demographic predictors with sex adjustment/stratification; SP: survival::coxph; IV/DV: race or birth country + sex handling / OCD; MF: Surv(...) ~ factor(race/birth_country) + sex handling; M: evaluate demographic associations when the required prepared variables and sufficient category support are available.

# 11. Omnibus likelihood-ratio and joint Wald tests.
# SM: Global hypothesis tests for multi-parameter predictor sets; SP: survival::coxph + stats; IV/DV: sets of trauma/demographic coefficients / model fit or joint coefficient null; MF: full versus reduced Cox model LRT and joint Wald test H0: beta1 = ... = betak = 0; M: test whether a predictor set contributes collectively rather than interpreting only individual coefficients.

# 12. VIF, pairwise-correlation, and condition-number multicollinearity diagnostics.
# SM: Multicollinearity diagnostics on the predictor design matrix; SP: base R / matrix diagnostics; IV/DV: mutually adjusted trauma predictors / collinearity measures; MF: VIF, pairwise correlations, and condition number of X; M: assess whether correlated trauma exposures make coefficient estimates unstable or difficult to interpret.

# 13. Benjamini-Hochberg FDR 
# SM: Multiple-testing adjustment; SP: stats::p.adjust(method = "BH"); IV/DV: nominal p-values / FDR-adjusted p-values; MF: p.adjust(p, method = "BH"); M: control the expected false-discovery proportion across the prespecified family of substantive trauma tests.
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

progress_log("RQ3", sprintf("helper API %d loaded from %s", SURVIVAL_HELPER_API, normalizePath(helper_path, mustWork = FALSE)))

INVENTORY <- c(
  "Trauma-variable eligibility and sparsity screening",
  "Inherited RQ2 sex PH decision and global sex handling",
  "Main age-at-latest-trauma Cox model using inherited sex handling",
  "Main age-at-worst-trauma Cox model using inherited sex handling",
  "Main race Cox model using inherited sex handling and omnibus LRT",
  "Main birth-country Cox model using inherited sex handling and omnibus LRT",
  "Main trauma-type-count Cox model using inherited sex handling",
  "Main mutually adjusted trauma-type Cox model using inherited sex handling",
  "Raw cox.zph PH diagnostics and scaled Schoenfeld residual plots for substantive RQ3 models",
  "Age-varying coefficient models for PH-flagged main exposures",
  "Age-specific hazard ratios with delta-method confidence intervals",
  "Omnibus likelihood-ratio and joint Wald tests",
  "VIF and multicollinearity diagnostics for the mutually adjusted trauma-type model",
  "BH-FDR for the prespecified primary trauma-type coefficient family only",
  "Figure-ready trauma HR table for Python/Matplotlib rendering"
)
analysis_inventory("RQ3: Risk and resilience factors", INVENTORY)

out_dir <- file.path(OUTPUT_ROOT, "RQ3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ph_plot_dir <- file.path(out_dir, "PH_residual_plots")
dir.create(ph_plot_dir, recursive = TRUE, showWarnings = FALSE)

progress_log("RQ3", "loading and validating the prepared dataset")
data <- read_survival_data()
progress_log("RQ3", sprintf("loaded N=%d participants; events=%d", nrow(data), sum(data$event == 1, na.rm = TRUE)))
require_columns(data, c("record_id", "sex", "sex_female", "event", "age_entry", "age_exit"), "RQ3")

TRAUMA_VARS <- c(
  "sexual_abuse",
  "physical_abuse",
  "emotional_abuse",
  "bullying",
  "accident_experienced_witnessed",
  "natural_disaster",
  "illness_injury",
  "kidnapping",
  "family_violence",
  "community_violence"
)

# These transparent defaults should match the final statistical-analysis plan.
MAX_MISSING_PROP <- as.numeric(Sys.getenv("TRAUMA_MAX_MISSING_PROP", "0.20"))
MIN_EXPOSED_N <- as.integer(Sys.getenv("TRAUMA_MIN_EXPOSED_N", "10"))
MIN_UNEXPOSED_N <- as.integer(Sys.getenv("TRAUMA_MIN_UNEXPOSED_N", "10"))
MIN_EVENTS_PER_LEVEL <- as.integer(Sys.getenv("TRAUMA_MIN_EVENTS_PER_LEVEL", "3"))

progress_log("RQ3", "screening trauma variables for missingness, sparsity, and event support")
screening <- screen_binary_exposures(
  data,
  TRAUMA_VARS,
  max_missing_prop = MAX_MISSING_PROP,
  min_exposed_n = MIN_EXPOSED_N,
  min_unexposed_n = MIN_UNEXPOSED_N,
  min_events_per_level = MIN_EVENTS_PER_LEVEL
)
eligible_trauma <- screening$Variable[screening$Eligible]
progress_log("RQ3", sprintf("eligible trauma variables: %d of %d", length(eligible_trauma), length(TRAUMA_VARS)))

term_tables <- list()
model_tables <- list()
ph_tables <- list()
omnibus_tables <- list()
lrt_tables <- list()
vif_tables <- list()
correlation_tables <- list()
multicollinearity_tables <- list()
extended_age_rows <- list()
extended_term_rows <- list()
sex_decision_tables <- list()
ph_plot_audits <- list()
errors <- list()
skipped <- list()

record_error <- function(model, error) {
  errors[[length(errors) + 1L]] <<- data.frame(Model = model, Error = conditionMessage(error), stringsAsFactors = FALSE)
}

record_skip <- function(analysis, reason) {
  skipped[[length(skipped) + 1L]] <<- data.frame(Analysis = analysis, Reason = reason, stringsAsFactors = FALSE)
}
# Adjustment pool used only for targeted extended sensitivity models.
extended_adjustment_pool <- unique(c(psych_available, family_history_available))
run_model <- function(
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
    record_error(model_name, e)
    NULL
  })
}

# Sex handling is inherited from the RQ2 one-time sex association/PH assessment.
# RQ3 does not refit or re-test sex; it only applies the RQ2 decision.
load_rq2_sex_decision <- function(context = "RQ3") {
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

rq2_sex_decision <- load_rq2_sex_decision("RQ3")
SEX_PH_P <- suppressWarnings(as.numeric(rq2_sex_decision$Sex_PH_p[[1L]]))
SEX_STRATIFY <- isTRUE(rq2_sex_decision$Sex_stratified[[1L]]) ||
  tolower(as.character(rq2_sex_decision$Sex_stratified[[1L]])) %in% c("true", "t", "1", "yes")
GLOBAL_SEX_STRATA <- if (SEX_STRATIFY) "sex" else character()
GLOBAL_SEX_COVARIATES <- if (SEX_STRATIFY) character() else "sex"
rq2_sex_decision$Model <- "RQ3 inherited RQ2 sex PH decision"
rq2_sex_decision$Final_formula <- if (SEX_STRATIFY) "RQ3 Cox models: ... + strata(sex)" else "RQ3 Cox models: ... + sex"
rq2_sex_decision$Final_sex_handling <- if (SEX_STRATIFY) "strata(sex) inherited from RQ2" else "sex covariate inherited from RQ2"
sex_decision_tables[[1L]] <- rq2_sex_decision
progress_log("RQ3 sex handling", paste0("inherited RQ2 sex PH p ", fmt_p(SEX_PH_P), "; RQ3 models use ", if (SEX_STRATIFY) "strata(sex)" else "sex as a covariate"))

# -----------------------------------------------------------------------------
# Primary trauma-type model
# -----------------------------------------------------------------------------
trauma_model <- NULL
if (length(eligible_trauma)) {
  progress_log("RQ3", "fitting mutually adjusted trauma-type model")
  trauma_model <- run_model(
    "RQ3 primary trauma-type model",
    covariates = c(eligible_trauma, GLOBAL_SEX_COVARIATES),
    strata = GLOBAL_SEX_STRATA,
    main_exposures = eligible_trauma,
    nuisance = GLOBAL_SEX_COVARIATES
  )
  if (!is.null(trauma_model)) {
    omnibus <- wald_joint_test(trauma_model$fit, eligible_trauma, "Primary trauma-type joint Wald test")
    omnibus$Model <- "RQ3 primary trauma-type model"
    omnibus_tables[[length(omnibus_tables) + 1L]] <- omnibus
    lrt <- likelihood_ratio_block_test(
      trauma_model$fit,
      eligible_trauma,
      "Primary trauma-type omnibus likelihood-ratio test"
    )
    lrt$Model <- "RQ3 primary trauma-type model"
    lrt_tables[[length(lrt_tables) + 1L]] <- lrt

    fit_data <- attr(trauma_model$fit, "analysis_data")
    vif <- numeric_vif_table(fit_data, eligible_trauma, "RQ3 primary trauma-type model")
    corr <- predictor_correlation_long(fit_data, eligible_trauma, "RQ3 primary trauma-type model")
    multi <- multicollinearity_summary(fit_data, eligible_trauma, vif, "RQ3 primary trauma-type model")
    vif_tables[[length(vif_tables) + 1L]] <- vif
    correlation_tables[[length(correlation_tables) + 1L]] <- corr
    multicollinearity_tables[[length(multicollinearity_tables) + 1L]] <- multi
  }
} else {
  record_skip("Primary trauma-type model", "No trauma variables passed the prespecified screening thresholds")
}

# -----------------------------------------------------------------------------
# Trauma-type-count model
# -----------------------------------------------------------------------------
count_model <- NULL
if ("trauma_type_count" %in% names(data)) {
  progress_log("RQ3", "fitting trauma-type-count model")
  count_model <- run_model(
    "RQ3 trauma-type-count model",
    covariates = c("trauma_type_count", GLOBAL_SEX_COVARIATES),
    strata = GLOBAL_SEX_STRATA,
    main_exposures = "trauma_type_count",
    nuisance = GLOBAL_SEX_COVARIATES
  )
  if (!is.null(count_model)) {
    omnibus <- wald_joint_test(count_model$fit, "trauma_type_count", "Trauma-type-count Wald test")
    omnibus$Model <- "RQ3 trauma-type-count model"
    omnibus_tables[[length(omnibus_tables) + 1L]] <- omnibus
    lrt <- likelihood_ratio_block_test(
      count_model$fit,
      "trauma_type_count",
      "Trauma-type-count likelihood-ratio test"
    )
    lrt$Model <- "RQ3 trauma-type-count model"
    lrt_tables[[length(lrt_tables) + 1L]] <- lrt
  }
} else {
  record_skip("Trauma-type-count model", "trauma_type_count is absent from the prepared dataset")
}

# -----------------------------------------------------------------------------
# Supplementary extended models for PH-flagged main trauma exposures
# -----------------------------------------------------------------------------
if (!is.null(trauma_model)) {
  final_ph <- trauma_model$final_ph %||% trauma_model$ph
  global_ph_concern <- any(final_ph$PH_term == "GLOBAL" & final_ph$PH_concern %in% TRUE)
  flagged <- if (global_ph_concern) {
    intersect(final_ph$PH_term[final_ph$PH_concern %in% TRUE], eligible_trauma)
  } else {
    character()
  }
  if (!global_ph_concern) {
    progress_log("RQ3 PH hierarchy", "primary trauma-type global PH test not concerning; no automatic term-specific time-varying models fitted")
  }
  for (exposure in flagged) {
    model_name <- paste("RQ3 age-varying trauma effect:", exposure)
    adjustment <- setdiff(eligible_trauma, exposure)
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
      extended_age_rows[[length(extended_age_rows) + 1L]] <- ext$age_specific
      extended_term_rows[[length(extended_term_rows) + 1L]] <- ext$interaction
      model_tables[[length(model_tables) + 1L]] <- extract_model_stats(ext$fit, model_name)
    }
  }
}

if (!is.null(count_model)) {
  count_final_ph <- count_model$final_ph %||% count_model$ph
  count_flag <- any(count_final_ph$PH_term == "trauma_type_count" & count_final_ph$PH_concern %in% TRUE)
  if (count_flag) {
    model_name <- "RQ3 age-varying trauma-type-count effect"
    ext <- tryCatch(
      fit_extended_age_effect(
        data,
        "trauma_type_count",
        adjustment_terms = GLOBAL_SEX_COVARIATES,
        strata = GLOBAL_SEX_STRATA,
        model_name = model_name
      ),
      error = function(e) { record_error(model_name, e); NULL }
    )
    if (!is.null(ext)) {
      extended_age_rows[[length(extended_age_rows) + 1L]] <- ext$age_specific
      extended_term_rows[[length(extended_term_rows) + 1L]] <- ext$interaction
      model_tables[[length(model_tables) + 1L]] <- extract_model_stats(ext$fit, model_name)
    }
  }
}

# -----------------------------------------------------------------------------
# Age at trauma models
# -----------------------------------------------------------------------------
for (exposure in c("age_at_most_recent_trauma", "age_at_worst_trauma")) {
  if (!(exposure %in% names(data))) {
    record_skip(paste("Model for", exposure), "Variable absent")
    next
  }
  model_name <- paste("RQ3", gsub("_", " ", exposure))
  obj <- run_model(model_name, c(exposure, GLOBAL_SEX_COVARIATES), strata = GLOBAL_SEX_STRATA, main_exposures = exposure, nuisance = GLOBAL_SEX_COVARIATES)
  final_ph <- if (!is.null(obj)) obj$final_ph %||% obj$ph else data.frame()
  if (!is.null(obj) && any(final_ph$PH_term == exposure & final_ph$PH_concern %in% TRUE)) {
    ext_name <- paste("RQ3 age-varying effect:", exposure)
    ext <- tryCatch(
      fit_extended_age_effect(
        data,
        exposure,
        adjustment_terms = GLOBAL_SEX_COVARIATES,
        strata = GLOBAL_SEX_STRATA,
        model_name = ext_name
      ),
      error = function(e) { record_error(ext_name, e); NULL }
    )
    if (!is.null(ext)) {
      extended_age_rows[[length(extended_age_rows) + 1L]] <- ext$age_specific
      extended_term_rows[[length(extended_term_rows) + 1L]] <- ext$interaction
      model_tables[[length(model_tables) + 1L]] <- extract_model_stats(ext$fit, ext_name)
    }
  }
}

# -----------------------------------------------------------------------------
# Conditional race and birth-country models
# -----------------------------------------------------------------------------
run_categorical_model2 <- function(candidates, label) {
  variable <- first_existing(data, candidates)
  if (is.na(variable)) {
    record_skip(label, paste("None of the candidate columns were present:", paste(candidates, collapse = ", ")))
    return(invisible(NULL))
  }
  data[[variable]] <<- droplevels(factor(data[[variable]]))
  if (nlevels(data[[variable]]) < 2L) {
    record_skip(label, paste(variable, "has fewer than two observed levels"))
    return(invisible(NULL))
  }
  model_name <- paste("RQ3", label)
  obj <- run_model(model_name, c(variable, GLOBAL_SEX_COVARIATES), strata = GLOBAL_SEX_STRATA, main_exposures = variable, nuisance = GLOBAL_SEX_COVARIATES)
  if (!is.null(obj)) {
    coefficient_names <- names(coef(obj$fit))[grepl(paste0("^", variable), names(coef(obj$fit)))]
    omnibus <- wald_joint_test(obj$fit, coefficient_names, paste(label, "joint Wald test"))
    omnibus$Model <- model_name
    omnibus_tables[[length(omnibus_tables) + 1L]] <<- omnibus
    lrt <- likelihood_ratio_block_test(obj$fit, variable, paste(label, "omnibus likelihood-ratio test"))
    lrt$Model <- model_name
    lrt_tables[[length(lrt_tables) + 1L]] <<- lrt
  }
  invisible(obj)
}
run_categorical_model <- function(candidates, label) {
  
  variable <- first_existing(data, candidates)
  
  if (is.na(variable)) {
    record_skip(
      label,
      paste(
        "None of the candidate columns were present:",
        paste(candidates, collapse = ", ")
      )
    )
    return(invisible(NULL))
  }
  
  # Convert to factor
  data[[variable]] <<- droplevels(factor(data[[variable]]))
  
  # ----------------------------------------------------------
  # Automatically use largest observed category as reference
  # ----------------------------------------------------------
  category_counts <- table(data[[variable]], useNA = "no")
  
  reference_level <- names(category_counts)[
    which.max(category_counts)
  ]
  
  data[[variable]] <<- relevel(
    data[[variable]],
    ref = reference_level
  )
  
  # Print reference and category counts
  cat(
    "\n",
    label,
    ": reference category = ",
    reference_level,
    "\n",
    sep = ""
  )
  
  print(
    sort(category_counts, decreasing = TRUE)
  )
  
  # Require at least two observed categories
  if (nlevels(data[[variable]]) < 2L) {
    record_skip(
      label,
      paste(variable, "has fewer than two observed levels")
    )
    return(invisible(NULL))
  }
  
  model_name <- paste("RQ3", label)
  
  obj <- run_model(
    model_name,
    c(variable, GLOBAL_SEX_COVARIATES),
    strata = GLOBAL_SEX_STRATA,
    main_exposures = variable,
    nuisance = GLOBAL_SEX_COVARIATES
  )
  
  if (!is.null(obj)) {
    
    coefficient_names <- names(coef(obj$fit))[
      grepl(
        paste0("^", variable),
        names(coef(obj$fit))
      )
    ]
    
    omnibus <- wald_joint_test(
      obj$fit,
      coefficient_names,
      paste(label, "joint Wald test")
    )
    
    omnibus$Model <- model_name
    omnibus_tables[[length(omnibus_tables) + 1L]] <<- omnibus
    
    lrt <- likelihood_ratio_block_test(
      obj$fit,
      variable,
      paste(label, "omnibus likelihood-ratio test")
    )
    
    lrt$Model <- model_name
    lrt_tables[[length(lrt_tables) + 1L]] <<- lrt
  }
  
  invisible(obj)
}

run_categorical_model <- function(candidates, label) {
  
  variable <- first_existing(data, candidates)
  
  if (is.na(variable)) {
    record_skip(
      label,
      paste(
        "None of the candidate columns were present:",
        paste(candidates, collapse = ", ")
      )
    )
    return(invisible(NULL))
  }
  
  # ----------------------------------------------------------
  # Clean categorical variable
  # ----------------------------------------------------------
  x <- droplevels(factor(data[[variable]]))
  # Exclude Canada from birth-country model because of extreme sparsity
  # (N = 3; 0 OCD events)
  if (variable == "birth_country") {
    x[x == "Canada"] <- NA
    x <- droplevels(x)
  }
  
  # Count participants in each observed category
  category_counts <- table(x, useNA = "no")
  
  # Largest category
  reference_level <- names(category_counts)[which.max(category_counts)]
  
  # ----------------------------------------------------------
  # Rebuild factor with largest category FIRST
  # This guarantees treatment coding uses it as reference
  # ----------------------------------------------------------
  new_levels <- c(
    reference_level,
    setdiff(levels(x), reference_level)
  )
  
  data[[variable]] <<- factor(
    as.character(x),
    levels = new_levels
  )
  
  # ----------------------------------------------------------
  # Verify what will actually be used
  # ----------------------------------------------------------
  cat("\n========================================\n")
  cat("CATEGORICAL MODEL:", label, "\n")
  cat("Variable:", variable, "\n")
  cat("========================================\n")
  
  cat("\nCategory counts:\n")
  print(sort(category_counts, decreasing = TRUE))
  
  cat("\nLargest category:", reference_level, "\n")
  cat("Factor reference:", levels(data[[variable]])[1], "\n")
  
  cat("\nFactor levels in model order:\n")
  print(levels(data[[variable]]))
  
  cat("\nTreatment contrasts:\n")
  print(contrasts(data[[variable]]))
  
  # ----------------------------------------------------------
  # Require at least two observed categories
  # ----------------------------------------------------------
  if (nlevels(data[[variable]]) < 2L) {
    record_skip(
      label,
      paste(variable, "has fewer than two observed levels")
    )
    return(invisible(NULL))
  }
  
  # ----------------------------------------------------------
  # Fit model
  # ----------------------------------------------------------
  model_name <- paste("RQ3", label)
  
  obj <- run_model(
    model_name,
    c(variable, GLOBAL_SEX_COVARIATES),
    strata = GLOBAL_SEX_STRATA,
    main_exposures = variable,
    nuisance = GLOBAL_SEX_COVARIATES
  )
  
  # ----------------------------------------------------------
  # Omnibus tests
  # ----------------------------------------------------------
  if (!is.null(obj)) {
    
    coefficient_names <- names(coef(obj$fit))[
      grepl(
        paste0("^", variable),
        names(coef(obj$fit))
      )
    ]
    
    # Print fitted coefficient names
    cat("\nCoefficients actually estimated:\n")
    print(coefficient_names)
    
    omnibus <- wald_joint_test(
      obj$fit,
      coefficient_names,
      paste(label, "joint Wald test")
    )
    
    omnibus$Model <- model_name
    omnibus_tables[[length(omnibus_tables) + 1L]] <<- omnibus
    
    lrt <- likelihood_ratio_block_test(
      obj$fit,
      variable,
      paste(label, "omnibus likelihood-ratio test")
    )
    
    lrt$Model <- model_name
    lrt_tables[[length(lrt_tables) + 1L]] <<- lrt
  }
  
  invisible(obj)
}

run_categorical_model(c("race_category", "race", "child_race", "cfhx_race"), "race model")
run_categorical_model(c("birth_country", "country_of_birth", "child_birth_country", "cfhx_birth_country"), "birth-country model")

# -----------------------------------------------------------------------------
# Final tables and substantive FDR
# -----------------------------------------------------------------------------
coefficient_results <- rbind_fill(c(term_tables, extended_term_rows))
model_results <- rbind_fill(model_tables)
ph_results <- rbind_fill(ph_tables)
omnibus_results <- rbind_fill(omnibus_tables)
lrt_results <- rbind_fill(lrt_tables)
vif_results <- rbind_fill(vif_tables)
correlation_results <- rbind_fill(correlation_tables)
multicollinearity_results <- rbind_fill(multicollinearity_tables)
extended_results <- rbind_fill(extended_age_rows)
error_results <- rbind_fill(errors)
skipped_results <- rbind_fill(skipped)
sex_decisions <- rbind_fill(sex_decision_tables)
ph_plot_audit <- rbind_fill(ph_plot_audits)
write_csv_safe(ph_plot_audit, file.path(out_dir, "RQ3_PH_residual_plot_audit.csv"))

# RQ3 inherits the sex PH diagnostic from RQ2. RQ3 PH plots are saved only for
# substantive RQ3 predictors; no RQ3 sex residual plot is generated or aliased.
inherited_sex_diagnostic_note <- data.frame(
  Model = "RQ3 inherited RQ2 sex PH diagnostic",
  RQ2_decision_file = rq2_sex_decision$Decision_source_file[[1L]],
  Sex_PH_p = SEX_PH_P,
  Sex_handling = if (SEX_STRATIFY) "strata(sex)" else "sex covariate",
  RQ3_refit_sex_model = FALSE,
  stringsAsFactors = FALSE
)
write_csv_safe(inherited_sex_diagnostic_note, file.path(out_dir, "RQ3_inherited_RQ2_sex_diagnostic.csv"))
progress_log("RQ3 sex handling", "sex diagnostic inherited from RQ2; RQ3 does not save a new sex Schoenfeld plot")

progress_log("RQ3", "saving Figure 3 data for Python/Matplotlib rendering")
trauma_display_labels <- c(
  sexual_abuse="Sexual abuse", physical_abuse="Physical abuse", emotional_abuse="Emotional abuse",
  bullying="Bullying", accident_experienced_witnessed="Accident experienced/witnessed",
  natural_disaster="Natural disaster", illness_injury="Illness/injury", kidnapping="Kidnapping",
  family_violence="Family violence", community_violence="Community violence",
  trauma_type_count="Number of trauma types"
)
ph_lookup_for_figure <- function(model, term) {
  required <- c("Model", "PH_term", "PH_concern")
  if (!is.data.frame(ph_results) || !all(required %in% names(ph_results))) {
    return(FALSE)
  }
  
  x <- ph_results[
    ph_results$Model == model & ph_results$PH_term == term,
    ,
    drop = FALSE
  ]
  
  if (!nrow(x)) {
    return(FALSE)
  }
  
  isTRUE(x$PH_concern[[nrow(x)]])
}

type_rows_fig <- data.frame()
if (
  is.data.frame(coefficient_results) &&
  all(c("Model", "Term", "HR", "CI_lower", "CI_upper") %in% names(coefficient_results))
) {
  type_rows_fig <- coefficient_results[
    coefficient_results$Model == "RQ3 primary trauma-type model" &
      coefficient_results$Term %in% eligible_trauma,
    ,
    drop = FALSE
  ]
  
  type_rows_fig <- type_rows_fig[
    is.finite(type_rows_fig$HR) &
      is.finite(type_rows_fig$CI_lower) &
      is.finite(type_rows_fig$CI_upper),
    ,
    drop = FALSE
  ]
}
type_rows_fig <- type_rows_fig[is.finite(type_rows_fig$HR) & is.finite(type_rows_fig$CI_lower) & is.finite(type_rows_fig$CI_upper),,drop=FALSE]
if (nrow(type_rows_fig)) {
  type_rows_fig$PH_concern <- vapply(type_rows_fig$Term, function(z) ph_lookup_for_figure("RQ3 primary trauma-type model",z), logical(1))
  type_rows_fig <- type_rows_fig[order(type_rows_fig$HR),,drop=FALSE]
}
count_rows_fig <- coefficient_results[coefficient_results$Model == "RQ3 trauma-type-count model" & coefficient_results$Term == "trauma_type_count",,drop=FALSE]
count_rows_fig <- count_rows_fig[is.finite(count_rows_fig$HR) & is.finite(count_rows_fig$CI_lower) & is.finite(count_rows_fig$CI_upper),,drop=FALSE]
if (nrow(count_rows_fig)) count_rows_fig$PH_concern <- ph_lookup_for_figure("RQ3 trauma-type-count model","trauma_type_count")
figure3_data <- rbind_fill(list(type_rows_fig,count_rows_fig))
if (nrow(figure3_data)) {
  figure3_data$Display <- unname(trauma_display_labels[figure3_data$Term])
  figure3_data$Display[is.na(figure3_data$Display)] <- gsub("_"," ",figure3_data$Term[is.na(figure3_data$Display)])
  figure3_data$Significance <- p_stars(figure3_data$p)
  figure3_data$Is_count <- figure3_data$Term == "trauma_type_count"
  figure3_data$Plot_order <- seq_len(nrow(figure3_data))
  figure3_data$Plot_color <- c(rep(SET2_PALETTE, length.out=max(0,nrow(figure3_data)-sum(figure3_data$Is_count))), rep("black",sum(figure3_data$Is_count)))
}
write_csv_safe(figure3_data, file.path(out_dir,"Figure3_trauma_predictors_data.csv"))
write_csv_safe(figure3_data, file.path(out_dir,"RQ3_trauma_HR_figure_data.csv"))

rq3_model_specifications <- data.frame(
  Analysis = c(
    "Inherited RQ2 sex PH decision",
    "Main age-at-latest/worst-trauma and demographic models",
    "Primary trauma-count model",
    "Primary mutually adjusted trauma-type model",
    "Age-varying coefficient extension",
    "PH diagnostic"
  ),
  R_package = c("base R", "survival", "survival", "survival", "survival", "survival"),
  R_function = c("read.csv", "survival::coxph", "survival::coxph", "survival::coxph", "survival::coxph", "survival::cox.zph"),
  Model_type = c(
    "Imported RQ2 one-time sex association/PH decision; no RQ3 sex model fitted",
    "Delayed-entry attained-age Cox with inherited sex handling",
    "Delayed-entry attained-age Cox with inherited sex handling",
    "Delayed-entry attained-age Cox with inherited sex handling",
    "Delayed-entry attained-age Cox with exposure x centered log attained-age term",
    "Scaled Schoenfeld residual PH test"
  ),
  Formula = c(
    "RQ2: Surv(age_entry, age_exit, event) ~ sex; RQ3 inherits the resulting sex handling",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ predictor + strata(sex)" else "Surv(age_entry, age_exit, event) ~ predictor + sex",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ trauma_type_count + strata(sex)" else "Surv(age_entry, age_exit, event) ~ trauma_type_count + sex",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ eligible trauma indicators + strata(sex)" else "Surv(age_entry, age_exit, event) ~ eligible trauma indicators + sex",
    if (SEX_STRATIFY) "Surv(age_entry, age_exit, event) ~ exposure + tt(exposure) + strata(sex)" else "Surv(age_entry, age_exit, event) ~ exposure + tt(exposure) + sex",
    "cox.zph(fit, transform='rank')"
  ),
  Time_scale = c(NA, rep("Attained age (years)", 4), "Rank-transformed event time for diagnostic"),
  Ties = c(NA, rep("Efron", 4), NA),
  Variance = c(NA, rep("Model-based", 4), NA),
  Figure_renderer = "Python / Matplotlib from saved CSV; no manuscript Figure 3 is drawn in R",
  stringsAsFactors = FALSE
)
write_csv_safe(rq3_model_specifications, file.path(out_dir, "RQ3_model_specifications.csv"))


sex_decision_lines <- character()
if (nrow(sex_decisions)) {
  sex_decision_lines <- apply(sex_decisions, 1, function(x) sprintf(
    "%s: sex PH test chi-square(%s) = %s, p %s; global handling: %s.",
    x[["Model"]], fmt_num(as.numeric(x[["Sex_PH_df"]]), 0),
    fmt_num(as.numeric(x[["Sex_PH_chisq"]])), fmt_p(as.numeric(x[["Sex_PH_p"]])),
    x[["Final_sex_handling"]]
  ))
}

screening_lines <- sprintf(
  "%s: available N = %s, exposed N = %s, unexposed N = %s, exposed events = %s, unexposed events = %s, missing = %s%%; eligible = %s (%s).",
  screening$Variable, screening$Available_N, screening$Exposed_N, screening$Unexposed_N,
  screening$Exposed_events, screening$Unexposed_events,
  fmt_num(100 * screening$Missing_prop), ifelse(screening$Eligible, "yes", "no"), screening$Reason
)

extended_lines <- character()
if (nrow(extended_results)) {
  extended_lines <- apply(extended_results, 1, function(x) sprintf(
    "%s at age %s: HR = %s, 95%% CI [%s, %s], z = %s, p %s.",
    x[["Model"]], fmt_num(as.numeric(x[["Age"]]), 1), fmt_num(as.numeric(x[["HR"]])),
    fmt_num(as.numeric(x[["CI_lower"]])), fmt_num(as.numeric(x[["CI_upper"]])),
    fmt_num(as.numeric(x[["z"]])), fmt_p(as.numeric(x[["p"]]))
  ))
}

lrt_lines <- character()
if (nrow(lrt_results)) {
  lrt_lines <- apply(lrt_results, 1, function(x) sprintf(
    "%s: omnibus likelihood-ratio chi-square(%s) = %s, p %s.",
    x[["Model"]], fmt_num(as.numeric(x[["df"]]), 0),
    fmt_num(as.numeric(x[["LR_chisq"]])), fmt_p(as.numeric(x[["p"]]))
  ))
}

multicollinearity_lines <- character()
if (nrow(multicollinearity_results)) {
  multicollinearity_lines <- apply(multicollinearity_results, 1, function(x) sprintf(
    "%s: maximum VIF = %s, mean VIF = %s, maximum absolute pairwise correlation = %s, condition number = %s; any VIF > 5: %s; any VIF > 10: %s.",
    x[["Model"]], fmt_num(as.numeric(x[["Maximum_VIF"]])),
    fmt_num(as.numeric(x[["Mean_VIF"]])),
    fmt_num(as.numeric(x[["Maximum_absolute_pairwise_correlation"]])),
    fmt_num(as.numeric(x[["Condition_number"]])),
    x[["Any_VIF_gt_5"]], x[["Any_VIF_gt_10"]]
  ))
}

save_apa_bundle(
  output_dir = out_dir,
  prefix = "RQ3",
  title = "RQ3: Risk and resilience factors",
  inventory = INVENTORY,
  model_terms = coefficient_results,
  model_stats = model_results,
  ph = ph_results,
  extra_tables = list(
    trauma_screening = screening,
    omnibus_Wald_tests = omnibus_results,
    omnibus_likelihood_ratio_tests = lrt_results,
    trauma_type_VIF = vif_results,
    trauma_type_pairwise_correlations = correlation_results,
    trauma_type_multicollinearity_summary = multicollinearity_results,
    trauma_HR_figure_data = figure3_data,
    extended_age_specific_HRs = extended_results,
    inherited_RQ2_sex_PH_decision = sex_decisions,
    inherited_RQ2_sex_diagnostic_note = inherited_sex_diagnostic_note,
    PH_residual_plot_audit = ph_plot_audit,
    skipped_analyses = skipped_results,
    model_errors = error_results,
    model_specifications = rq3_model_specifications
  ),
  additional_lines = c(
    "SCALED SCHOENFELD RESIDUAL PLOTS",
    paste0("Residual plots were saved under ", normalizePath(ph_plot_dir, mustWork = FALSE), ". Sex residual diagnostics are inherited from RQ2; RQ3 saves plots only for substantive RQ3 predictors."),
    "",
    "INHERITED RQ2 SEX PH DECISION AND GLOBAL HANDLING",
    if (length(sex_decision_lines)) sex_decision_lines else "The RQ2 sex PH decision was not available.",
    paste0("RQ3 imports the RQ2 sex decision and does not fit a new sex-association model. Current handling: ", if (SEX_STRATIFY) "strata(sex), allowing sex-specific baseline hazards while estimating common exposure HRs." else "sex retained as a covariate."),
    "No female-only or male-only trauma models are fitted.",
    "",
    "",
    "OMNIBUS LIKELIHOOD-RATIO TESTS",
    if (length(lrt_lines)) lrt_lines else "No omnibus likelihood-ratio test was estimable.",
    "",
    "MULTICOLLINEARITY DIAGNOSTICS",
    if (length(multicollinearity_lines)) multicollinearity_lines else "Multicollinearity diagnostics were not estimable.",
    "",
    "TRAUMA SCREENING",
    screening_lines,
    "",
    "EXTENDED-MODEL AGE-SPECIFIC RESULTS",
    if (length(extended_lines)) extended_lines else "No main exposure was PH-flagged or no extended model was estimable.",
    "",
    paste0(
      "Screening thresholds: maximum missing proportion = ", MAX_MISSING_PROP,
      "; minimum exposed N = ", MIN_EXPOSED_N,
      "; minimum unexposed N = ", MIN_UNEXPOSED_N,
      "; minimum events per exposure level = ", MIN_EVENTS_PER_LEVEL, "."
    ),
    "Raw PH p-values are reported without FDR adjustment. BH-FDR is applied only to substantive coefficient tests. Sex is tested for PH once in RQ2 and the resulting handling is inherited throughout RQ3. For the multivariable trauma-type model, automatic time-varying coefficients require both a concerning global cox.zph test and a concerning term-specific test; residual plots are saved for visual review. Primary proportional Cox estimates remain in the output even when a supplementary age-varying model is fitted."
  )
)

progress_log("RQ3 complete", paste0("outputs: ", normalizePath(out_dir, mustWork = FALSE), "; elapsed=", format_elapsed(SCRIPT_STARTED)))


# =============================================================================
# BH-FDR: PRIMARY TRAUMA-TYPE COEFFICIENTS
# =============================================================================

fdr_results <- coefficient_results[
  coefficient_results$Model == "RQ3 primary trauma-type model" &
    coefficient_results$Term %in% TRAUMA_VARS,
  c("Term", "p"),
  drop = FALSE
]

if (nrow(fdr_results) > 0L) {

  fdr_results$p_BH_FDR <- stats::p.adjust(
    fdr_results$p,
    method = "BH"
  )

  cat("\n")
  cat("============================================================\n")
  cat("BH-FDR CORRECTION: PRIMARY TRAUMA-TYPE COEFFICIENTS\n")
  cat("============================================================\n")
  cat("Number of tests:", nrow(fdr_results), "\n\n")

  print(
    fdr_results,
    row.names = FALSE,
    digits = 4
  )

} else {
  cat("\nNo primary trauma-type coefficients found for BH-FDR correction.\n")
}

# =============================================================================
# RQ3 SENSITIVITY ANALYSES (prespecified)
# =============================================================================

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
# 2. Baseline psychiatric-comorbidity adjustment.
# SM: Covariate-adjusted delayed-entry Cox sensitivity models; SP: survival::coxph; IV/DV: trauma predictors + baseline MINI psychiatric categories / OCD; MF: Surv(age_entry, age_exit, event) ~ trauma predictors + baseline_psychiatric_categories + inherited sex handling; M: assess whether trauma associations are robust to baseline psychiatric comorbidity.

# 3. First-degree-family psychiatric-history adjustment.
# SM: Covariate-adjusted delayed-entry Cox sensitivity models; SP: survival::coxph; IV/DV: trauma predictors + first-degree psychiatric history / OCD; MF: Surv(age_entry, age_exit, event) ~ trauma predictors + first_degree_psychiatric_history_any + inherited sex handling; M: assess whether trauma associations are robust to parent-or-sibling psychiatric history.
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
helper_path <- file.path(script_dir, "custom_functions.R")
if (!file.exists(helper_path)) {
  stop("Cannot find custom_functions.R beside this script: ", helper_path, call. = FALSE)
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
  "Baseline MINI psychiatric-comorbidity adjustment",
  "First-degree-family psychiatric-history adjustment",
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
# 2. Baseline psychiatric-comorbidity adjustment.
psych_available <- usable_adjustments(data, BASELINE_PSYCH)

if (length(psych_available)) {
  run_sensitivity_pair(
    "Psychiatric-comorbidity adjusted",
    data,
    adjustments = psych_available
  )
} else {
  record_skip(
    "Psychiatric-comorbidity adjusted",
    "No usable baseline psychiatric-comorbidity adjustment variables were available"
  )
}

# 3. First-degree-family psychiatric-history adjustment.
family_history_adjustment <- if (
  !is.na(FIRST_DEGREE_HISTORY_VAR) &&
  nzchar(FIRST_DEGREE_HISTORY_VAR)
) {
  FIRST_DEGREE_HISTORY_VAR
} else {
  character()
}

family_history_available <- usable_adjustments(data, family_history_adjustment)

if (length(family_history_available)) {
  run_sensitivity_pair(
    "First-degree-family-history adjusted",
    data,
    adjustments = family_history_available
  )
} else {
  record_skip(
    "First-degree-family-history adjusted",
    "No usable first-degree-family psychiatric-history adjustment variable was available"
  )
}

# 3. Site-clustered robust standard errors only; no site-fixed-effects model.
site_var <- first_existing(data, c("site", "site_id", "record_dag_name", "dag", "data_access_group"))
site_var <- first_existing(
  data,
  c(
    "site",
    "site_id",
    "study_site",
    "redcap_data_access_group",
    "record_dag_name",
    "dag",
    "data_access_group"
  )
)
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
      
      if (exposure %in% eligible_trauma) {
        adjustment <- c(setdiff(eligible_trauma, exposure), extended_adjustment_pool)
      }
      
      if (exposure == "trauma_type_count") {
        adjustment <- extended_adjustment_pool
      }
      
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


psychiatric_formula_term <- if (length(psych_available)) {
  "baseline psychiatric categories"
} else {
  "baseline psychiatric categories not available"
}

family_formula_term <- if (length(family_history_available)) {
  family_history_available[[1L]]
} else {
  "first-degree family psychiatric history not available"
}

sensitivity_model_specifications <- data.frame(
  Analysis = c(
    "Inherited RQ2 sex PH decision",
    "Early-case exclusion trauma-type sensitivity models",
    "Early-case exclusion trauma-count sensitivity models",
    "Psychiatric-comorbidity adjusted trauma models",
    "First-degree-family-history adjusted trauma models",
    "Site-clustered robust sensitivity where available",
    "Age-varying coefficient extension for repeated PH concern"
  ),
  R_package = c("base R", rep("survival", 6)),
  R_function = c("read.csv", rep("survival::coxph", 6)),
  Model_type = c(
    "Imported RQ2 one-time sex association/PH decision; no sensitivity-specific sex model fitted",
    "Restricted-sample delayed-entry attained-age Cox using the inherited sex decision",
    "Restricted-sample delayed-entry attained-age Cox using the inherited sex decision",
    "Covariate-adjusted delayed-entry attained-age Cox using the inherited sex decision",
    "Covariate-adjusted delayed-entry attained-age Cox using the inherited sex decision",
    "Delayed-entry attained-age Cox with clustered sandwich variance and inherited sex decision",
    "Delayed-entry attained-age Cox with exposure x centered log attained-age term and inherited sex decision"
  ),
  Formula = c(
    "RQ2: Surv(age_entry, age_exit, event) ~ sex; sensitivity scripts inherit the resulting sex handling",
    if (SEX_STRATIFY) {
      "Surv(age_entry, age_exit, event) ~ eligible trauma indicators + strata(sex)"
    } else {
      "Surv(age_entry, age_exit, event) ~ eligible trauma indicators + sex"
    },
    if (SEX_STRATIFY) {
      "Surv(age_entry, age_exit, event) ~ trauma_type_count + strata(sex)"
    } else {
      "Surv(age_entry, age_exit, event) ~ trauma_type_count + sex"
    },
    if (SEX_STRATIFY) {
      paste0(
        "Surv(age_entry, age_exit, event) ~ trauma predictors + ",
        psychiatric_formula_term,
        " + strata(sex)"
      )
    } else {
      paste0(
        "Surv(age_entry, age_exit, event) ~ trauma predictors + ",
        psychiatric_formula_term,
        " + sex"
      )
    },
    if (SEX_STRATIFY) {
      paste0(
        "Surv(age_entry, age_exit, event) ~ trauma predictors + ",
        family_formula_term,
        " + strata(sex)"
      )
    } else {
      paste0(
        "Surv(age_entry, age_exit, event) ~ trauma predictors + ",
        family_formula_term,
        " + sex"
      )
    },
    if (SEX_STRATIFY) {
      "Surv(age_entry, age_exit, event) ~ trauma predictors + strata(sex) + cluster(site)"
    } else {
      "Surv(age_entry, age_exit, event) ~ trauma predictors + sex + cluster(site)"
    },
    if (SEX_STRATIFY) {
      "Surv(age_entry, age_exit, event) ~ exposure + tt(exposure) + strata(sex)"
    } else {
      "Surv(age_entry, age_exit, event) ~ exposure + tt(exposure) + sex"
    }
  ),
  Ties = c(NA, rep("Efron", 6)),
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
  c("Exposure", "Sensitivity_models_with_PH_concern", "Targeted_extended_model")
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
  "These supplementary analyses evaluated whether the primary trauma findings were sensitive to early-case exclusion, psychiatric-comorbidity adjustment, first-degree-family-history adjustment, site clustering when available, multiplicity control, and non-proportional hazards.",
  "",
  "Methods",
  "All fixed-coefficient sensitivity analyses used delayed-entry Cox regression with attained age as the time scale and Efron handling of tied event times. Participants entered the risk set at age at study entry and were censored at their final completed OCD assessment.",
  sex_handling_sentence,
  "Sensitivity analyses excluded OCD cases first detected within one month; separately adjusted for available baseline MINI psychiatric categories and first-degree psychiatric history (parent or sibling); and evaluated site clustering when a prepared site field was available.",
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
  #doc <- add_para(doc, "These supplementary analyses evaluated whether the primary trauma findings were sensitive to early-case exclusion, combined psychiatric/first-degree-family confounder adjustment, site clustering when available, multiplicity control, and non-proportional hazards.")
  doc <- add_para(doc, "These supplementary analyses evaluated whether the primary trauma findings were sensitive to early-case exclusion, psychiatric-comorbidity adjustment, first-degree-family-history adjustment, site clustering when available, multiplicity control, and non-proportional hazards.")
  
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
