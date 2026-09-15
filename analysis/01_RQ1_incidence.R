# 01_RQ1_incidence.R
#
# ANALYSES INCLUDED
# 01_RQ1_incidence.R
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

# 1. Analysis-data audit for each time origin.
# SM: Descriptive eligibility/risk-set audit; SP: base R; IV/DV: time origin / valid survival interval and event status; MF: none; M: document N, events, delayed entry, and exclusions before each KM analysis.

# 2. Delayed-entry Kaplan-Meier cumulative OCD risk from study entry.
# SM: Kaplan-Meier survival estimator; SP: survival::survfit; IV/DV: study-entry time / first detected OCD; MF: Surv(entry, duration, event) ~ 1 with entry = 0; M: estimate cumulative OCD risk from study entry without assuming proportional hazards.

# 3. Delayed-entry Kaplan-Meier cumulative OCD risk from most recent trauma.
# SM: Left-truncated Kaplan-Meier survival estimator; SP: survival::survfit; IV/DV: time since most recent trauma / first detected OCD; MF: Surv(entry, duration, event) ~ 1; M: estimate post-trauma cumulative risk while accounting for observation beginning after trauma.

# 4. Delayed-entry Kaplan-Meier cumulative OCD risk from worst trauma.
# SM: Left-truncated Kaplan-Meier survival estimator; SP: survival::survfit; IV/DV: time since worst trauma / first detected OCD; MF: Surv(entry, duration, event) ~ 1; M: estimate cumulative risk since worst trauma while correctly handling delayed entry.

# 5. Overall and sex-specific cumulative risk at 1 and 2 years.
# SM: Kaplan-Meier predictions overall and within sex groups; SP: survival::survfit; IV/DV: sex/time origin / cumulative first-detected OCD risk; MF: separate Surv(entry, duration, event) ~ 1 fits for Overall, Female, and Male; M: report clinically interpretable 1- and 2-year cumulative risks without estimating a sex HR.

# 6. Pointwise and participant-bootstrap 95% confidence intervals.
# SM: KM pointwise CI plus nonparametric participant bootstrap; SP: survival::survfit + base R bootstrap code; IV/DV: participant resamples / cumulative OCD risk; MF: repeated Surv(entry, duration, event) ~ 1; M: quantify uncertainty while preserving participant-level survival histories.

# 7. Numbers at risk, event counts, full curve data, and time to 10% risk.
# SM: Derived Kaplan-Meier/risk-set summaries; SP: survival + base R; IV/DV: time / risk-set size, cumulative risk, and events; MF: derived from fitted Surv(entry, duration, event) ~ 1 curves; M: provide transparent curve diagnostics and interpretable timing of reaching 10% cumulative risk.

# 8. Crude participant-level proportion ever diagnosed with incident OCD during follow-up.
# SM: Binomial descriptive proportion; SP: base R; IV/DV: participant event indicator / proportion with incident OCD; MF: sum(event == 1) / N; M: provide the observed participant-level incidence proportion separately from time-to-event KM estimates.
#
# RQ1 is nonparametric; proportional-hazards testing is not applicable.

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
assert_helper_api(21L, c("progress_log", "format_elapsed", "read_survival_data", "bootstrap_km", "save_apa_bundle"))
SCRIPT_STARTED <- Sys.time()

# =============================================================================
# RQ1 STATISTICAL ENGINE -- core model code is intentionally defined here
# =============================================================================
# Package: survival
# Estimator: survival::survfit
# Response: survival::Surv(entry, duration, event)
# Formula:  Surv(entry, duration, event) ~ 1
# Delayed entry: yes (entry argument is embedded in the counting-process Surv object)
# CI type: log-log pointwise confidence intervals from survfit
# PH assumption: not applicable; this is nonparametric Kaplan-Meier estimation.
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("RQ1 requires the R package 'survival'.", call. = FALSE)
}
progress_log("RQ1 statistical engine", paste0("survival ", as.character(utils::packageVersion("survival")), "; survival::survfit with delayed entry"))

rq1_make_origin_data <- function(data, origin) {
  if (origin == "Study entry") {
    out <- data.frame(
      record_id = data$record_id, sex = data$sex, event = data$event,
      entry = 0, duration = data$time_since_start / 12
    )
  } else if (origin == "Most recent trauma") {
    out <- data.frame(
      record_id = data$record_id, sex = data$sex, event = data$event,
      entry = data$months_since_most_recent_at_baseline / 12,
      duration = data$time_since_most_recent_trauma / 12
    )
  } else if (origin == "Worst trauma") {
    out <- data.frame(
      record_id = data$record_id, sex = data$sex, event = data$event,
      entry = data$months_since_worst_at_baseline / 12,
      duration = data$time_since_worst_trauma / 12
    )
  } else {
    stop("Unknown RQ1 time origin: ", origin, call. = FALSE)
  }
  out <- out[stats::complete.cases(out[, c("event", "entry", "duration")]), , drop = FALSE]
  out[out$event %in% c(0, 1) & out$entry >= 0 & out$duration > out$entry, , drop = FALSE]
}

rq1_fit_km <- function(data) {
  survival::survfit(
    survival::Surv(entry, duration, event) ~ 1,
    data = data,
    conf.type = "log-log"
  )
}

rq1_risk_at_times <- function(fit, data, times = REPORT_TIMES, origin = "", group = "Overall") {
  sm <- summary(fit, times = times, extend = TRUE)
  data.frame(
    Time_origin = origin, Group = group, Time_years = times,
    Risk = 1 - sm$surv, CI_lower = 1 - sm$upper, CI_upper = 1 - sm$lower,
    N_at_risk = sm$n.risk,
    Events_by_time = vapply(times, function(t) sum(data$event == 1 & data$duration <= t), integer(1)),
    Analysis_N = nrow(data), Total_events = sum(data$event),
    Delayed_entries = sum(data$entry > 0), stringsAsFactors = FALSE
  )
}

rq1_curve_table <- function(fit, origin = "", group = "Overall") {
  data.frame(
    Time_origin = origin, Group = group, Time_years = fit$time,
    Risk = 1 - fit$surv, CI_lower = 1 - fit$upper, CI_upper = 1 - fit$lower,
    N_at_risk = fit$n.risk, N_events_at_time = fit$n.event,
    N_censored_at_time = fit$n.censor, stringsAsFactors = FALSE
  )
}

rq1_numbers_at_risk <- function(data, times = seq(0, 2, by = .5), origin = "", group = "Overall") {
  data.frame(
    Time_origin = origin, Group = group, Time_years = times,
    N_at_risk = vapply(times, function(t) sum(data$entry <= t & data$duration >= t), integer(1)),
    stringsAsFactors = FALSE
  )
}


# Dense plotting grid for the Python renderer. No model is fitted here: values
# are step evaluations of the R survival::survfit object above. Python may apply
# display-only rolling smoothing to these R-derived estimates.
rq1_step_value <- function(time, value, grid, initial) {
  out <- rep(initial, length(grid))
  if (!length(time) || !length(value)) return(out)
  keep <- is.finite(time) & is.finite(value)
  time <- time[keep]
  value <- value[keep]
  if (!length(time)) return(out)
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  for (j in seq_along(grid)) {
    idx <- which(time <= grid[[j]])
    if (length(idx)) out[[j]] <- value[[max(idx)]]
  }
  out
}

rq1_python_rolling_grid <- function(fit, data, origin, group, step = 0.1) {
  max_years <- if (identical(origin, "Worst trauma")) 5 else 2
  grid <- seq(0, max_years, by = step)
  risk <- rq1_step_value(fit$time, 1 - fit$surv, grid, 0)
  ci_lower <- rq1_step_value(fit$time, 1 - fit$upper, grid, 0)
  ci_upper <- rq1_step_value(fit$time, 1 - fit$lower, grid, 0)
  n_at_risk <- vapply(
    grid,
    function(t) sum(data$entry <= t & data$duration >= t),
    integer(1)
  )
  data.frame(
    Time_origin = origin,
    Group = group,
    Time_years = grid,
    Risk = risk,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    N_at_risk = n_at_risk,
    Grid_step_years = step,
    R_estimator = "survival::survfit delayed-entry Kaplan-Meier",
    stringsAsFactors = FALSE
  )
}

rq1_time_to_risk <- function(fit, target = .10, origin = "", group = "Overall") {
  idx <- which((1 - fit$surv) >= target)
  data.frame(
    Time_origin = origin, Group = group, Target_risk = target,
    Time_years = if (length(idx)) fit$time[min(idx)] else NA_real_,
    Threshold_reached = length(idx) > 0L, stringsAsFactors = FALSE
  )
}

INVENTORY <- c(
  "Analysis-data audit by time origin",
  "Delayed-entry cumulative OCD risk from study entry",
  "Delayed-entry cumulative OCD risk from most recent trauma",
  "Delayed-entry cumulative OCD risk from worst trauma",
  "Overall and sex-specific risks at 1 and 2 years",
  "Pointwise and bootstrap 95% confidence intervals",
  "Numbers at risk, event counts, curve data, and time to 10% risk",
  "Crude participant-level proportion ever diagnosed with incident OCD during follow-up",
  "Figure-ready Kaplan-Meier curve table for Python rendering",
  "Dense 0.1-year R-derived KM grid for rolling-style Python rendering"
)
analysis_inventory("RQ1: OCD incidence", INVENTORY)

out_dir <- file.path(OUTPUT_ROOT, "RQ1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

progress_log("RQ1", "loading and validating the prepared dataset")
data <- read_survival_data()
progress_log("RQ1", sprintf("loaded N=%d participants; events=%d", nrow(data), sum(data$event == 1, na.rm = TRUE)))
require_columns(
  data,
  c(
    "record_id", "sex", "event", "time_since_start",
    "time_since_most_recent_trauma", "time_since_worst_trauma",
    "months_since_most_recent_at_baseline", "months_since_worst_at_baseline"
  ),
  "RQ1"
)


baseline_descriptive_variables <- c(
  "sex",
  "event",
  "age_entry",
  "age_exit",
  "time_since_start",
  "time_since_most_recent_trauma",
  "time_since_worst_trauma",
  "months_since_most_recent_at_baseline",
  "months_since_worst_at_baseline"
)
baseline_descriptives <- save_descriptive_table(
  data = data,
  variables = baseline_descriptive_variables,
  output_path = file.path(out_dir, "RQ1_baseline_descriptive_summary.csv"),
  group_var = NULL,
  label = "RQ1 baseline descriptive summary"
)

# Crude observed incident proportion across all completed follow-up assessments.
# This is the participant-level proportion ever observed OCD-positive and does not
# account for differential censoring; the KM estimates below remain the primary risk estimates.
crude_groups <- list(Overall = data, Female = data[data$sex == "Female", , drop = FALSE], Male = data[data$sex == "Male", , drop = FALSE])
crude_incident_rows <- lapply(names(crude_groups), function(group_name) {
  d <- crude_groups[[group_name]]
  d <- d[d$event %in% c(0, 1), , drop = FALSE]
  n <- nrow(d)
  events <- sum(d$event == 1)
  bt <- if (n > 0L) binom.test(events, n) else NULL
  data.frame(
    Group = group_name,
    N = n,
    Incident_OCD_cases = events,
    Crude_proportion = if (n > 0L) events / n else NA_real_,
    Crude_percent = if (n > 0L) 100 * events / n else NA_real_,
    Exact_CI_lower = if (!is.null(bt)) unname(bt$conf.int[[1L]]) else NA_real_,
    Exact_CI_upper = if (!is.null(bt)) unname(bt$conf.int[[2L]]) else NA_real_,
    Measure = "Ever-detected incident OCD during observed follow-up",
    stringsAsFactors = FALSE
  )
})
crude_incident <- do.call(rbind, crude_incident_rows)
crude_apa <- vapply(
  seq_len(nrow(crude_incident)),
  function(i) {
    x <- crude_incident[i, ]
    
    sprintf(
      "%s: %d of %d participants (%.2f%%) were ever observed with incident OCD during follow-up, exact 95%% CI [%.2f%%, %.2f%%].",
      x$Group,
      x$Incident_OCD_cases,
      x$N,
      x$Crude_percent,
      100 * x$Exact_CI_lower,
      100 * x$Exact_CI_upper
    )
  },
  character(1)
)

origins <- c("Study entry", "Most recent trauma", "Worst trauma")
groups <- list(Overall = NULL, Female = "Female", Male = "Male")

risk_rows <- list()
curve_rows <- list()
risk_set_rows <- list()
rolling_plot_rows <- list()
threshold_rows <- list()
audit_rows <- list()
apa_lines <- character()
row_id <- 0L

for (origin_index in seq_along(origins)) {
  origin <- origins[[origin_index]]
  progress_log("RQ1", paste0("processing time origin: ", origin))
  origin_data <- rq1_make_origin_data(data, origin)

  audit_rows[[length(audit_rows) + 1L]] <- data.frame(
    Time_origin = origin,
    Source_rows = nrow(data),
    Analysis_N = nrow(origin_data),
    Excluded_invalid_or_missing = nrow(data) - nrow(origin_data),
    Events = sum(origin_data$event),
    Delayed_entries = sum(origin_data$entry > 0),
    Minimum_entry_years = if (nrow(origin_data)) min(origin_data$entry) else NA_real_,
    Maximum_entry_years = if (nrow(origin_data)) max(origin_data$entry) else NA_real_,
    Maximum_followup_years = if (nrow(origin_data)) max(origin_data$duration) else NA_real_,
    stringsAsFactors = FALSE
  )

  for (group_name in names(groups)) {
    sex_value <- groups[[group_name]]
    d <- if (is.null(sex_value)) origin_data else origin_data[origin_data$sex == sex_value, , drop = FALSE]
    if (nrow(d) == 0L || sum(d$event) == 0L) {
      apa_lines <- c(apa_lines, sprintf("%s, %s: not estimable because no rows or no events were available.", origin, group_name))
      next
    }

    fit <- rq1_fit_km(d)
    risks <- rq1_risk_at_times(fit, d, REPORT_TIMES, origin, group_name)
    boot <- bootstrap_km(
      d,
      times = REPORT_TIMES,
      B = BOOT_B,
      seed = RNG_SEED + 1000L * origin_index + match(group_name, names(groups)),
      label = paste("RQ1 risk", origin, group_name, sep = " | ")
    )
    boot_ci <- t(apply(boot$risk, 2, percentile_ci))
    risks$Bootstrap_CI_lower <- boot_ci[, 1]
    risks$Bootstrap_CI_upper <- boot_ci[, 2]
    risks$Bootstrap_replicates_requested <- BOOT_B
    risks$Bootstrap_replicates_valid <- colSums(is.finite(boot$risk))

    risk_rows[[length(risk_rows) + 1L]] <- risks
    curve_rows[[length(curve_rows) + 1L]] <- rq1_curve_table(fit, origin, group_name)
    risk_set_rows[[length(risk_set_rows) + 1L]] <- rq1_numbers_at_risk(d, origin = origin, group = group_name)
    rolling_plot_rows[[length(rolling_plot_rows) + 1L]] <- rq1_python_rolling_grid(fit, d, origin, group_name, step = 0.1)
    threshold_rows[[length(threshold_rows) + 1L]] <- rq1_time_to_risk(fit, .10, origin, group_name)

    for (j in seq_len(nrow(risks))) {
      x <- risks[j, ]
      apa_lines <- c(apa_lines, sprintf(
        "%s, %s, at %s year(s): cumulative OCD risk = %s%%, bootstrap 95%% CI [%s%%, %s%%], pointwise 95%% CI [%s%%, %s%%], n at risk = %s, cumulative events = %s (analysis N = %s; total events = %s).",
        origin, group_name, fmt_num(x$Time_years, 1),
        fmt_num(100 * x$Risk), fmt_num(100 * x$Bootstrap_CI_lower), fmt_num(100 * x$Bootstrap_CI_upper),
        fmt_num(100 * x$CI_lower), fmt_num(100 * x$CI_upper),
        x$N_at_risk, x$Events_by_time, x$Analysis_N, x$Total_events
      ))
    }
  }
}

risk_results <- rbind_fill(risk_rows)
curve_results <- rbind_fill(curve_rows)
risk_sets <- rbind_fill(risk_set_rows)
rolling_plot_data <- rbind_fill(rolling_plot_rows)
thresholds <- rbind_fill(threshold_rows)
audits <- rbind_fill(audit_rows)

if (nrow(thresholds)) {
  for (i in seq_len(nrow(thresholds))) {
    x <- thresholds[i, ]
    apa_lines <- c(apa_lines, sprintf(
      "%s, %s: 10%% cumulative risk was %s%s.",
      x$Time_origin,
      x$Group,
      if (isTRUE(x$Threshold_reached)) "reached at " else "not reached",
      if (isTRUE(x$Threshold_reached)) paste0(fmt_num(x$Time_years), " years") else " during observed follow-up"
    ))
  }
}

progress_log("RQ1", "saving figure-ready Kaplan-Meier curve data for Python")
figure1_data_path <- file.path(out_dir, "Figure1_delayed_entry_cumulative_risk_data.csv")
write_csv_safe(curve_results, figure1_data_path)
write_csv_safe(rolling_plot_data, file.path(out_dir, "Figure1_rolling_KM_plot_data.csv"))

rq1_model_specifications <- data.frame(
  Analysis = c("Delayed-entry cumulative risk from study entry", "Delayed-entry cumulative risk from most recent trauma", "Delayed-entry cumulative risk from worst trauma"),
  R_package = "survival",
  R_function = "survival::survfit",
  Model_type = "Delayed-entry Kaplan-Meier estimator",
  Formula = "survival::Surv(entry, duration, event) ~ 1",
  Time_scale = "Years from selected time origin",
  Delayed_entry = TRUE,
  CI_method = "survfit log-log pointwise 95% CI; participant bootstrap also reported at 1 and 2 years",
  PH_assumption = "Not applicable",
  Figure_renderer = "Python / Matplotlib supplied rolling-KM style; R saves both event-time curve data and a dense 0.1-year plot grid",
  stringsAsFactors = FALSE
)
write_csv_safe(rq1_model_specifications, file.path(out_dir, "RQ1_model_specifications.csv"))

save_apa_bundle(
  output_dir = out_dir,
  prefix = "RQ1",
  title = "RQ1: OCD incidence",
  inventory = INVENTORY,
  extra_tables = list(
    analysis_audit = audits,
    cumulative_risk = risk_results,
    curve_data_for_python = curve_results,
    rolling_KM_plot_data = rolling_plot_data,
    numbers_at_risk = risk_sets,
    time_to_10_percent_risk = thresholds,
    crude_ever_detected_incident_OCD = crude_incident,
    baseline_descriptive_summary = baseline_descriptives,
    model_specifications = rq1_model_specifications
  ),
  additional_lines = c(
    "RQ1 uses Kaplan-Meier estimation with delayed entry; the proportional-hazards assumption does not apply.",
    paste0("Participant-level bootstrap replicates requested per estimate: ", BOOT_B, "."),
    "CRUDE EVER-DETECTED INCIDENT OCD PROPORTION",
    crude_apa,
    "These crude proportions do not account for unequal follow-up or censoring; Kaplan-Meier cumulative-risk estimates are the primary incidence estimates.",
    "",
    apa_lines
  )
)

progress_log("RQ1 complete", paste0("outputs: ", normalizePath(out_dir, mustWork = FALSE), "; elapsed=", format_elapsed(SCRIPT_STARTED)))
