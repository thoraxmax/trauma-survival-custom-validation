# Demonstration of the three most custom survival analyses.
# Requires custom_functions_minimal.R in the same directory.

suppressPackageStartupMessages(library(survival))
source("custom_functions_minimal.R")

stopifnot(exists("km_fit", mode = "function"))
stopifnot(exists("rmst_from_fit", mode = "function"))
stopifnot(exists("bootstrap_km", mode = "function"))
stopifnot(exists("percentile_ci", mode = "function"))
stopifnot(exists("build_midpoint_dataset", mode = "function"))
stopifnot(exists("build_cutoff_long", mode = "function"))

set.seed(20260812)
n <- 2700L
n_events <- 250L
event_id <- sample.int(n, n_events, replace = FALSE)

# Simulated study cohort: ages in years; trauma/onset timing in months.
age_at_trauma <- runif(n, 10, 16)
entry_gap <- runif(n, 0.02, 0.20)
age_at_most_recent_trauma <- age_at_trauma
age_at_worst_trauma <- pmax(8, age_at_trauma - runif(n, 0, 0.25))

event_gap <- numeric(n)
censor_gap <- numeric(n)

early_events <- sample(event_id, 50L, replace = FALSE)
later_events <- setdiff(event_id, early_events)
non_events <- setdiff(seq_len(n), event_id)

# Keep the event-time distribution broad so the cutoff effect is not manufactured.
event_gap[early_events] <- runif(length(early_events), 0.25, 0.50)
event_gap[later_events] <- runif(length(later_events), 0.60, 4.50)
censor_gap[event_id] <- event_gap[event_id] + 0.25

event_gap[non_events] <- runif(length(non_events), 0.50, 4.50)
censor_gap[non_events] <- runif(length(non_events), 0.75, 4.00)

follow_gap <- pmin(event_gap, censor_gap)

dat <- data.frame(
  record_id = seq_len(n),
  sex = factor(sample(c("Male", "Female"), n, TRUE), levels = c("Male", "Female")),
  event = as.integer(seq_len(n) %in% event_id),
  age_entry = age_at_trauma + entry_gap,
  age_exit = age_at_trauma + entry_gap + follow_gap,
  age_at_most_recent_trauma = age_at_most_recent_trauma,
  age_at_worst_trauma = age_at_worst_trauma
)

dat$event_gap <- dat$age_exit - dat$age_entry
stopifnot(nrow(dat) == 2700L)
stopifnot(sum(dat$event) == 250L)
stopifnot(all(is.finite(dat$age_entry)))
stopifnot(all(is.finite(dat$age_exit)))
stopifnot(all(dat$age_exit > dat$age_entry))
stopifnot(all(dat$event %in% c(0, 1)))

# 1. Delayed-entry KM + RMST + participant bootstrap
km_dat <- make_origin_data(
  transform(
    dat,
    time_since_start = event_gap * 12,
    months_since_most_recent_at_baseline = 0,
    time_since_most_recent_trauma = event_gap * 12,
    months_since_worst_at_baseline = 0,
    time_since_worst_trauma = event_gap * 12
  ),
  "Most recent trauma"
)

stopifnot(nrow(km_dat) == nrow(dat))
stopifnot(all(km_dat$duration > km_dat$entry))

fit <- km_fit(km_dat)
km <- km_risk_at_times(fit, km_dat, times = c(1, 2))
rmst <- rmst_from_fit(fit, tau = 4)
boot <- bootstrap_km(
  km_dat, times = c(1, 2), tau = 4, B = 200, seed = 1,
  label = "Reviewer demo KM/RMST"
)
ci <- t(apply(boot$risk, 2, percentile_ci))

stopifnot(nrow(km) == 2L)
stopifnot(all(is.finite(km$Risk)))
stopifnot(all(km$Risk >= 0 & km$Risk <= 1))
stopifnot(is.finite(rmst) && rmst > 0 && rmst <= 4)
stopifnot(all(dim(boot$risk) == c(200L, 2L)))
stopifnot(all(is.finite(ci)))
stopifnot(all(ci[, 1] <= ci[, 2]))

cat("\nSimulated cohort: ", nrow(dat), " participants; ", sum(dat$event), " events.\n", sep = "")

cat("\n1. Delayed-entry KM/RMST\n")
print(km[, c("Time_years", "Risk", "CI_lower", "CI_upper", "N_at_risk")], row.names = FALSE)
cat("RMST(4 years):", round(rmst, 3), "\n")
cat("Bootstrap percentile CI (type 6):\n")
print(round(ci, 3))
cat("Validation: delayed-entry risk sets, RMST bounds, and bootstrap output passed.\n")

# 2. Midpoint-onset reconstruction
mid_dat <- dat
mid_dat$first_positive_months <- NA_real_
mid_dat$last_negative_before_first_positive_months <- NA_real_
event_rows <- which(dat$event == 1)
event_follow_months <- dat$event_gap[event_rows] * 12
mid_dat$first_positive_months[event_rows] <- 0.75 * event_follow_months
mid_dat$last_negative_before_first_positive_months[event_rows] <- 0.50 * event_follow_months
mid_dat$time_since_start <- (mid_dat$age_exit - mid_dat$age_entry) * 12
mid_dat$time_since_most_recent_trauma <- mid_dat$time_since_start
mid_dat$time_since_worst_trauma <- mid_dat$time_since_start

mid <- build_midpoint_dataset(mid_dat)
audit <- attr(mid, "midpoint_audit")
valid <- attr(mid, "valid_midpoint")

# Recover the audit from returned columns if attributes are unavailable.
if (is.null(valid)) {
  valid <- mid_dat$event == 1 &
    is.finite(mid_dat$first_positive_months) &
    is.finite(mid_dat$last_negative_before_first_positive_months) &
    mid_dat$first_positive_months > mid_dat$last_negative_before_first_positive_months &
    is.finite(mid$midpoint_months) &
    is.finite(mid$age_exit_midpoint) &
    mid$age_exit_midpoint > mid$age_entry &
    mid$age_exit_midpoint <= mid$age_exit
}

if (!is.data.frame(audit)) {
  audit <- data.frame(
    Events = sum(mid_dat$event == 1, na.rm = TRUE),
    Valid_midpoint_events = sum(valid, na.rm = TRUE),
    Event_cases_without_valid_midpoint = sum(mid_dat$event == 1, na.rm = TRUE) - sum(valid, na.rm = TRUE)
  )
}

stopifnot(length(valid) == n)
stopifnot(audit$Events[[1]] == 250L)
stopifnot(audit$Valid_midpoint_events[[1]] == 250L)
if (audit$Valid_midpoint_events[[1]] > 0L) {
  stopifnot(all(is.finite(mid$midpoint_months[valid])))
  expected_mid <- (
    mid$first_positive_months[valid] +
      mid$last_negative_before_first_positive_months[valid]
  ) / 2
  stopifnot(max(abs(mid$midpoint_months[valid] - expected_mid)) < 1e-12)
  stopifnot(all(mid$age_exit_midpoint[valid] > mid$age_entry[valid]))
  stopifnot(all(mid$age_exit_midpoint[valid] <= mid$age_exit[valid]))
}

cat("\n2. Midpoint-onset reconstruction\n")
cat("Valid midpoint events:", audit$Valid_midpoint_events[[1]], "/", audit$Events[[1]], "\n")
valid_mid <- is.finite(mid$midpoint_months) &
  is.finite(mid$age_exit_midpoint) &
  mid$event == 1
show <- which(valid_mid)[seq_len(min(5L, sum(valid_mid)))]
print(mid[show, c(
  "first_positive_months",
  "last_negative_before_first_positive_months",
  "midpoint_months",
  "age_exit_midpoint"
)], row.names = FALSE)
cat("Validation: midpoint arithmetic and reconstructed event ages passed.\n")

# 3. Post-trauma 6-month cutoff Cox
cutoff_data <- build_cutoff_long(
  dat,
  "age_at_most_recent_trauma",
  cutoff_months = 6L
)

stopifnot(nrow(cutoff_data) > 0L)
stopifnot(all(cutoff_data$stop > cutoff_data$start))
stopifnot(all(cutoff_data$within_cutoff %in% c(0, 1)))
stopifnot(all(tapply(cutoff_data$event, cutoff_data$record_id, sum) <= 1))
stopifnot(sum(cutoff_data$event) == sum(dat$event))

rows_per_id <- table(cutoff_data$record_id)
stopifnot(all(rows_per_id <= 2L))

split_ids <- names(rows_per_id[rows_per_id == 2L])
if (length(split_ids)) {
  for (id in split_ids) {
    x <- cutoff_data[as.character(cutoff_data$record_id) == id, , drop = FALSE]
    x <- x[order(x$start), , drop = FALSE]
    stopifnot(abs(x$stop[1] - x$start[2]) < 1e-12)
    stopifnot(identical(as.integer(x$within_cutoff), c(1L, 0L)))
  }
}

fit_cutoff <- survival::coxph(
  Surv(start, stop, event) ~
    within_cutoff + age_at_trauma + strata(sex) + cluster(record_id),
  data = cutoff_data,
  ties = "efron"
)

cf <- summary(fit_cutoff)$coefficients["within_cutoff", , drop = FALSE]
stopifnot(nrow(cf) == 1L)
stopifnot(all(is.finite(cf[, c("coef", "exp(coef)", "robust se", "z", "Pr(>|z|)")])))
stopifnot(is.finite(exp(cf[, "coef"])))

cat("\n3. Post-trauma 6-month cutoff Cox\n")
print(cf)
cat("Validation: interval splitting, event preservation, continuity, and Cox statistics passed.\n")

cat("\nCompleted: all three custom analyses and validation checks passed.\n")
