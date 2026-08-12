# Minimal helper set for reviewer validation of the three custom analyses.
# No plotting, reporting, screening, or unrelated RQ functions.

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The 'survival' package is required.", call. = FALSE)
}

# -----------------------------
# RQ1/RQ2: delayed-entry KM/RMST
# -----------------------------
km_fit <- function(data) {
  survival::survfit(
    survival::Surv(entry, duration, event) ~ 1,
    data = data,
    conf.type = "log-log"
  )
}

rmst_from_fit <- function(fit, tau) {
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  times <- fit$time[fit$time < tau]
  survs <- fit$surv[fit$time < tau]
  cuts <- c(0, times, tau)
  interval_surv <- c(1, survs)
  sum(diff(cuts) * interval_surv)
}

fast_delayed_km <- function(entry, duration, event, times = numeric(), tau = NULL) {
  keep <- is.finite(entry) & is.finite(duration) & event %in% c(0, 1) &
    entry >= 0 & duration > entry
  entry <- entry[keep]
  duration <- duration[keep]
  event <- event[keep]

  if (!length(duration)) {
    return(list(risk = rep(NA_real_, length(times)), rmst = NA_real_))
  }

  event_times <- sort(unique(duration[event == 1]))
  max_target <- suppressWarnings(max(c(times, tau), na.rm = TRUE))
  if (is.finite(max_target)) event_times <- event_times[event_times <= max_target]

  if (length(event_times)) {
    n_risk <- vapply(event_times, function(t) sum(entry <= t & duration >= t), integer(1))
    n_event <- vapply(event_times, function(t) sum(event == 1 & duration == t), integer(1))
    survival_after <- cumprod(pmax(0, 1 - n_event / n_risk))
  } else {
    survival_after <- numeric()
  }

  risk <- vapply(times, function(t) {
    idx <- which(event_times <= t)
    if (!length(idx)) 0 else 1 - survival_after[max(idx)]
  }, numeric(1))

  rmst <- NA_real_
  if (!is.null(tau) && is.finite(tau) && tau > 0) {
    before_tau <- which(event_times < tau)
    cuts <- c(0, event_times[before_tau], tau)
    interval_survival <- c(1, survival_after[before_tau])
    rmst <- sum(diff(cuts) * interval_survival)
  }

  list(risk = risk, rmst = rmst)
}

bootstrap_km <- function(data, times = c(1, 2), tau = 4, B = 200L,
                         seed = 1L, label = NULL) {
  set.seed(seed)
  n <- nrow(data)
  risk_mat <- matrix(NA_real_, nrow = B, ncol = length(times))
  rmst_vec <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    est <- fast_delayed_km(
      data$entry[idx], data$duration[idx], data$event[idx],
      times = times, tau = tau
    )
    risk_mat[b, ] <- est$risk
    rmst_vec[b] <- est$rmst
  }

  list(risk = risk_mat, rmst = rmst_vec)
}

percentile_ci <- function(x, probs = c(.025, .975)) {
  x <- x[is.finite(x)]
  if (length(x) < 20L) return(c(NA_real_, NA_real_))
  unname(stats::quantile(
    x, probs = probs, na.rm = TRUE, names = FALSE, type = 6
  ))
}

# -----------------------------
# RQ2: midpoint-onset sensitivity
# -----------------------------
build_midpoint_dataset <- function(data) {
  required <- c(
    "event", "age_entry", "age_exit",
    "first_positive_months",
    "last_negative_before_first_positive_months"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Midpoint sensitivity missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  out <- data
  event_case <- out$event == 1

  valid <- event_case &
    is.finite(out$age_entry) &
    is.finite(out$age_exit) &
    out$age_exit > out$age_entry &
    is.finite(out$first_positive_months) &
    is.finite(out$last_negative_before_first_positive_months) &
    out$first_positive_months >= 0 &
    out$last_negative_before_first_positive_months >= 0 &
    out$first_positive_months > out$last_negative_before_first_positive_months

  out$midpoint_months <- NA_real_
  out$midpoint_months[valid] <- (
    out$first_positive_months[valid] +
      out$last_negative_before_first_positive_months[valid]
  ) / 2

  shift <- out$first_positive_months - out$midpoint_months

  out$age_exit_midpoint <- out$age_exit
  out$age_exit_midpoint[valid] <-
    out$age_entry[valid] + out$midpoint_months[valid] / 12

  invalid_age <- valid & (
    !is.finite(out$age_exit_midpoint) |
      out$age_exit_midpoint <= out$age_entry |
      out$age_exit_midpoint > out$age_exit
  )
  valid[invalid_age] <- FALSE

  out$midpoint_months[event_case & !valid] <- NA_real_
  out$age_exit_midpoint[event_case & !valid] <- NA_real_

  duration_cols <- intersect(
    c("time_since_start", "time_since_most_recent_trauma", "time_since_worst_trauma"),
    names(out)
  )

  for (col in duration_cols) {
    new_col <- paste0(col, "_midpoint")
    out[[new_col]] <- out[[col]]
    shifted <- out[[col]] - shift
    invalid_shifted <- valid & (!is.finite(shifted) | shifted < 0)
    out[[new_col]][valid] <- shifted[valid]
    out[[new_col]][invalid_shifted] <- NA_real_
    out[[new_col]][event_case & !valid] <- NA_real_
  }

  attr(out, "valid_midpoint") <- valid
  out
}

# -----------------------------
# RQ2: post-trauma start-stop cutoff
# -----------------------------
build_cutoff_long <- function(data, trauma_age_col,
                              exit_col = "age_exit", cutoff_months = 6L) {
  required <- c("record_id", "sex", "event", "age_entry", exit_col, trauma_age_col)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Cutoff model missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  d <- data[stats::complete.cases(data[, required, drop = FALSE]), required, drop = FALSE]
  d <- d[d[[exit_col]] > d$age_entry & d$event %in% c(0, 1), , drop = FALSE]

  rows <- vector("list", nrow(d) * 2L)
  k <- 0L

  for (i in seq_len(nrow(d))) {
    start <- d$age_entry[i]
    stop <- d[[exit_col]][i]
    cut_age <- d[[trauma_age_col]][i] + cutoff_months / 12

    if (start < cut_age) {
      k <- k + 1L
      first_stop <- min(stop, cut_age)
      rows[[k]] <- data.frame(
        record_id = d$record_id[i],
        sex = as.character(d$sex[i]),
        start = start,
        stop = first_stop,
        event = as.integer(d$event[i] == 1 && stop <= cut_age),
        within_cutoff = 1L,
        age_at_trauma = d[[trauma_age_col]][i],
        stringsAsFactors = FALSE
      )
    }

    if (stop > cut_age) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        record_id = d$record_id[i],
        sex = as.character(d$sex[i]),
        start = max(start, cut_age),
        stop = stop,
        event = as.integer(d$event[i] == 1),
        within_cutoff = 0L,
        age_at_trauma = d[[trauma_age_col]][i],
        stringsAsFactors = FALSE
      )
    }
  }

  if (!k) return(data.frame())
  out <- do.call(rbind, rows[seq_len(k)])
  out <- out[out$stop > out$start, , drop = FALSE]
  out$sex <- factor(as.character(out$sex), levels = c("Male", "Female"))
  out
}
