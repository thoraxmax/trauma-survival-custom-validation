# 00_survival_helpers.R
# Shared utilities for the R survival analyses.
# Dependencies: survival (recommended R >= 4.2).

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))


SURVIVAL_HELPER_API <- 21L

progress_log <- function(stage, detail = NULL) {
  suffix <- if (!is.null(detail) && length(detail) && !is.na(detail) && nzchar(as.character(detail))) {
    paste0(" - ", detail)
  } else {
    ""
  }
  message(sprintf("[%s] %s%s", format(Sys.time(), "%H:%M:%S"), stage, suffix))
}

format_elapsed <- function(started) {
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (!is.finite(seconds)) return("NA")
  if (seconds < 60) return(sprintf("%.1f s", seconds))
  sprintf("%d min %.1f s", floor(seconds / 60), seconds %% 60)
}

p_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
}

assert_helper_api <- function(required_api, required_functions = character()) {
  current_api <- get0("SURVIVAL_HELPER_API", ifnotfound = 0L, inherits = TRUE)
  missing_functions <- required_functions[!vapply(
    required_functions,
    function(fn) exists(fn, mode = "function", inherits = TRUE),
    logical(1)
  )]
  if (current_api < required_api || length(missing_functions)) {
    stop(
      "The sourced 00_survival_helpers.R is outdated or incomplete. ",
      "Required helper API: ", required_api, "; loaded API: ", current_api, ".",
      if (length(missing_functions)) paste0(" Missing functions: ", paste(missing_functions, collapse = ", "), ".") else "",
      " Replace 00_survival_helpers.R with the copy distributed beside this script.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

PACKAGE_DIR <- normalizePath(
  Sys.getenv("SURVIVAL_PACKAGE_DIR", getwd()),
  mustWork = FALSE
)
INPUT_CSV_ENV <- trimws(Sys.getenv("SURVIVAL_DATA_PATH", ""))
if (nzchar(INPUT_CSV_ENV)) {
  INPUT_CSV <- normalizePath(INPUT_CSV_ENV, mustWork = FALSE)
} else {
  INPUT_CANDIDATES <- c(
    file.path(PACKAGE_DIR, "cox_survival_dataset_full.csv"),
    file.path(dirname(PACKAGE_DIR), "cox_survival_dataset_full.csv")
  )
  existing_input <- INPUT_CANDIDATES[file.exists(INPUT_CANDIDATES)]
  INPUT_CSV <- normalizePath(if (length(existing_input)) existing_input[[1L]] else INPUT_CANDIDATES[[1L]], mustWork = FALSE)
}
OUTPUT_ENV <- trimws(Sys.getenv("SURVIVAL_OUTPUT_DIR", ""))
if (nzchar(OUTPUT_ENV)) {
  output_candidate <- if (grepl("^(/|[A-Za-z]:[/\\])", OUTPUT_ENV)) OUTPUT_ENV else file.path(PACKAGE_DIR, OUTPUT_ENV)
  OUTPUT_ROOT <- normalizePath(output_candidate, mustWork = FALSE)
} else {
  OUTPUT_ROOT <- normalizePath(file.path(PACKAGE_DIR, "r_results"), mustWork = FALSE)
}
BOOT_B <- as.integer(Sys.getenv("SURVIVAL_BOOTSTRAPS", "1000"))
RNG_SEED <- as.integer(Sys.getenv("SURVIVAL_SEED", "20260806"))
PH_ALPHA <- as.numeric(Sys.getenv("PH_ALPHA", "0.05"))
REPORT_TIMES <- c(1, 2)
REPORT_AGES <- c(10, 12, 15, 18, 20)

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

fmt_num <- function(x, digits = 2L) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p))))
}

first_existing <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit)) hit[[1L]] else NA_character_
}

rbind_fill <- function(items) {
  items <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L, items)
  if (!length(items)) return(data.frame())
  all_names <- unique(unlist(lapply(items, names), use.names = FALSE))
  items <- lapply(items, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[all_names]
  })
  rownames_out <- NULL
  out <- do.call(rbind, items)
  rownames(out) <- rownames_out
  out
}

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
  message("Saved: ", path, " (", nrow(x), " rows)")
}

require_columns <- function(data, cols, context = "analysis") {
  missing <- setdiff(cols, names(data))
  if (length(missing)) {
    stop(context, " requires missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

coerce_numeric <- function(data, cols) {
  for (col in intersect(cols, names(data))) {
    data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
  }
  data
}

standardize_analysis_names <- function(data) {
  aliases <- c(
    time_since_most_recent_trauma = "time_since_latest",
    time_since_worst_trauma = "time_since_worst",
    months_since_most_recent_at_baseline = "months_since_latest_at_baseline",
    months_since_worst_at_baseline = "months_since_worst_at_baseline",
    age_at_most_recent_trauma = "age_at_latest",
    age_at_worst_trauma = "age_at_worst",
    baseline_most_recent_trauma_date = "trauma_latest_date",
    baseline_worst_trauma_date = "trauma_worst_date"
  )
  for (new_name in names(aliases)) {
    old_name <- aliases[[new_name]]
    if (!(new_name %in% names(data)) && old_name %in% names(data)) {
      data[[new_name]] <- data[[old_name]]
    }
  }
  data
}

read_survival_data <- function(path = INPUT_CSV) {
  if (!file.exists(path)) {
    stop(
      "Analysis dataset not found: ", path,
      "\nSet SURVIVAL_DATA_PATH or place cox_survival_dataset_full.csv in the working directory.",
      call. = FALSE
    )
  }
  data <- read.csv(path, check.names = FALSE, na.strings = c("", "NA", "NaN"))
  data <- standardize_analysis_names(data)

  if (!("record_id" %in% names(data))) data$record_id <- seq_len(nrow(data))
  if ("sex" %in% names(data)) {
    data$sex <- factor(trimws(as.character(data$sex)), levels = c("Male", "Female"))
    data$sex_female <- ifelse(is.na(data$sex), NA_real_, as.numeric(data$sex == "Female"))
  }

  numeric_candidates <- c(
    "event", "age_entry", "age_exit", "time_since_start",
    "time_since_most_recent_trauma", "time_since_worst_trauma",
    "months_since_most_recent_at_baseline", "months_since_worst_at_baseline",
    "age_at_most_recent_trauma", "age_at_worst_trauma",
    "first_positive_months", "last_negative_before_first_positive_months",
    "assessment_interval_months", "trauma_type_count",
    "new_worst_trauma_flag", "new_worst_trauma_months_since_start",
    "followup_interval_worst_trauma_flag",
    "followup_interval_worst_trauma_months_since_start"
  )
  numeric_candidates <- unique(c(numeric_candidates, grep("^(sexual_abuse|physical_abuse|emotional_abuse|bullying|accident_experienced_witnessed|natural_disaster|illness_injury|kidnapping|family_violence|community_violence)$", names(data), value = TRUE)))
  data <- coerce_numeric(data, numeric_candidates)
  data
}

analysis_inventory <- function(title, analyses) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
  for (i in seq_along(analyses)) cat(sprintf("%2d. %s\n", i, analyses[[i]]))
  cat("\n")
}

complete_model_data <- function(data, vars) {
  vars <- unique(vars)
  vars <- vars[vars %in% names(data)]
  data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
}

make_formula <- function(response, covariates = character(), strata = character(), cluster = NULL) {
  rhs <- c(covariates, sprintf("strata(%s)", strata))
  if (!is.null(cluster) && nzchar(cluster)) rhs <- c(rhs, sprintf("cluster(%s)", cluster))
  if (!length(rhs)) rhs <- "1"
  as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}

# Core model fitting moved to the individual RQ scripts for transparency.

extract_cox_terms <- function(fit, model_name = attr(fit, "model_name") %||% "Cox model") {
  s <- summary(fit)
  co <- as.data.frame(s$coefficients, check.names = FALSE)
  if (!nrow(co)) return(data.frame())
  beta <- co[["coef"]]
  se_name <- if ("robust se" %in% names(co)) "robust se" else "se(coef)"
  se <- co[[se_name]]
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  analysis_data <- attr(fit, "analysis_data")
  participants <- if (is.data.frame(analysis_data) && "record_id" %in% names(analysis_data)) {
    length(unique(analysis_data$record_id))
  } else {
    fit$n
  }
  variance_estimator <- if (identical(se_name, "robust se")) {
    "Cluster-robust sandwich"
  } else {
    "Model-based"
  }
  out <- data.frame(
    Model = model_name,
    Term = rownames(co),
    B = beta,
    SE = se,
    HR = exp(beta),
    CI_lower = exp(beta - qnorm(.975) * se),
    CI_upper = exp(beta + qnorm(.975) * se),
    z = z,
    p = p,
    N = participants,
    Analysis_rows = fit$n,
    Events = fit$nevent,
    Variance_estimator = variance_estimator,
    Robust_SE = identical(se_name, "robust se"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out
}

cox_model_loglik <- function(fit) {
  value <- tryCatch(
    suppressWarnings(as.numeric(stats::logLik(fit))[[1L]]),
    error = function(e) NA_real_
  )
  if (is.finite(value)) return(unname(value))

  values <- suppressWarnings(as.numeric(fit$loglik))
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  unname(values[[length(values)]])
}

cox_loglik_values <- function(fit) {
  values <- suppressWarnings(as.numeric(fit$loglik))
  values <- values[is.finite(values)]
  model_value <- cox_model_loglik(fit)
  null_value <- if (length(values)) values[[1L]] else model_value
  c(null = unname(null_value), model = unname(model_value))
}

extract_model_stats <- function(fit, model_name = attr(fit, "model_name") %||% "Cox model") {
  s <- summary(fit)
  get_test <- function(x, field) {
    if (is.null(x) || !(field %in% names(x))) return(NA_real_)
    unname(x[[field]])
  }
  concordance <- s$concordance %||% c(concordance = NA_real_, se = NA_real_)
  data.frame(
    Model = model_name,
    N = fit$n,
    Events = fit$nevent,
    Coefficients = length(stats::coef(fit)),
    Log_likelihood_null = unname(cox_loglik_values(fit)[["null"]]),
    Log_likelihood_model = unname(cox_loglik_values(fit)[["model"]]),
    AIC = tryCatch(AIC(fit), error = function(e) NA_real_),
    Concordance = unname(concordance[[1L]]),
    Concordance_SE = if (length(concordance) >= 2L) unname(concordance[[2L]]) else NA_real_,
    LR_chisq = get_test(s$logtest, "test"),
    LR_df = get_test(s$logtest, "df"),
    LR_p = get_test(s$logtest, "pvalue"),
    Wald_chisq = get_test(s$waldtest, "test"),
    Wald_df = get_test(s$waldtest, "df"),
    Wald_p = get_test(s$waldtest, "pvalue"),
    Score_chisq = get_test(s$sctest, "test"),
    Score_df = get_test(s$sctest, "df"),
    Score_p = get_test(s$sctest, "pvalue"),
    Formula = paste(deparse(formula(fit)), collapse = " "),
    stringsAsFactors = FALSE
  )
}

cox_ph_table <- function(
  fit,
  model_name = attr(fit, "model_name") %||% "Cox model",
  main_exposures = character(),
  nuisance_terms = character(),
  transform = "rank"
) {
  out <- tryCatch({
    zph <- cox.zph(fit, transform = transform, terms = TRUE, singledf = FALSE)
    tab <- as.data.frame(zph$table, check.names = FALSE)
    data.frame(
      Model = model_name,
      PH_term = rownames(tab),
      PH_chisq = tab[["chisq"]],
      PH_df = tab[["df"]],
      PH_p = tab[["p"]],
      PH_concern = tab[["p"]] < PH_ALPHA,
      PH_transform = transform,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) {
    data.frame(
      Model = model_name,
      PH_term = "ERROR",
      PH_chisq = NA_real_,
      PH_df = NA_real_,
      PH_p = NA_real_,
      PH_concern = NA,
      PH_transform = transform,
      PH_error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  out$Planned_action <- "No action unless supported by residual pattern and scientific role"
  is_global <- out$PH_term == "GLOBAL"
  is_main <- out$PH_term %in% main_exposures
  is_nuisance <- out$PH_term %in% nuisance_terms
  out$Planned_action[is_global] <- "Review variable-level diagnostics and residual plots"
  out$Planned_action[is_main & out$PH_concern %in% TRUE] <- "Retain exposure and fit age-varying coefficient model"
  out$Planned_action[is_nuisance & out$PH_concern %in% TRUE] <- "Consider clinically coherent stratification or cautious interpretation"
  out
}

safe_path_component <- function(x, fallback = "model", max_chars = 100L) {
  value <- trimws(as.character(x %||% ""))
  converted <- iconv(value, to = "ASCII//TRANSLIT", sub = "")
  if (!is.na(converted) && nzchar(converted)) value <- converted
  value <- gsub("[^A-Za-z0-9._-]+", "_", value)
  value <- gsub("_+", "_", value)
  value <- gsub("^[_-]+|[_-]+$", "", value)
  if (!nzchar(value)) fallback else substr(value, 1L, max(1L, as.integer(max_chars)))
}

ph_plot_device <- function(path, type = c("png", "svg"), width_cm = 16, height_cm = 12) {
  type <- match.arg(type)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  width_in <- width_cm / 2.54
  height_in <- height_cm / 2.54
  if (type == "png") {
    grDevices::png(path, width = width_in, height = height_in, units = "in", res = 600, pointsize = 9, bg = "white")
  } else {
    grDevices::svg(path, width = width_in, height = height_in, pointsize = 9, bg = "white")
  }
  invisible(path)
}

write_ph_plot_pair <- function(path_base, width_cm, height_cm, draw) {
  rows <- list()
  for (type in c("png", "svg")) {
    path <- paste0(path_base, ".", type)
    opened <- FALSE
    error_message <- ""
    tryCatch({
      ph_plot_device(path, type = type, width_cm = width_cm, height_cm = height_cm)
      opened <- TRUE
      draw()
      grDevices::dev.off()
      opened <- FALSE
    }, error = function(e) {
      error_message <<- conditionMessage(e)
    }, finally = {
      if (opened && grDevices::dev.cur() > 1L) {
        try(grDevices::dev.off(), silent = TRUE)
      }
    })
    exists_nonempty <- file.exists(path) && is.finite(file.info(path)$size) && file.info(path)$size > 0
    rows[[length(rows) + 1L]] <- data.frame(
      File_type = type,
      File_path = normalizePath(path, mustWork = FALSE),
      Saved = exists_nonempty && !nzchar(error_message),
      Plot_error = error_message,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

save_cox_zph_plots <- function(
  fit,
  output_dir,
  model_name = attr(fit, "model_name") %||% "Cox model",
  model_phase = "Final model",
  transform = "rank",
  ph_table = NULL,
  save_individual_terms = TRUE
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  model_stem <- safe_path_component(model_name, "cox_model", max_chars = 100L)
  phase_stem <- safe_path_component(model_phase, "phase", max_chars = 30L)
  stem <- paste(model_stem, phase_stem, sep = "__")

  zph <- tryCatch(
    cox.zph(fit, transform = transform, terms = TRUE, singledf = FALSE),
    error = function(e) e
  )
  if (inherits(zph, "error")) {
    return(data.frame(
      Model = model_name,
      Model_phase = model_phase,
      PH_term = "ERROR",
      PH_p = NA_real_,
      PH_concern = NA,
      All_terms_PNG = "",
      All_terms_SVG = "",
      Individual_PNG = "",
      Individual_SVG = "",
      All_terms_saved = FALSE,
      Individual_saved = FALSE,
      Plot_error = conditionMessage(zph),
      stringsAsFactors = FALSE
    ))
  }

  term_names <- colnames(zph$y)
  if (is.null(term_names) || !length(term_names)) {
    term_names <- rownames(zph$table)
    term_names <- term_names[term_names != "GLOBAL"]
  }
  n_terms <- length(term_names)
  if (!n_terms) {
    return(data.frame(
      Model = model_name, Model_phase = model_phase, PH_term = "NONE",
      PH_p = NA_real_, PH_concern = NA,
      All_terms_PNG = "", All_terms_SVG = "",
      Individual_PNG = "", Individual_SVG = "",
      All_terms_saved = FALSE, Individual_saved = FALSE,
      Plot_error = "cox.zph returned no plottable coefficient terms",
      stringsAsFactors = FALSE
    ))
  }

  tab <- as.data.frame(zph$table, check.names = FALSE)
  get_term_p <- function(term) {
    idx <- which(rownames(tab) == term)
    if (!length(idx)) idx <- which(startsWith(rownames(tab), term))
    if (!length(idx) || !("p" %in% names(tab))) return(NA_real_)
    value <- suppressWarnings(as.numeric(tab[["p"]][idx[[1L]]]))
    if (length(value) && is.finite(value)) value else NA_real_
  }

  ncol_plot <- min(3L, max(1L, ceiling(sqrt(n_terms))))
  nrow_plot <- ceiling(n_terms / ncol_plot)
  all_base <- file.path(output_dir, paste0(stem, "__all_terms"))
  all_files <- write_ph_plot_pair(
    all_base,
    width_cm = max(16, 6.2 * ncol_plot),
    height_cm = max(11, 5.4 * nrow_plot),
    draw = function() {
      old <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old), add = TRUE)
      graphics::par(mfrow = c(nrow_plot, ncol_plot), mar = c(4.2, 4.4, 3.1, 1.0), oma = c(0, 0, 2.0, 0))
      for (i in seq_len(n_terms)) {
        p_value <- get_term_p(term_names[[i]])
        plot(
          zph, var = i, resid = TRUE, se = TRUE,
          xlab = "Attained age (years)",
          ylab = expression(beta(t)),
          main = paste0(term_names[[i]], "; PH p ", fmt_p(p_value))
        )
        graphics::abline(h = 0, lty = 3, lwd = 1)
      }
      if (n_terms < nrow_plot * ncol_plot) {
        for (i in seq_len(nrow_plot * ncol_plot - n_terms)) graphics::plot.new()
      }
      graphics::mtext(paste0(model_name, " — scaled Schoenfeld residuals"), outer = TRUE, cex = 1.05, font = 2)
    }
  )
  all_png <- all_files$File_path[all_files$File_type == "png"] %||% ""
  all_svg <- all_files$File_path[all_files$File_type == "svg"] %||% ""
  all_saved <- all(all_files$Saved)
  all_error <- paste(unique(all_files$Plot_error[nzchar(all_files$Plot_error)]), collapse = "; ")

  audit <- vector("list", n_terms)
  for (i in seq_len(n_terms)) {
    term <- term_names[[i]]
    p_value <- get_term_p(term)
    individual_png <- ""
    individual_svg <- ""
    individual_saved <- NA
    individual_error <- ""
    if (isTRUE(save_individual_terms)) {
      term_base <- file.path(output_dir, paste0(stem, "__term_", safe_path_component(term, paste0("term_", i), max_chars = 80L)))
      individual_files <- write_ph_plot_pair(
        term_base,
        width_cm = 13,
        height_cm = 9,
        draw = function() {
          old <- graphics::par(no.readonly = TRUE)
          on.exit(graphics::par(old), add = TRUE)
          graphics::par(mar = c(4.6, 4.8, 3.6, 1.2))
          plot(
            zph, var = i, resid = TRUE, se = TRUE,
            xlab = "Attained age (years)",
            ylab = expression(beta(t)),
            main = paste0(model_name, "\n", term, "; PH p ", fmt_p(p_value))
          )
          graphics::abline(h = 0, lty = 3, lwd = 1)
        }
      )
      individual_png <- individual_files$File_path[individual_files$File_type == "png"] %||% ""
      individual_svg <- individual_files$File_path[individual_files$File_type == "svg"] %||% ""
      individual_saved <- all(individual_files$Saved)
      individual_error <- paste(unique(individual_files$Plot_error[nzchar(individual_files$Plot_error)]), collapse = "; ")
    }
    audit[[i]] <- data.frame(
      Model = model_name,
      Model_phase = model_phase,
      PH_term = term,
      PH_p = p_value,
      PH_concern = if (is.finite(p_value)) p_value < PH_ALPHA else NA,
      All_terms_PNG = all_png,
      All_terms_SVG = all_svg,
      Individual_PNG = individual_png,
      Individual_SVG = individual_svg,
      All_terms_saved = all_saved,
      Individual_saved = individual_saved,
      Plot_error = paste(c(all_error, individual_error)[nzchar(c(all_error, individual_error))], collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, audit)
  progress_log(
    "PH residual plots",
    sprintf("%s [%s]: %d term(s), all-term pair saved=%s", model_name, model_phase, n_terms, all_saved)
  )
  out
}


ph_term_p <- function(ph_table, term) {
  if (is.null(ph_table) || !nrow(ph_table) || !("PH_term" %in% names(ph_table))) return(NA_real_)
  idx <- which(ph_table$PH_term == term)
  if (!length(idx)) idx <- grep(paste0("^", term), ph_table$PH_term)
  values <- suppressWarnings(as.numeric(ph_table$PH_p[idx]))
  values <- values[is.finite(values)]
  if (length(values)) values[[1L]] else NA_real_
}

scalar_numeric_or_na <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  value <- value[is.finite(value)]
  if (length(value)) value[[1L]] else NA_real_
}

sex_effect_from_fit <- function(fit) {
  terms <- extract_cox_terms(fit, attr(fit, "model_name") %||% "Sex PH screening model")
  empty <- data.frame(
    HR = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    coefficient_p = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!nrow(terms) || !("Term" %in% names(terms))) return(empty)

  idx <- grep("^sex($|[^[:alnum:]_])|^sex", as.character(terms$Term), ignore.case = TRUE)
  if (!length(idx)) return(empty)

  i <- idx[[1L]]
  data.frame(
    HR = scalar_numeric_or_na(terms$HR[i]),
    CI_lower = scalar_numeric_or_na(terms$CI_lower[i]),
    CI_upper = scalar_numeric_or_na(terms$CI_upper[i]),
    coefficient_p = scalar_numeric_or_na(terms$p[i]),
    stringsAsFactors = FALSE
  )
}

# Core model fitting moved to the individual RQ scripts for transparency.

# Core model fitting moved to the individual RQ scripts for transparency.

extract_pooled_sex_interaction <- function(
  fit,
  main_term,
  interaction_term,
  effect_name,
  unit_label,
  center_value = NA_real_
) {
  b <- stats::coef(fit)
  V <- stats::vcov(fit)
  analysis_data <- attr(fit, "analysis_data")
  participants <- if (is.data.frame(analysis_data) && "record_id" %in% names(analysis_data)) {
    length(unique(analysis_data$record_id))
  } else {
    fit$n
  }
  coefficient_table <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  robust_variance <- !is.null(coefficient_table) && "robust se" %in% colnames(coefficient_table)
  variance_estimator <- if (robust_variance) "Cluster-robust sandwich" else "Model-based"
  required_terms <- c(main_term, interaction_term)
  if (!all(required_terms %in% names(b))) {
    stop(
      "Could not identify the required interaction coefficients: ",
      paste(setdiff(required_terms, names(b)), collapse = ", "),
      call. = FALSE
    )
  }
  if (!all(is.finite(b[required_terms]))) {
    stop("The pooled sex-interaction coefficients are not finite.", call. = FALSE)
  }

  male_log_hr <- unname(b[[main_term]])
  female_log_hr <- unname(b[[main_term]] + b[[interaction_term]])
  male_var <- unname(V[main_term, main_term])
  female_var <- unname(
    V[main_term, main_term] + V[interaction_term, interaction_term] +
      2 * V[main_term, interaction_term]
  )
  male_se <- sqrt(max(0, male_var))
  female_se <- sqrt(max(0, female_var))
  interaction_se <- sqrt(max(0, unname(V[interaction_term, interaction_term])))
  interaction_beta <- unname(b[[interaction_term]])

  make_row <- function(sex_label, log_hr, se) {
    z <- if (is.finite(se) && se > 0) log_hr / se else NA_real_
    data.frame(
      Model = attr(fit, "model_name") %||% "Pooled sex-interaction model",
      Effect = effect_name,
      Sex = sex_label,
      Unit = unit_label,
      Center_value = center_value,
      log_HR = log_hr,
      SE = se,
      HR = exp(log_hr),
      CI_lower = exp(log_hr - stats::qnorm(.975) * se),
      CI_upper = exp(log_hr + stats::qnorm(.975) * se),
      z = z,
      p = if (is.finite(z)) 2 * stats::pnorm(abs(z), lower.tail = FALSE) else NA_real_,
      N = participants,
      Analysis_rows = fit$n,
      Events = fit$nevent,
      Variance_estimator = variance_estimator,
      stringsAsFactors = FALSE
    )
  }

  sex_specific <- rbind(
    make_row("Male", male_log_hr, male_se),
    make_row("Female", female_log_hr, female_se)
  )
  interaction_z <- if (is.finite(interaction_se) && interaction_se > 0) interaction_beta / interaction_se else NA_real_
  interaction_test <- data.frame(
    Model = attr(fit, "model_name") %||% "Pooled sex-interaction model",
    Effect = effect_name,
    Contrast = "Female versus male difference in the effect",
    Unit = unit_label,
    Center_value = center_value,
    Interaction_term = interaction_term,
    B_interaction = interaction_beta,
    SE_interaction = interaction_se,
    HR_ratio_Female_vs_Male = exp(interaction_beta),
    CI_lower = exp(interaction_beta - stats::qnorm(.975) * interaction_se),
    CI_upper = exp(interaction_beta + stats::qnorm(.975) * interaction_se),
    z = interaction_z,
    p = if (is.finite(interaction_z)) 2 * stats::pnorm(abs(interaction_z), lower.tail = FALSE) else NA_real_,
    N = fit$n,
    Events = fit$nevent,
    stringsAsFactors = FALSE
  )
  list(sex_specific = sex_specific, interaction_test = interaction_test)
}

# Core model fitting moved to the individual RQ scripts for transparency.

# Core model fitting moved to the individual RQ scripts for transparency.

wald_joint_test <- function(fit, terms, label = "Joint Wald test") {
  b <- coef(fit)
  idx <- which(names(b) %in% terms)
  if (!length(idx)) {
    idx <- unique(unlist(lapply(terms, function(term) grep(paste0("^", term), names(b)))))
  }
  if (!length(idx)) {
    return(data.frame(Test = label, Chi_square = NA_real_, df = 0L, p = NA_real_))
  }
  V <- vcov(fit)[idx, idx, drop = FALSE]
  b_sub <- b[idx]
  stat <- tryCatch(as.numeric(t(b_sub) %*% solve(V, b_sub)), error = function(e) NA_real_)
  data.frame(
    Test = label,
    Chi_square = stat,
    df = length(idx),
    p = if (is.finite(stat)) pchisq(stat, df = length(idx), lower.tail = FALSE) else NA_real_,
    Terms = paste(names(b_sub), collapse = "; "),
    stringsAsFactors = FALSE
  )
}


likelihood_ratio_block_test <- function(fit, block_terms, label = "Omnibus likelihood-ratio test") {
  fit_data <- attr(fit, "analysis_data")
  if (is.null(fit_data) || !is.data.frame(fit_data)) {
    return(data.frame(
      Test = label, LR_chisq = NA_real_, df = NA_integer_, p = NA_real_,
      Full_log_likelihood = NA_real_, Reduced_log_likelihood = NA_real_,
      Block_terms = paste(block_terms, collapse = "; "),
      Error = "Analysis data were not attached to the fitted model",
      stringsAsFactors = FALSE
    ))
  }

  full_formula <- formula(fit)
  term_labels <- attr(terms(full_formula), "term.labels")
  remove_labels <- term_labels[term_labels %in% block_terms]
  if (!length(remove_labels)) {
    remove_labels <- unique(unlist(lapply(block_terms, function(term) {
      term_labels[grepl(paste0("^", term, "($|:|\\()"), term_labels)]
    })))
  }
  if (!length(remove_labels)) {
    return(data.frame(
      Test = label, LR_chisq = NA_real_, df = 0L, p = NA_real_,
      Full_log_likelihood = cox_model_loglik(fit), Reduced_log_likelihood = NA_real_,
      Block_terms = paste(block_terms, collapse = "; "),
      Error = "No requested block term was present in the full model",
      stringsAsFactors = FALSE
    ))
  }

  keep_labels <- setdiff(term_labels, remove_labels)
  response_text <- paste(deparse(full_formula[[2L]]), collapse = "")
  reduced_formula <- if (length(keep_labels)) {
    as.formula(paste(response_text, "~", paste(keep_labels, collapse = " + ")))
  } else {
    as.formula(paste(response_text, "~ 1"))
  }

  progress_log("Omnibus LRT", paste0("fitting reduced model: ", label))
  reduced_fit <- tryCatch(
    coxph(
      reduced_formula,
      data = fit_data,
      ties = "efron",
      x = TRUE,
      y = TRUE,
      model = TRUE,
      singular.ok = TRUE
    ),
    error = function(e) e
  )
  if (inherits(reduced_fit, "error")) {
    return(data.frame(
      Test = label, LR_chisq = NA_real_, df = NA_integer_, p = NA_real_,
      Full_log_likelihood = cox_model_loglik(fit), Reduced_log_likelihood = NA_real_,
      Block_terms = paste(remove_labels, collapse = "; "),
      Error = conditionMessage(reduced_fit),
      stringsAsFactors = FALSE
    ))
  }

  full_df <- sum(is.finite(coef(fit)))
  reduced_df <- sum(is.finite(coef(reduced_fit)))
  df_difference <- full_df - reduced_df
  full_ll <- cox_model_loglik(fit)
  reduced_ll <- cox_model_loglik(reduced_fit)
  statistic <- 2 * (full_ll - reduced_ll)
  progress_log("Omnibus LRT", paste0("completed: ", label, "; LR chi-square=", if (is.finite(statistic)) sprintf("%.3f", statistic) else "NA", ", df=", df_difference))
  if (is.finite(statistic) && statistic < 0 && abs(statistic) < 1e-8) statistic <- 0

  data.frame(
    Test = label,
    LR_chisq = statistic,
    df = df_difference,
    p = if (is.finite(statistic) && df_difference > 0L) {
      pchisq(statistic, df = df_difference, lower.tail = FALSE)
    } else NA_real_,
    Full_log_likelihood = full_ll,
    Reduced_log_likelihood = reduced_ll,
    Block_terms = paste(remove_labels, collapse = "; "),
    Reduced_formula = paste(deparse(reduced_formula), collapse = " "),
    Error = NA_character_,
    stringsAsFactors = FALSE
  )
}

numeric_vif_table <- function(data, predictors, model_name = "Model") {
  predictors <- unique(predictors[predictors %in% names(data)])
  if (length(predictors) < 2L) return(data.frame())
  d <- data[complete.cases(data[, predictors, drop = FALSE]), predictors, drop = FALSE]
  for (v in predictors) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d <- d[complete.cases(d), , drop = FALSE]
  if (nrow(d) < length(predictors) + 2L) return(data.frame())

  rows <- lapply(predictors, function(variable) {
    others <- setdiff(predictors, variable)
    y <- d[[variable]]
    if (!is.finite(sd(y)) || sd(y) == 0 || !length(others)) {
      return(data.frame(
        Model = model_name, Variable = variable, N = nrow(d), R_squared = NA_real_,
        Tolerance = NA_real_, VIF = NA_real_, VIF_gt_5 = NA, VIF_gt_10 = NA,
        stringsAsFactors = FALSE
      ))
    }
    aux <- tryCatch(lm(reformulate(others, response = variable), data = d), error = function(e) NULL)
    r2 <- if (is.null(aux)) NA_real_ else summary(aux)$r.squared
    tolerance <- if (is.finite(r2)) 1 - r2 else NA_real_
    vif <- if (is.finite(tolerance) && tolerance > 0) 1 / tolerance else Inf
    data.frame(
      Model = model_name,
      Variable = variable,
      N = nrow(d),
      R_squared = r2,
      Tolerance = tolerance,
      VIF = vif,
      VIF_gt_5 = is.finite(vif) && vif > 5,
      VIF_gt_10 = is.finite(vif) && vif > 10,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

predictor_correlation_long <- function(data, predictors, model_name = "Model") {
  predictors <- unique(predictors[predictors %in% names(data)])
  if (length(predictors) < 2L) return(data.frame())
  d <- data[complete.cases(data[, predictors, drop = FALSE]), predictors, drop = FALSE]
  for (v in predictors) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  cmat <- suppressWarnings(cor(d, use = "pairwise.complete.obs"))
  idx <- which(upper.tri(cmat), arr.ind = TRUE)
  data.frame(
    Model = model_name,
    Variable_1 = rownames(cmat)[idx[, 1]],
    Variable_2 = colnames(cmat)[idx[, 2]],
    Correlation = cmat[idx],
    Absolute_correlation = abs(cmat[idx]),
    stringsAsFactors = FALSE
  )
}

multicollinearity_summary <- function(data, predictors, vif_table, model_name = "Model") {
  predictors <- unique(predictors[predictors %in% names(data)])
  if (length(predictors) < 2L) return(data.frame())
  d <- data[complete.cases(data[, predictors, drop = FALSE]), predictors, drop = FALSE]
  for (v in predictors) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  x <- as.matrix(d)
  x <- x[, apply(x, 2, function(z) is.finite(sd(z)) && sd(z) > 0), drop = FALSE]
  condition_number <- NA_real_
  if (nrow(x) > 1L && ncol(x) > 1L) {
    z <- scale(x)
    sv <- svd(z, nu = 0, nv = 0)$d
    condition_number <- if (length(sv) && min(sv) > .Machine$double.eps) max(sv) / min(sv) else Inf
  }
  corr <- predictor_correlation_long(d, predictors, model_name)
  max_corr <- if (nrow(corr)) max(corr$Absolute_correlation, na.rm = TRUE) else NA_real_
  finite_vif <- vif_table$VIF[is.finite(vif_table$VIF)]
  data.frame(
    Model = model_name,
    N = nrow(d),
    Predictors = length(predictors),
    Maximum_VIF = if (length(finite_vif)) max(finite_vif) else NA_real_,
    Mean_VIF = if (length(finite_vif)) mean(finite_vif) else NA_real_,
    Maximum_absolute_pairwise_correlation = max_corr,
    Condition_number = condition_number,
    Any_VIF_gt_5 = any(vif_table$VIF_gt_5 %in% TRUE),
    Any_VIF_gt_10 = any(vif_table$VIF_gt_10 %in% TRUE),
    stringsAsFactors = FALSE
  )
}

build_new_worst_long <- function(d) {
  required <- c(
    "record_id", "sex", "event", "age_entry", "age_exit",
    "new_worst_trauma_flag", "new_worst_trauma_months_since_start"
  )
  require_columns(d, required, "New worst trauma time-varying model")
  keep <- complete.cases(d[, c("record_id", "sex", "event", "age_entry", "age_exit")])
  d <- d[keep & d$age_exit > d$age_entry & d$event %in% c(0, 1), , drop = FALSE]

  flagged_missing_date <- d$new_worst_trauma_flag == 1 & !is.finite(d$new_worst_trauma_months_since_start)
  d <- d[!flagged_missing_date, , drop = FALSE]
  switch_age <- d$age_entry + d$new_worst_trauma_months_since_start / 12
  valid_switch <- d$new_worst_trauma_flag == 1 &
    is.finite(switch_age) & switch_age > d$age_entry & switch_age < d$age_exit
  at_entry <- d$new_worst_trauma_flag == 1 & is.finite(switch_age) & switch_age <= d$age_entry

  rows <- vector("list", nrow(d) * 2L)
  k <- 0L
  for (i in seq_len(nrow(d))) {
    if (valid_switch[[i]]) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        record_id = d$record_id[[i]], sex = as.character(d$sex[[i]]),
        start = d$age_entry[[i]], stop = switch_age[[i]], event = 0L,
        new_worst_trauma_tv = 0, stringsAsFactors = FALSE
      )
      k <- k + 1L
      rows[[k]] <- data.frame(
        record_id = d$record_id[[i]], sex = as.character(d$sex[[i]]),
        start = switch_age[[i]], stop = d$age_exit[[i]], event = d$event[[i]],
        new_worst_trauma_tv = 1, stringsAsFactors = FALSE
      )
    } else {
      k <- k + 1L
      rows[[k]] <- data.frame(
        record_id = d$record_id[[i]], sex = as.character(d$sex[[i]]),
        start = d$age_entry[[i]], stop = d$age_exit[[i]], event = d$event[[i]],
        new_worst_trauma_tv = as.numeric(at_entry[[i]]), stringsAsFactors = FALSE
      )
    }
  }
  long <- do.call(rbind, rows[seq_len(k)])
  long <- long[long$stop > long$start, , drop = FALSE]
  long$sex <- factor(long$sex, levels = c("Male", "Female"))

  carry <- setdiff(
    intersect(c("trauma_type_count"), names(d)),
    names(long)
  )
  if (length(carry)) {
    carry_data <- unique(d[, c("record_id", carry), drop = FALSE])
    long <- merge(long, carry_data, by = "record_id", all.x = TRUE, sort = FALSE)
  }
  attr(long, "audit") <- data.frame(
    Participants = length(unique(long$record_id)),
    Interval_rows = nrow(long),
    Events = sum(long$event),
    Participants_switching_exposure = length(unique(long$record_id[long$new_worst_trauma_tv == 1])),
    Person_years_unexposed = sum((long$stop - long$start)[long$new_worst_trauma_tv == 0]),
    Person_years_exposed = sum((long$stop - long$start)[long$new_worst_trauma_tv == 1]),
    stringsAsFactors = FALSE
  )
  long
}

# Core model fitting moved to the individual RQ scripts for transparency.

make_origin_data <- function(data, origin) {
  if (origin == "Study entry") {
    out <- data.frame(
      record_id = data$record_id,
      sex = data$sex,
      event = data$event,
      entry = 0,
      duration = data$time_since_start / 12
    )
  } else if (origin == "Most recent trauma") {
    out <- data.frame(
      record_id = data$record_id,
      sex = data$sex,
      event = data$event,
      entry = data$months_since_most_recent_at_baseline / 12,
      duration = data$time_since_most_recent_trauma / 12
    )
  } else if (origin == "Worst trauma") {
    out <- data.frame(
      record_id = data$record_id,
      sex = data$sex,
      event = data$event,
      entry = data$months_since_worst_at_baseline / 12,
      duration = data$time_since_worst_trauma / 12
    )
  } else {
    stop("Unknown time origin: ", origin, call. = FALSE)
  }
  out <- out[complete.cases(out[, c("event", "entry", "duration")]), , drop = FALSE]
  out <- out[out$event %in% c(0, 1) & out$entry >= 0 & out$duration > out$entry, , drop = FALSE]
  out
}

km_fit <- function(data) {
  survfit(Surv(entry, duration, event) ~ 1, data = data, conf.type = "log-log")
}

km_risk_at_times <- function(fit, data, times = REPORT_TIMES, origin = "", group = "Overall") {
  s <- summary(fit, times = times, extend = TRUE)
  data.frame(
    Time_origin = origin,
    Group = group,
    Time_years = times,
    Risk = 1 - s$surv,
    CI_lower = 1 - s$upper,
    CI_upper = 1 - s$lower,
    N_at_risk = s$n.risk,
    Events_by_time = vapply(times, function(t) sum(data$event == 1 & data$duration <= t), integer(1)),
    Analysis_N = nrow(data),
    Total_events = sum(data$event),
    Delayed_entries = sum(data$entry > 0),
    stringsAsFactors = FALSE
  )
}

km_curve_table <- function(fit, origin = "", group = "Overall") {
  data.frame(
    Time_origin = origin,
    Group = group,
    Time_years = fit$time,
    Risk = 1 - fit$surv,
    CI_lower = 1 - fit$upper,
    CI_upper = 1 - fit$lower,
    N_at_risk = fit$n.risk,
    N_events_at_time = fit$n.event,
    N_censored_at_time = fit$n.censor,
    stringsAsFactors = FALSE
  )
}

km_numbers_at_risk <- function(data, times = seq(0, 2, by = .5), origin = "", group = "Overall") {
  data.frame(
    Time_origin = origin,
    Group = group,
    Time_years = times,
    N_at_risk = vapply(times, function(t) sum(data$entry <= t & data$duration >= t), integer(1)),
    stringsAsFactors = FALSE
  )
}

km_time_to_risk <- function(fit, target = .10, origin = "", group = "Overall") {
  idx <- which((1 - fit$surv) >= target)
  data.frame(
    Time_origin = origin,
    Group = group,
    Target_risk = target,
    Time_years = if (length(idx)) fit$time[min(idx)] else NA_real_,
    Threshold_reached = length(idx) > 0L,
    stringsAsFactors = FALSE
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

# The production bootstrap_km() implementation is defined below using a direct
# delayed-entry risk-set calculation. The earlier survfit-based implementation
# was removed to avoid signature/version ambiguity.

percentile_ci <- function(x, probs = c(.025, .975)) {
  x <- x[is.finite(x)]
  if (length(x) < 20L) return(c(NA_real_, NA_real_))
  unname(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 6))
}

latency_descriptives <- function(x, scenario, origin) {
  x <- x[is.finite(x)]
  mode_value <- NA_real_
  if (length(x)) {
    tab <- table(round(x, 6))
    mode_value <- as.numeric(names(tab)[which.max(tab)])
  }
  data.frame(
    Scenario = scenario,
    Time_origin = origin,
    N_events = length(x),
    Mean_years = if (length(x)) mean(x) else NA_real_,
    SD_years = if (length(x) > 1L) sd(x) else NA_real_,
    Median_years = if (length(x)) median(x) else NA_real_,
    Mode_years = mode_value,
    Minimum_years = if (length(x)) min(x) else NA_real_,
    Maximum_years = if (length(x)) max(x) else NA_real_,
    stringsAsFactors = FALSE
  )
}
build_midpoint_dataset <- function(data) {
  require_columns(
    data,
    c(
      "event",
      "age_entry",
      "age_exit",
      "first_positive_months",
      "last_negative_before_first_positive_months"
    ),
    "Midpoint sensitivity"
  )
  
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
  out$age_exit_midpoint[valid] <- out$age_entry[valid] + out$midpoint_months[valid] / 12
  
  invalid_midpoint_age <- valid & (
    !is.finite(out$age_exit_midpoint) |
      out$age_exit_midpoint <= out$age_entry |
      out$age_exit_midpoint > out$age_exit
  )
  
  valid[invalid_midpoint_age] <- FALSE
  
  out$midpoint_months[event_case & !valid] <- NA_real_
  out$age_exit_midpoint[event_case & !valid] <- NA_real_
  
  duration_cols <- c(
    "time_since_start",
    "time_since_most_recent_trauma",
    "time_since_worst_trauma"
  )
  
  for (col in duration_cols[duration_cols %in% names(out)]) {
    new_col <- paste0(col, "_midpoint")
    out[[new_col]] <- out[[col]]
    
    shifted <- out[[col]] - shift
    invalid_shifted <- valid & (!is.finite(shifted) | shifted < 0)
    
    out[[new_col]][valid] <- shifted[valid]
    out[[new_col]][invalid_shifted] <- NA_real_
    out[[new_col]][event_case & !valid] <- NA_real_
  }
  
  attr(out, "valid_midpoint") <- valid
  attr(out, "midpoint_audit") <- data.frame(
    Events = sum(event_case, na.rm = TRUE),
    Valid_midpoint_events = sum(valid, na.rm = TRUE),
    Event_cases_without_valid_midpoint = sum(event_case & !valid, na.rm = TRUE),
    Invalid_midpoint_age = sum(invalid_midpoint_age, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  out
}

build_cutoff_long <- function(data, trauma_age_col, exit_col = "age_exit", cutoff_months = 6L) {
  required <- c("record_id", "sex", "event", "age_entry", exit_col, trauma_age_col)
  require_columns(data, required, "Cutoff model")
  d <- complete_model_data(data, required)
  d <- d[d[[exit_col]] > d$age_entry & d$event %in% c(0, 1), , drop = FALSE]
  rows <- vector("list", nrow(d) * 2L)
  k <- 0L
  for (i in seq_len(nrow(d))) {
    start <- d$age_entry[[i]]
    stop <- d[[exit_col]][[i]]
    cut_age <- d[[trauma_age_col]][[i]] + cutoff_months / 12
    if (start < cut_age) {
      k <- k + 1L
      first_stop <- min(stop, cut_age)
      rows[[k]] <- data.frame(
        record_id = d$record_id[[i]], sex = as.character(d$sex[[i]]),
        start = start, stop = first_stop,
        event = as.integer(d$event[[i]] == 1 && stop <= cut_age),
        within_cutoff = 1,
        age_at_trauma = d[[trauma_age_col]][[i]],
        stringsAsFactors = FALSE
      )
    }
    if (stop > cut_age) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        record_id = d$record_id[[i]], sex = as.character(d$sex[[i]]),
        start = max(start, cut_age), stop = stop,
        event = as.integer(d$event[[i]] == 1),
        within_cutoff = 0,
        age_at_trauma = d[[trauma_age_col]][[i]],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows[seq_len(k)])
  out <- out[out$stop > out$start, , drop = FALSE]
  out$sex <- factor(as.character(out$sex), levels = c("Male", "Female"))
  out
}

# Core model fitting moved to the individual RQ scripts for transparency.

screen_binary_exposures <- function(
  data,
  variables,
  max_missing_prop = .20,
  min_exposed_n = 10L,
  min_unexposed_n = 10L,
  min_events_per_level = 3L
) {
  rows <- lapply(variables, function(variable) {
    if (!(variable %in% names(data))) {
      return(data.frame(
        Variable = variable, Available_N = 0L, Missing_prop = 1,
        Exposed_N = 0L, Unexposed_N = 0L, Exposed_events = 0L, Unexposed_events = 0L,
        Eligible = FALSE, Reason = "Variable absent", stringsAsFactors = FALSE
      ))
    }
    x <- suppressWarnings(as.numeric(data[[variable]]))
    observed <- is.finite(x) & x %in% c(0, 1) & data$event %in% c(0, 1)
    exposed <- observed & x == 1
    unexposed <- observed & x == 0
    missing_prop <- 1 - mean(observed)
    reason <- character()
    if (missing_prop > max_missing_prop) reason <- c(reason, "excess missingness")
    if (sum(exposed) < min_exposed_n) reason <- c(reason, "too few exposed")
    if (sum(unexposed) < min_unexposed_n) reason <- c(reason, "too few unexposed")
    if (sum(data$event[exposed]) < min_events_per_level) reason <- c(reason, "too few exposed events")
    if (sum(data$event[unexposed]) < min_events_per_level) reason <- c(reason, "too few unexposed events")
    data.frame(
      Variable = variable,
      Available_N = sum(observed),
      Missing_prop = missing_prop,
      Exposed_N = sum(exposed),
      Unexposed_N = sum(unexposed),
      Exposed_events = sum(data$event[exposed]),
      Unexposed_events = sum(data$event[unexposed]),
      Eligible = !length(reason),
      Reason = if (length(reason)) paste(reason, collapse = "; ") else "Eligible",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_apa_text <- function(title, inventory, model_terms = data.frame(), model_stats = data.frame(), ph = data.frame(), additional_lines = character()) {
  lines <- c(title, strrep("=", nchar(title)), "", "Analyses included:")
  lines <- c(lines, sprintf("%d. %s", seq_along(inventory), inventory), "")

  if (nrow(model_stats)) {
    lines <- c(lines, "MODEL-LEVEL RESULTS")
    for (i in seq_len(nrow(model_stats))) {
      x <- model_stats[i, ]
      lines <- c(lines, sprintf(
        "%s: N = %s, events = %s, likelihood-ratio chi-square(%s) = %s, p %s; Wald chi-square(%s) = %s, p %s; score chi-square(%s) = %s, p %s; concordance = %s (SE = %s); AIC = %s.",
        x$Model, x$N, x$Events,
        fmt_num(x$LR_df, 0), fmt_num(x$LR_chisq), fmt_p(x$LR_p),
        fmt_num(x$Wald_df, 0), fmt_num(x$Wald_chisq), fmt_p(x$Wald_p),
        fmt_num(x$Score_df, 0), fmt_num(x$Score_chisq), fmt_p(x$Score_p),
        fmt_num(x$Concordance, 3), fmt_num(x$Concordance_SE, 3), fmt_num(x$AIC)
      ))
    }
    lines <- c(lines, "")
  }

  if (nrow(model_terms)) {
    lines <- c(lines, "COEFFICIENT RESULTS")
    for (i in seq_len(nrow(model_terms))) {
      x <- model_terms[i, ]
      lines <- c(lines, sprintf(
        "%s, %s: HR = %s, 95%% CI [%s, %s], B = %s, SE = %s, z = %s, p %s (N = %s; events = %s).",
        x$Model, x$Term, fmt_num(x$HR), fmt_num(x$CI_lower), fmt_num(x$CI_upper),
        fmt_num(x$B, 3), fmt_num(x$SE, 3), fmt_num(x$z), fmt_p(x$p), x$N, x$Events
      ))
    }
    lines <- c(lines, "")
  }

  if (nrow(ph)) {
    lines <- c(lines, "PROPORTIONAL-HAZARDS DIAGNOSTICS")
    for (i in seq_len(nrow(ph))) {
      x <- ph[i, ]
      if (x$PH_term == "ERROR") {
        err <- if ("PH_error" %in% names(ph)) x$PH_error else "Unknown error"
        lines <- c(lines, sprintf("%s: cox.zph failed: %s", x$Model, err))
      } else {
        lines <- c(lines, sprintf(
          "%s, %s: rank-transformed Schoenfeld test chi-square(%s) = %s, p %s; PH concern = %s. %s.",
          x$Model, x$PH_term, fmt_num(x$PH_df, 0), fmt_num(x$PH_chisq), fmt_p(x$PH_p),
          ifelse(isTRUE(x$PH_concern), "yes", "no"), x$Planned_action
        ))
      }
    }
    lines <- c(lines, "")
  }

  c(lines, additional_lines)
}

save_apa_bundle <- function(
  output_dir,
  prefix,
  title,
  inventory,
  model_terms = data.frame(),
  model_stats = data.frame(),
  ph = data.frame(),
  extra_tables = list(),
  additional_lines = character()
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (nrow(model_terms)) write_csv_safe(model_terms, file.path(output_dir, paste0(prefix, "_coefficient_results.csv")))
  if (nrow(model_stats)) write_csv_safe(model_stats, file.path(output_dir, paste0(prefix, "_model_statistics.csv")))
  if (nrow(ph)) write_csv_safe(ph, file.path(output_dir, paste0(prefix, "_PH_diagnostics.csv")))
  for (nm in names(extra_tables)) {
    write_csv_safe(extra_tables[[nm]], file.path(output_dir, paste0(prefix, "_", nm, ".csv")))
  }
  lines <- make_apa_text(title, inventory, model_terms, model_stats, ph, additional_lines)
  txt_path <- file.path(output_dir, paste0(prefix, "_APA_results_summary.txt"))
  writeLines(lines, txt_path, useBytes = TRUE)
  message("Saved: ", txt_path)
  invisible(txt_path)
}


# =============================================================================
# Fast execution, progress reporting, and publication-style figures
# =============================================================================

SURVIVAL_CORES <- suppressWarnings(as.integer(Sys.getenv(
  "SURVIVAL_CORES",
  if (.Platform$OS.type == "windows") "1" else as.character(max(1L, min(4L, parallel::detectCores(logical = TRUE) - 1L)))
)))
if (!is.finite(SURVIVAL_CORES) || SURVIVAL_CORES < 1L) SURVIVAL_CORES <- 1L
FIGURE_FONT <- Sys.getenv("SURVIVAL_FIGURE_FONT", "serif")
SET2_PALETTE <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
SET2_COLORS <- setNames(SET2_PALETTE[c(1L, 5L, 8L)], c("Overall", "Female", "Male"))

open_figure_device <- function(path, width_cm, height_cm, type = c("png", "svg")) {
  type <- match.arg(type)
  width_in <- width_cm / 2.54
  height_in <- height_cm / 2.54
  if (type == "png") {
    grDevices::png(path, width = width_in, height = height_in, units = "in", res = 600, pointsize = 9, bg = "white")
  } else {
    grDevices::svg(path, width = width_in, height = height_in, pointsize = 9, family = FIGURE_FONT, onefile = FALSE)
  }
}

with_figure_devices <- function(path_base, width_cm, height_cm, draw) {
  outputs <- character()
  expected <- paste0(path_base, c(".png", ".svg"))
  for (ext in c("png", "svg")) {
    path <- paste0(path_base, ".", ext)
    if (file.exists(path)) unlink(path)
    ok <- tryCatch({
      open_figure_device(path, width_cm, height_cm, ext)
      on.exit(grDevices::dev.off(), add = TRUE)
      draw()
      grDevices::dev.off()
      on.exit(NULL, add = FALSE)
      TRUE
    }, error = function(e) {
      if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
      progress_log("Figure error", sprintf("%s: %s", path, conditionMessage(e)))
      FALSE
    })
    valid_file <- isTRUE(ok) && file.exists(path) && is.finite(file.info(path)$size) && file.info(path)$size > 0
    if (valid_file) {
      outputs <- c(outputs, path)
      progress_log("Figure saved", sprintf("%s [%.1f x %.1f cm]", normalizePath(path, mustWork = FALSE), width_cm, height_cm))
    }
  }
  missing <- expected[!file.exists(expected) | is.na(file.info(expected)$size) | file.info(expected)$size <= 0]
  if (length(missing)) {
    stop("Required figure output was not saved: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(outputs)
}

add_minor_ticks <- function(side, at, tcl = -0.18) {
  axis(side, at = at, labels = FALSE, tcl = tcl)
}

# Fast delayed-entry Kaplan-Meier summaries used only inside bootstrap loops.
# The published curve and pointwise confidence limits still come from survival::survfit.
fast_delayed_km <- function(entry, duration, event, times = numeric(), tau = NULL) {
  keep <- is.finite(entry) & is.finite(duration) & event %in% c(0, 1) & entry >= 0 & duration > entry
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
    survival_step <- ifelse(n_risk > 0L, pmax(0, 1 - n_event / n_risk), 1)
    survival_after <- cumprod(survival_step)
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

bootstrap_chunk_fast <- function(data, times, tau, repetitions, seed) {
  set.seed(seed)
  n <- nrow(data)
  risk_mat <- matrix(NA_real_, nrow = repetitions, ncol = length(times))
  rmst_vec <- rep(NA_real_, repetitions)
  entry <- data$entry
  duration <- data$duration
  event <- data$event
  for (b in seq_len(repetitions)) {
    idx <- sample.int(n, n, replace = TRUE)
    estimate <- fast_delayed_km(entry[idx], duration[idx], event[idx], times = times, tau = tau)
    if (length(times)) risk_mat[b, ] <- estimate$risk
    if (!is.null(tau)) rmst_vec[b] <- estimate$rmst
  }
  list(risk = risk_mat, rmst = rmst_vec)
}

# Overrides the earlier bootstrap helper with a faster direct risk-set implementation.
bootstrap_km <- function(data, times = REPORT_TIMES, tau = NULL, B = BOOT_B, seed = RNG_SEED, label = NULL) {
  if (B <= 0L) return(list(risk = matrix(numeric(), 0, length(times)), rmst = numeric()))
  label <- label %||% sprintf("KM bootstrap; N=%d", nrow(data))
  workers <- min(SURVIVAL_CORES, B)
  if (.Platform$OS.type == "windows") workers <- 1L
  started <- Sys.time()
  progress_log("Bootstrap", sprintf("starting %s; B=%d, cores=%d", label, B, workers))

  if (workers <= 1L) {
    set.seed(seed)
    risk_mat <- matrix(NA_real_, nrow = B, ncol = length(times))
    rmst_vec <- rep(NA_real_, B)
    n <- nrow(data)
    step <- max(1L, floor(B / 10L))
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      estimate <- fast_delayed_km(data$entry[idx], data$duration[idx], data$event[idx], times = times, tau = tau)
      if (length(times)) risk_mat[b, ] <- estimate$risk
      if (!is.null(tau)) rmst_vec[b] <- estimate$rmst
      if (b %% step == 0L || b == B) progress_log("Bootstrap progress", sprintf("%s: %d/%d", label, b, B))
    }
    result <- list(risk = risk_mat, rmst = rmst_vec)
  } else {
    counts <- rep(B %/% workers, workers)
    counts[seq_len(B %% workers)] <- counts[seq_len(B %% workers)] + 1L
    counts <- counts[counts > 0L]
    chunks <- parallel::mclapply(
      seq_along(counts),
      function(i) bootstrap_chunk_fast(data, times, tau, counts[[i]], seed + 100003L * i),
      mc.cores = length(counts),
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
    risk_mat <- if (length(times)) do.call(rbind, lapply(chunks, `[[`, "risk")) else matrix(numeric(), nrow = B, ncol = 0L)
    rmst_vec <- unlist(lapply(chunks, `[[`, "rmst"), use.names = FALSE)
    result <- list(risk = risk_mat, rmst = rmst_vec)
  }

  progress_log("Bootstrap", sprintf("completed %s; elapsed=%s", label, format_elapsed(started)))
  result
}

step_values <- function(time, value, x_grid) {
  keep <- is.finite(time) & is.finite(value)
  time <- time[keep]
  value <- value[keep]
  if (!length(time)) return(rep(0, length(x_grid)))
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  # Retain the last value at duplicated times.
  last <- !duplicated(time, fromLast = TRUE)
  time <- time[last]
  value <- value[last]
  fn <- stats::stepfun(time, c(0, value), right = TRUE)
  as.numeric(fn(x_grid))
}

# Manuscript Figure 1-3 rendering was removed from the shared R helper.
# R scripts save figure-ready CSV files; 05_make_figures.py renders them with Matplotlib.


# =============================================================================
# Compatibility wrappers and additional requested helper functions
# =============================================================================

# Requested alias for the RQ1 threshold summary.
rq1_time_to_risk <- function(fit, target = .10, origin = "", group = "Overall") {
  km_time_to_risk(fit, target = target, origin = origin, group = group)
}

extract_hr_table <- function(fit, analysis_name = attr(fit, "analysis_name") %||% attr(fit, "model_name") %||% "Cox model") {
  term_df <- extract_cox_terms(fit, model_name = analysis_name)
  stats_df <- extract_model_stats(fit, model_name = analysis_name)
  if (!nrow(term_df)) {
    return(data.frame())
  }
  term_df$analysis_name <- analysis_name
  term_df$formula <- stats_df$Formula[1]
  term_df$N_model <- stats_df$N[1]
  term_df$Events_model <- stats_df$Events[1]
  term_df$concordance <- stats_df$Concordance[1]
  term_df$AIC <- stats_df$AIC[1]
  term_df$log_likelihood <- stats_df$Log_likelihood_model[1]
  term_df
}

save_model_results <- function(
  fit,
  output_dir,
  prefix = attr(fit, "model_name") %||% "model",
  analysis_name = attr(fit, "analysis_name") %||% attr(fit, "model_name") %||% "Cox model",
  ph_transform = "rank",
  ph_main_exposures = character(),
  ph_nuisance_terms = character(),
  extra_tables = list()
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  term_table <- extract_hr_table(fit, analysis_name = analysis_name)
  model_stats <- extract_model_stats(fit, model_name = analysis_name)
  ph_results <- tryCatch(
    run_ph_test(fit, transform = ph_transform, main_exposures = ph_main_exposures, nuisance_terms = ph_nuisance_terms),
    error = function(e) NULL
  )
  ph_table <- if (is.list(ph_results) && !is.null(ph_results$table)) ph_results$table else data.frame()

  term_path <- file.path(output_dir, paste0(prefix, "_coefficient_results.csv"))
  stats_path <- file.path(output_dir, paste0(prefix, "_model_statistics.csv"))
  ph_path <- file.path(output_dir, paste0(prefix, "_PH_diagnostics.csv"))
  write_csv_safe(term_table, term_path)
  write_csv_safe(model_stats, stats_path)
  if (nrow(ph_table)) write_csv_safe(ph_table, ph_path)
  for (nm in names(extra_tables)) {
    write_csv_safe(extra_tables[[nm]], file.path(output_dir, paste0(prefix, "_", nm, ".csv")))
  }
  invisible(list(term_table = term_table, model_stats = model_stats, ph_table = ph_table,
                 term_path = term_path, stats_path = stats_path, ph_path = ph_path))
}

run_ph_test <- function(
  fit,
  transform = "rank",
  terms = TRUE,
  singledf = FALSE,
  main_exposures = character(),
  nuisance_terms = character()
) {
  zph <- survival::cox.zph(fit, transform = transform, terms = terms, singledf = singledf)
  table <- extract_ph_results(zph, model_name = attr(fit, "model_name") %||% attr(fit, "analysis_name") %||% "Cox model")
  structure(
    list(object = zph, table = table, main_exposures = main_exposures, nuisance_terms = nuisance_terms),
    class = "ph_test_result"
  )
}

extract_ph_results <- function(ph_object, model_name = NULL) {
  zph <- if (inherits(ph_object, "ph_test_result") && !is.null(ph_object$object)) {
    ph_object$object
  } else {
    ph_object
  }
  model_name <- model_name %||% attr(ph_object, "model_name") %||% "Cox model"

  if (inherits(zph, "error")) {
    return(data.frame(
      Model = model_name,
      PH_term = "ERROR",
      PH_chisq = NA_real_,
      PH_df = NA_real_,
      PH_p = NA_real_,
      PH_concern = NA,
      PH_transform = NA_character_,
      PH_error = conditionMessage(zph),
      stringsAsFactors = FALSE
    ))
  }

  tab <- as.data.frame(zph$table, check.names = FALSE)
  if (!nrow(tab)) {
    return(data.frame(
      Model = model_name,
      PH_term = "NONE",
      PH_chisq = NA_real_,
      PH_df = NA_real_,
      PH_p = NA_real_,
      PH_concern = NA,
      PH_transform = attr(zph, "transform") %||% NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  out <- data.frame(
    Model = model_name,
    PH_term = rownames(tab),
    PH_chisq = tab[["chisq"]],
    PH_df = tab[["df"]],
    PH_p = tab[["p"]],
    PH_concern = tab[["p"]] < PH_ALPHA,
    PH_transform = attr(zph, "transform") %||% NA_character_,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out$Planned_action <- "No action unless supported by residual pattern and scientific role"
  out$Planned_action[out$PH_term == "GLOBAL"] <- "Review variable-level diagnostics and residual plots"
  out$Planned_action[out$PH_concern %in% TRUE & out$PH_term != "GLOBAL"] <- "Consider extended Cox model or clinically coherent stratification"
  out
}

save_descriptive_table <- function(
  data,
  variables,
  output_path = NULL,
  group_var = NULL,
  label = "Descriptive table"
) {
  if (!is.data.frame(data)) stop("save_descriptive_table() requires a data frame.", call. = FALSE)
  if (missing(variables) || !length(variables)) variables <- names(data)
  variables <- variables[variables %in% names(data)]
  if (!length(variables)) stop("No requested variables were found in the data.", call. = FALSE)

  summarise_one <- function(df, var_name, group_value = NA_character_) {
    x <- df[[var_name]]
    n_total <- length(x)
    n_missing <- sum(is.na(x))
    non_missing <- x[!is.na(x)]
    is_numeric_like <- is.numeric(non_missing) || is.integer(non_missing) || is.logical(non_missing)
    binary_numeric <- is_numeric_like && length(unique(non_missing)) <= 2L && all(unique(non_missing) %in% c(0, 1, FALSE, TRUE, NA))
    if (is_numeric_like && !binary_numeric) {
      numeric_x <- suppressWarnings(as.numeric(non_missing))
      data.frame(
        Label = label,
        Group = group_value,
        Variable = var_name,
        Variable_type = "continuous",
        Level = NA_character_,
        N = n_total,
        Missing = n_missing,
        Mean = if (length(numeric_x)) mean(numeric_x) else NA_real_,
        Median = if (length(numeric_x)) median(numeric_x) else NA_real_,
        SD = if (length(numeric_x) > 1L) sd(numeric_x) else NA_real_,
        Min = if (length(numeric_x)) min(numeric_x) else NA_real_,
        Max = if (length(numeric_x)) max(numeric_x) else NA_real_,
        Count = length(numeric_x),
        Percent = if (n_total > 0L) 100 * length(numeric_x) / n_total else NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      if (binary_numeric) {
        f <- factor(as.character(non_missing), levels = c("0", "1", "FALSE", "TRUE"))
        f <- droplevels(f)
      } else {
        f <- factor(non_missing)
      }
      lvls <- levels(f)
      if (!length(lvls)) lvls <- character(0)
      counts <- if (length(lvls)) as.integer(tabulate(as.integer(f), nbins = length(lvls))) else integer(0)
      if (!length(lvls)) {
        data.frame(
          Label = label,
          Group = group_value,
          Variable = var_name,
          Variable_type = "categorical",
          Level = NA_character_,
          N = n_total,
          Missing = n_missing,
          Mean = NA_real_, Median = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_,
          Count = 0L, Percent = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          Label = label,
          Group = group_value,
          Variable = var_name,
          Variable_type = "categorical",
          Level = lvls,
          N = n_total,
          Missing = n_missing,
          Mean = NA_real_, Median = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_,
          Count = counts,
          Percent = if (n_total > 0L) 100 * counts / n_total else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  pieces <- list()
  if (!is.null(group_var) && nzchar(group_var)) {
    if (!(group_var %in% names(data))) stop("Grouping variable not found: ", group_var, call. = FALSE)
    groups <- unique(data[[group_var]])
    groups <- groups[!is.na(groups)]
    for (g in groups) {
      df <- data[data[[group_var]] == g, , drop = FALSE]
      for (v in variables) pieces[[length(pieces) + 1L]] <- summarise_one(df, v, as.character(g))
    }
  } else {
    for (v in variables) pieces[[length(pieces) + 1L]] <- summarise_one(data, v, NA_character_)
  }
  out <- rbind_fill(pieces)
  if (!is.null(output_path) && nzchar(output_path)) write_csv_safe(out, output_path)
  out
}
write_sample_summary_apa_txt <- function(
    data,
    out_file,
    age_col_candidates = c("age_entry", "age_at_study_start", "baseline_age", "age_start"),
    sex_col_candidates = c("sex", "gender"),
    baseline_ocd_candidates = c("baseline_ocd", "baseline_OCD", "ocd_baseline", "ocd_at_baseline"),
    baseline_only_candidates = c("baseline_only", "baseline_data_only"),
    sex_missing_candidates = c("sex_missing", "missing_sex"),
    followup_col_candidates = c("followup_years", "time_since_start", "duration_years", "duration"),
    recent_trauma_col_candidates = c("time_since_most_recent_trauma", "months_since_most_recent_trauma", "time_since_latest_trauma"),
    worst_trauma_col_candidates = c("time_since_worst_trauma", "months_since_worst_trauma"),
    recent_baseline_col_candidates = c("months_since_most_recent_at_baseline", "months_since_latest_at_baseline"),
    worst_baseline_col_candidates = c("months_since_worst_at_baseline")
) {
  stopifnot(is.data.frame(data), is.character(out_file), length(out_file) == 1L)
  
  pick_col <- function(cands) {
    hit <- cands[cands %in% names(data)]
    if (length(hit)) hit[[1L]] else NA_character_
  }
  
  num <- function(x, digits = 2L) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
  nint <- function(x) ifelse(is.na(x), "NA", formatC(round(x), format = "f", digits = 0L))
  smean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  ssd <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
  smed <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  smin <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  smax <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  
  age_col <- pick_col(age_col_candidates)
  sex_col <- pick_col(sex_col_candidates)
  if (is.na(age_col) || is.na(sex_col)) {
    stop("Need at least an age column and a sex column.", call. = FALSE)
  }
  
  baseline_ocd_col <- pick_col(baseline_ocd_candidates)
  baseline_only_col <- pick_col(baseline_only_candidates)
  sex_missing_col <- pick_col(sex_missing_candidates)
  followup_col <- pick_col(followup_col_candidates)
  recent_trauma_col <- pick_col(recent_trauma_col_candidates)
  worst_trauma_col <- pick_col(worst_trauma_col_candidates)
  recent_baseline_col <- pick_col(recent_baseline_col_candidates)
  worst_baseline_col <- pick_col(worst_baseline_col_candidates)
  
  d <- data
  for (nm in na.omit(c(age_col, followup_col, recent_trauma_col, worst_trauma_col,
                       recent_baseline_col, worst_baseline_col,
                       baseline_ocd_col, baseline_only_col, sex_missing_col))) {
    d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  }
  
  exclude <- rep(FALSE, nrow(d))
  if (!is.na(baseline_ocd_col)) exclude <- exclude | d[[baseline_ocd_col]] == 1
  if (!is.na(baseline_only_col)) exclude <- exclude | d[[baseline_only_col]] == 1
  if (!is.na(sex_missing_col)) exclude <- exclude | d[[sex_missing_col]] == 1
  exclude <- exclude | is.na(d[[sex_col]])
  
  final <- d[!exclude, , drop = FALSE]
  
  sex_tab <- table(final[[sex_col]], useNA = "no")
  female_n <- if ("Female" %in% names(sex_tab)) unname(sex_tab[["Female"]]) else NA_integer_
  male_n <- if ("Male" %in% names(sex_tab)) unname(sex_tab[["Male"]]) else NA_integer_
  sex_total <- sum(sex_tab)
  
  lines <- c(
    "Results",
    "",
    sprintf(
      "The analytic sample included %s trauma-exposed youth (mean age at study start = %s years, SD = %s).",
      nint(nrow(d)),
      num(smean(d[[age_col]])),
      num(ssd(d[[age_col]]))
    ),
    if (!is.na(baseline_ocd_col)) {
      sprintf(
        "At baseline, %s participants (%s%%) met criteria for OCD and were excluded from prospective analyses.",
        nint(sum(d[[baseline_ocd_col]] == 1, na.rm = TRUE)),
        num(100 * mean(d[[baseline_ocd_col]] == 1, na.rm = TRUE), 1)
      )
    } else {
      "At baseline, participants meeting OCD criteria were excluded from prospective analyses."
    },
    if (!is.na(baseline_only_col)) {
      sprintf("An additional %s participants had only baseline data.", nint(sum(d[[baseline_only_col]] == 1, na.rm = TRUE)))
    } else {
      "An additional subset had only baseline data."
    },
    if (!is.na(sex_missing_col)) {
      sprintf(
        "%s participants were missing information on sex, leaving %s youth for longitudinal follow-up.",
        nint(sum(d[[sex_missing_col]] == 1, na.rm = TRUE)),
        nint(nrow(final))
      )
    } else {
      sprintf("After excluding participants with missing sex information, %s youth remained for longitudinal follow-up.", nint(nrow(final)))
    },
    sprintf(
      "The final longitudinal sample had a baseline mean age of %s years (SD = %s, median = %s, range = %s\u2013%s).",
      num(smean(final[[age_col]])),
      num(ssd(final[[age_col]])),
      num(smed(final[[age_col]])),
      num(smin(final[[age_col]])),
      num(smax(final[[age_col]]))
    ),
    if (length(sex_tab)) {
      sprintf(
        "By sex, %s participants (%s%%) were female and %s (%s%%) were male.",
        nint(female_n),
        num(100 * female_n / sex_total, 1),
        nint(male_n),
        num(100 * male_n / sex_total, 1)
      )
    } else {
      "By sex, the distribution was not available."
    },
    if (!is.na(followup_col)) {
      sprintf(
        "Included participants were followed for a mean of %s years (SD = %s; median = %s, range = %s\u2013%s).",
        num(smean(final[[followup_col]])),
        num(ssd(final[[followup_col]])),
        num(smed(final[[followup_col]])),
        num(smin(final[[followup_col]])),
        num(smax(final[[followup_col]]))
      )
    } else {
      "Follow-up duration was not available."
    },
    if (!is.na(recent_trauma_col)) {
      sprintf(
        "Time since the most recent trauma had a mean of %s years (SD = %s; median = %s, range = %s\u2013%s).",
        num(smean(final[[recent_trauma_col]])),
        num(ssd(final[[recent_trauma_col]])),
        num(smed(final[[recent_trauma_col]])),
        num(smin(final[[recent_trauma_col]])),
        num(smax(final[[recent_trauma_col]]))
      )
    } else {
      "Time since the most recent trauma was not available."
    },
    if (!is.na(worst_trauma_col)) {
      sprintf(
        "Time since the worst trauma had a mean of %s years (SD = %s; median = %s, range = %s\u2013%s).",
        num(smean(final[[worst_trauma_col]])),
        num(ssd(final[[worst_trauma_col]])),
        num(smed(final[[worst_trauma_col]])),
        num(smin(final[[worst_trauma_col]])),
        num(smax(final[[worst_trauma_col]]))
      )
    } else {
      "Time since the worst trauma was not available."
    },
    if (!is.na(recent_baseline_col)) {
      sprintf(
        "At baseline, participants had experienced their most recent trauma a median of %s years earlier (mean = %s, SD = %s; n = %s).",
        num(smed(final[[recent_baseline_col]])),
        num(smean(final[[recent_baseline_col]])),
        num(ssd(final[[recent_baseline_col]])),
        nint(sum(is.finite(final[[recent_baseline_col]])))
      )
    } else {
      "Baseline timing for most recent trauma was not available."
    },
    if (!is.na(worst_baseline_col)) {
      sprintf(
        "At baseline, participants had experienced their worst trauma a median of %s years earlier (mean = %s, SD = %s; n = %s).",
        num(smed(final[[worst_baseline_col]])),
        num(smean(final[[worst_baseline_col]])),
        num(ssd(final[[worst_baseline_col]])),
        nint(sum(is.finite(final[[worst_baseline_col]])))
      )
    } else {
      "Baseline timing for worst trauma was not available."
    }
  )
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, con = out_file, useBytes = TRUE)
  invisible(normalizePath(out_file, mustWork = FALSE))
}
