# Shared IPCW Cox + figures for trajectory variants (exploratory).
# Requires: source("00_constants.R") in the caller so DIR_* and PATH_* exist.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(tibble)
  library(rlang)
})

#' @param traj_col Name of logical column in feat: TRUE = deteriorating, FALSE = stable
#' @param suffix Output file suffix, e.g. "tight", "equipoise", "equipoise_tight"
#' @param equipoise_only Restrict to equipoise (SOFA 8–14) cohort
#' @param title_short Short label for plot subtitles (ignored when legacy_main_03)
#' @param legacy_main_03 If TRUE, write the original 03 filenames and titles (cox_hr_by_trajectory_stratum.csv, etc.)
run_trajectory_analysis <- function(
    feat,
    baseline,
    ccw,
    traj_col,
    suffix,
    equipoise_only = FALSE,
    title_short = "",
    legacy_main_03 = FALSE
) {
  if (!traj_col %in% names(feat)) {
    stop("Column not found in feat: ", traj_col)
  }

  fj <- feat %>%
    dplyr::select(stay_id, traj = !!sym(traj_col)) %>%
    filter(!is.na(traj))

  ccw_tr <- ccw %>%
    distinct(clone_id, .keep_all = TRUE) %>%
    left_join(fj, by = "stay_id") %>%
    filter(!is.na(traj)) %>%
    mutate(
      trt = as.integer(strategy == "A_early"),
      det = as.integer(traj)
    )

  if (equipoise_only) {
    ccw_tr <- ccw_tr %>% filter(equipoise %in% TRUE)
  }

  n_pat <- n_distinct(ccw_tr$stay_id)
  log_path <- if (legacy_main_03) {
    file.path(DIR_LOG, "03_cox_output.txt")
  } else {
    file.path(DIR_LOG, paste0("cox_", suffix, ".txt"))
  }
  sink(log_path)
  if (legacy_main_03) {
    cat("Patients with non-NA trajectory label:", n_pat, "\n")
    cat("Clones:", nrow(ccw_tr), "\n\n")
  } else {
    cat("Variant:", suffix, "| traj_col:", traj_col, "| equipoise_only:", equipoise_only, "\n")
    cat("Patients with non-NA trajectory:", n_pat, "| Clones:", nrow(ccw_tr), "\n\n")
  }

  fit_int <- tryCatch(
    coxph(
      Surv(follow_days, event) ~ trt * det,
      data = ccw_tr,
      weights = ipcw_w_trim,
      cluster = stay_id,
      robust = TRUE
    ),
    error = function(e) {
      cat("Interaction model failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(fit_int)) {
    print(summary(fit_int))
  }

  strat_one <- function(is_det) {
    sub <- ccw_tr %>% filter(traj == is_det)
    if (nrow(sub) < 20L) {
      return(tibble(
        stratum = ifelse(is_det, "Deteriorating trajectory", "Stable trajectory"),
        hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, n_pat = 0L
      ))
    }
    fit <- coxph(
      Surv(follow_days, event) ~ trt,
      data = sub,
      weights = ipcw_w_trim,
      cluster = stay_id,
      robust = TRUE
    )
    s <- summary(fit)
    tibble(
      stratum = ifelse(is_det, "Deteriorating trajectory", "Stable trajectory"),
      hr = s$coef["trt", "exp(coef)"],
      lo = s$conf.int["trt", "lower .95"],
      hi = s$conf.int["trt", "upper .95"],
      p = s$coefficients["trt", "Pr(>|z|)"],
      n_pat = n_distinct(sub$stay_id)
    )
  }

  strat_res <- bind_rows(strat_one(FALSE), strat_one(TRUE))

  fit_pool <- coxph(
    Surv(follow_days, event) ~ trt,
    data = ccw_tr,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE
  )
  sp <- summary(fit_pool)
  overall_row <- tibble(
    stratum = "All (with trajectory label)",
    hr = sp$coef["trt", "exp(coef)"],
    lo = sp$conf.int["trt", "lower .95"],
    hi = sp$conf.int["trt", "upper .95"],
    p = sp$coefficients["trt", "Pr(>|z|)"],
    n_pat = n_distinct(ccw_tr$stay_id)
  )

  forest_df <- bind_rows(overall_row, strat_res)
  csv_name <- if (legacy_main_03) {
    "cox_hr_by_trajectory_stratum.csv"
  } else {
    paste0("cox_hr_", suffix, ".csv")
  }
  write_csv(forest_df, file.path(DIR_DATA, csv_name))

  if (!is.null(fit_int)) {
    int_name <- if (legacy_main_03) {
      "03_cox_interaction_full.txt"
    } else {
      paste0("cox_interaction_", suffix, ".txt")
    }
    capture.output(
      print(summary(fit_int)),
      file = file.path(DIR_LOG, int_name)
    )
  }

  sink()
  print(forest_df)

  forest_plot_df <- forest_df %>%
    filter(!is.na(hr)) %>%
    mutate(stratum = factor(stratum, levels = rev(unique(stratum))))

  xmax <- max(2.6, max(forest_plot_df$hi, na.rm = TRUE) * 1.05, na.rm = TRUE)
  brk_hi <- max(3, ceiling(xmax * 4) / 4)

  if (legacy_main_03) {
    ft_title <- "Exploratory: IPCW-weighted HR (Early vs Deferred) by 6h trajectory"
    ft_sub <- "Non-prespecified; clone-weighted Cox, cluster-robust SE"
    forest_png <- "forest_hr_by_trajectory.png"
    x_breaks <- seq(0.5, brk_hi, 0.25)
  } else {
    ft_title <- paste("Exploratory HR by trajectory:", suffix)
    ft_sub <- paste(
      "Non-prespecified | IPCW Cox, cluster-robust SE",
      if (nzchar(title_short)) title_short else suffix
    )
    forest_png <- paste0("forest_hr_", suffix, ".png")
    x_breaks <- seq(0.5, ceiling(xmax * 4) / 4, 0.25)
  }

  p_forest <- ggplot(forest_plot_df, aes(x = hr, y = stratum)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, color = "#2C3E50") +
    geom_point(size = 3, color = "#2980B9") +
    scale_x_continuous(limits = c(0.35, xmax), breaks = x_breaks) +
    labs(
      title = ft_title,
      subtitle = ft_sub,
      x = "Hazard ratio",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    file.path(DIR_FIG, forest_png),
    p_forest, width = 8, height = 4, dpi = 200
  )

  bl <- baseline %>%
    dplyr::select(stay_id, early_crrt) %>%
    distinct(stay_id, .keep_all = TRUE)

  cnt <- feat %>%
    filter(!is.na(!!sym(traj_col))) %>%
    inner_join(bl, by = "stay_id")

  if (equipoise_only) {
    cnt <- cnt %>% filter(sofa_equipoise %in% TRUE)
  }

  cnt <- cnt %>%
    count(!!sym(traj_col), early_crrt)

  p_bar <- ggplot(cnt, aes(
    x = ifelse(!!sym(traj_col), "Deteriorating", "Stable"),
    y = n,
    fill = factor(early_crrt, labels = c("Deferred", "Early"))
  )) +
    geom_col(position = "stack", width = 0.65) +
    scale_fill_manual(values = c("Deferred" = "#E67E22", "Early" = "#3498DB"), name = NULL) +
    labs(
      title = if (legacy_main_03) {
        "Exploratory: patients by 6h trajectory and treatment arm"
      } else {
        paste("N by trajectory:", suffix)
      },
      x = if (legacy_main_03) {
        "Trajectory (first 6h post-t0)"
      } else {
        "Trajectory (6h post-t0)"
      },
      y = "N patients"
    ) +
    theme_minimal(base_size = 12)

  bar_png <- if (legacy_main_03) {
    "bar_n_by_trajectory_treatment.png"
  } else {
    paste0("bar_n_", suffix, ".png")
  }
  ggsave(
    file.path(DIR_FIG, bar_png),
    p_bar, width = 7, height = 4, dpi = 200
  )

  if (legacy_main_03) {
    message("Figures written to ", DIR_FIG)
  } else {
    message("Saved variant: ", suffix, " → ", DIR_FIG)
  }
  invisible(list(ccw_tr = ccw_tr, forest_df = forest_df, fit_int = fit_int))
}
