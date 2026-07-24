# Exploratory: 6h post-t₀ trajectory — paths & thresholds (not part of main pre-spec)

Sys.setenv(LC_ALL = "en_US.UTF-8")

ROOT <- "analysis"
DIR_OUT <- file.path(ROOT, "exploratory_trajectory_6h", "output")
DIR_DATA <- file.path(DIR_OUT, "data")
DIR_FIG <- file.path(DIR_OUT, "figures")
DIR_LOG <- file.path(DIR_OUT, "logs")

for (d in c(DIR_OUT, DIR_DATA, DIR_FIG, DIR_LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

PATH_BASELINE <- file.path(ROOT, "data", "baseline.rds")
PATH_TV <- file.path(ROOT, "data", "tv_hourly.rds")
PATH_CCW <- file.path(ROOT, "data", "ccw_final.rds")

# Clinical thresholds (exploratory)
UO_THRESH_MLKGH <- 0.3
LACTATE_THRESH <- 4
NEE_INCREASE_EPS <- 1e-6

# Hour windows: hour_index 0 = [t0, t0+1h), …, 5 = [t0+5h, t0+6h)
HOURS_POST_T0 <- 0:5L
