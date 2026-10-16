library(data.table)
library(ggplot2)
source("R/palette.R")

#' Standard atmosphere pressure from elevation
#'
#' Uses the international barometric formula (ISA) to compute expected
#' mean sea-level-adjusted station pressure.
#'
#' @param elev numeric, elevation in metres above sea level
#' @return numeric, expected pressure in hPa
standard_pressure <- function(elev) {
  P0 <- 1013.25   # sea-level standard pressure (hPa)
  L  <- 0.0065    # temperature lapse rate (K/m)
  T0 <- 288.15    # sea-level standard temperature (K)
  g  <- 9.80665   # gravitational acceleration (m/s2)
  M  <- 0.0289644 # molar mass of dry air (kg/mol)
  R  <- 8.31447   # universal gas constant (J/(mol·K))
  P0 * (1 - L * elev / T0)^(g * M / (R * L))
}

#' Assign hydrological year (July 1 start)
#'
#' @param dates Date or POSIXct vector
#' @return integer vector of hydrological years (e.g. 2024 = Jul 2024 – Jun 2025)
hydro_year <- function(dates) {
  m <- as.integer(format(dates, "%m"))
  y <- as.integer(format(dates, "%Y"))
  ifelse(m >= 7, y, y - 1L)
}

#' Run barometric pressure QC for a single lake
#'
#' @param lake_name character
#' @param rds_dir character
#' @return data.table with DateTime, PrBaro, hydro_year, expected_hpa, lake
load_baro <- function(lake_name, rds_dir = "data/rds") {
  path <- file.path(rds_dir, paste0("full_", lake_name, ".rds"))
  x <- readRDS(path)

  meta <- x[["meta"]]
  elev <- as.numeric(meta$value[meta$attribute == "masl"])

  met <- as.data.table(x[["raw"]][["met"]])
  met <- met[!is.na(PrBaro) & PrBaro > 0, .(DateTime, PrBaro)]
  met[, `:=`(
    hydro_year   = hydro_year(DateTime),
    expected_hpa = standard_pressure(elev),
    lake         = lake_name
  )]
  met
}

#' Run barometric pressure QC across all lakes
#'
#' @param rds_dir character
#' @return data.table with baro data for all lakes
load_all_baro <- function(rds_dir = "data/rds") {
  files <- list.files(rds_dir, pattern = "^full_.*\\.rds$", full.names = FALSE)
  lake_names <- sub("^full_", "", sub("\\.rds$", "", files))

  rbindlist(lapply(lake_names, function(lk) {
    cat("Loading", lk, "...\n")
    tryCatch(load_baro(lk, rds_dir), error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      data.table()
    })
  }), fill = TRUE)
}

#' Boxplot of barometric pressure by hydrological year, per lake
#'
#' @param baro_dt data.table from load_all_baro()
#' @return ggplot object
plot_baro_boxplots <- function(baro_dt) {
  # Summarise expected pressure per lake for reference line
  expected <- unique(baro_dt[, .(lake, expected_hpa)])

  baro_dt[, hy_label := factor(hydro_year)]

  ggplot(baro_dt, aes(x = hy_label, y = PrBaro)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, fill = "steelblue",
                 alpha = 0.4, linewidth = 0.4) +
    geom_hline(data = expected, aes(yintercept = expected_hpa),
               linetype = "dashed", colour = "firebrick", linewidth = 0.6) +
    facet_wrap(~lake, scales = "free_y", ncol = 4) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(
      x = "Hydrological year (Jul-Jun)",
      y = "Barometric pressure (hPa)",
      title = "Barometric pressure QC",
      subtitle = "Annual boxplots vs expected standard atmosphere pressure (red dashed)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8)
    )
}

#' Monthly anomaly bar plot — all lakes stacked vertically
#'
#' Shows monthly mean PrBaro anomaly (observed - expected) as bars,
#' one panel per lake in a single column.
#'
#' @param baro_dt data.table from load_all_baro()
#' @return ggplot object
plot_baro_monthly_anomaly <- function(baro_dt, min_days = 15) {
  # Pre-compute date and year_month columns once
  dt <- baro_dt[, .(PrBaro, expected_hpa, lake,
                     date = as.IDate(DateTime),
                     ym = as.IDate(format(DateTime, "%Y-%m-01")))]

  # Fast aggregation
  monthly <- dt[, .(mean_obs = mean(PrBaro),
                     expected = expected_hpa[1],
                     n_days = uniqueN(date)),
                by = .(lake, ym)]
  monthly <- monthly[n_days >= min_days]
  monthly[, anomaly := mean_obs - expected]

  ggplot(monthly, aes(x = as.Date(ym), y = anomaly, colour = lake)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.6) +
    scale_colour_limno() +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      x = NULL, y = "Pressure anomaly (hPa)",
      title = "Barometric pressure monthly anomaly",
      subtitle = "Monthly mean observed minus expected standard atmosphere pressure",
      colour = "Lake"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
