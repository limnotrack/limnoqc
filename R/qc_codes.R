# =========================================================================
# NEMS Quality Code System
# =========================================================================
#
# NZ National Environmental Monitoring Standards (NEMS) quality codes
# QC 0–600, with grading functions based on the accuracy chain:
#   model_specs → device_calibrations → sensor_reference
#
# References:
#   NEMS National Quality Code Schema v3.0.1 (April 2026)
#   NEMS Water Temperature v3.0.0 (February 2026)
#   NEMS Dissolved Oxygen v3.0.0 (February 2026)
# =========================================================================


# -- Variables with a NEMS standard ----------------------------------------
.nems_vars <- c("t_wtr", "c_do", "c_dosat", "n_turb")


# =========================================================================
# Registry access
# =========================================================================

#' NEMS quality code registry
#'
#' Returns the 7-row code table (QC 0 through QC 600) with labels,
#' descriptions, colours, and QARTOD mapping.
#'
#' @return A data.frame with columns: code, label, description,
#'   colour_nems, colour_nonnems, qartod_map.
#' @export
qc_nems_codes <- local({
  cache <- NULL
  function(refresh = FALSE) {
    if (!is.null(cache) && !refresh) return(cache)
    if (requireNamespace("ltvocab", quietly = TRUE)) {
      cache <<- tryCatch(
        ltvocab::vocab_read("qc_nems"),
        error = function(e) .qc_nems_fallback()
      )
    } else {
      cache <<- .qc_nems_fallback()
    }
    cache$code <- as.integer(cache$code)
    cache$qartod_map <- as.integer(cache$qartod_map)
    cache
  }
})


#' NEMS quality code colours
#'
#' @param nems Logical. If \code{TRUE} (default), return the official NEMS
#'   colour palette. If \code{FALSE}, return the non-NEMS palette for
#'   internally-graded variables.
#' @return Named character vector: names are QC codes as strings, values
#'   are hex colours.
#' @export
qc_nems_colours <- function(nems = TRUE) {
  codes <- qc_nems_codes()
  col <- if (nems) "colour_nems" else "colour_nonnems"
  stats::setNames(codes[[col]], as.character(codes$code))
}


#' Select QC colour palette for a variable
#'
#' Returns the NEMS palette if the variable has a NEMS standard, otherwise
#' the non-NEMS palette.
#'
#' @param var Character. Variable code (e.g. \code{"t_wtr"}, \code{"c_ph"}).
#' @return Named character vector of QC code colours.
#' @export
qc_colours_for_var <- function(var) {
  is_nems <- var %in% .nems_vars
  qc_nems_colours(nems = is_nems)
}


#' ggplot2 colour scale for QC codes
#'
#' @param nems Logical. Use NEMS palette (\code{TRUE}) or non-NEMS
#'   (\code{FALSE}).
#' @param ... Passed to \code{\link[ggplot2]{scale_colour_manual}}.
#' @return A ggplot2 scale object.
#' @export
scale_colour_qc <- function(nems = TRUE, ...) {
  ggplot2::scale_colour_manual(values = qc_nems_colours(nems), ...)
}


#' ggplot2 fill scale for QC codes
#'
#' @inheritParams scale_colour_qc
#' @return A ggplot2 scale object.
#' @export
scale_fill_qc <- function(nems = TRUE, ...) {
  ggplot2::scale_fill_manual(values = qc_nems_colours(nems), ...)
}


# =========================================================================
# Grading functions
# =========================================================================

#' Grade sensor capability against NEMS thresholds
#'
#' Looks up the manufacturer-stated accuracy for a device model and
#' variable in \code{model_specs.csv}, then compares against the
#' \code{sensor_accuracy} criterion in \code{nems_thresholds.csv}.
#'
#' The device can be identified by either \code{make_model} (the
#' \code{models.csv} key, e.g. \code{"sns_aq_optod_"}) or by
#' \code{manufacturer} + \code{model} (the \code{model_specs.csv} key,
#' e.g. \code{"aqualabo"} + \code{"optod"}).
#'
#' @param var Character. Variable code (e.g. \code{"t_wtr"}).
#' @param make_model Character. Device model ID from \code{models.csv}
#'   (e.g. \code{"sns_aq_optod_"}). If supplied, manufacturer and model
#'   are looked up via \code{models.csv}.
#' @param manufacturer Character. Manufacturer name (e.g. \code{"aqualabo"}).
#'   Used with \code{model} as an alternative to \code{make_model}.
#' @param model Character. Model name (e.g. \code{"optod"}).
#' @param model_specs Data frame. If \code{NULL}, read via
#'   \code{ltvocab::vocab_metadata("model_specs")}.
#' @param models Data frame. If \code{NULL} and \code{make_model} is
#'   supplied, read via \code{ltvocab::vocab_metadata("models")}.
#' @param thresholds Data frame. If \code{NULL}, read via
#'   \code{ltvocab::vocab_read("nems_thresholds")}.
#' @return Integer QC code (600, 500, 400, or 200).
#' @export
qc_grade_sensor <- function(var, make_model = NULL,
                            manufacturer = NULL, model = NULL,
                            model_specs = NULL, models = NULL,
                            thresholds = NULL) {

  if (is.null(model_specs)) model_specs <- .read_model_specs()
  if (is.null(thresholds))  thresholds  <- .read_thresholds()

  # Resolve manufacturer + model from make_model if needed
  if (!is.null(make_model) && (is.null(manufacturer) || is.null(model))) {
    if (is.null(models)) models <- .read_models()
    m <- models[models$make_model == make_model, , drop = FALSE]
    if (nrow(m) == 0L) return(200L)
    manufacturer <- m$manufacturer[1L]
    model        <- m$model[1L]
  }

  if (is.null(manufacturer) || is.null(model)) return(200L)

  # Look up device accuracy in model_specs
  spec <- model_specs[model_specs$manufacturer == manufacturer &
                      model_specs$model == model &
                      model_specs$var == var, , drop = FALSE]
  if (nrow(spec) == 0L || is.na(spec$accuracy[1L])) return(200L)

  accuracy <- as.numeric(spec$accuracy[1L])

  # Look up threshold for this variable + sensor_accuracy criterion
  thresh <- thresholds[thresholds$variable == var &
                       thresholds$criterion == "sensor_accuracy", , drop = FALSE]
  if (nrow(thresh) == 0L) return(200L)

  .grade_value(accuracy, thresh)
}


#' Grade calibration currency
#'
#' Checks time since last calibration for a device. Returns QC code based
#' on whether calibration frequency meets NEMS requirements.
#'
#' @param var Character. Variable code.
#' @param device_id Character. Full device ID.
#' @param datetime POSIXct. Timestamp to grade.
#' @param calibrations Data frame. If \code{NULL}, read via
#'   \code{ltvocab::vocab_metadata("device_calibrations")}.
#' @param thresholds Data frame. If \code{NULL}, read via
#'   \code{ltvocab::vocab_read("nems_thresholds")}.
#' @return Integer QC code.
#' @export
qc_grade_calibration <- function(var, device_id, datetime,
                                 calibrations = NULL, thresholds = NULL) {

  if (is.null(calibrations)) calibrations <- .read_calibrations()
  if (is.null(thresholds))   thresholds   <- .read_thresholds()

  # Find most recent calibration before datetime
  cals <- calibrations[calibrations$device_id == device_id, , drop = FALSE]
  if (nrow(cals) == 0L) return(200L)

  cal_dates <- as.POSIXct(cals$cal_date, tz = "Etc/GMT-12")
  prior <- cal_dates[cal_dates <= datetime]
  if (length(prior) == 0L) return(200L)

  days_since <- as.numeric(difftime(datetime, max(prior), units = "days"))

  # Look up calibration interval threshold
  # NEMS requires at least every 12 months (365 days) for QC 600
  # Use verification_interval_days as proxy (cal interval is typically longer)
  thresh <- thresholds[thresholds$variable == var &
                       thresholds$criterion == "verification_interval_days",
                       , drop = FALSE]
  if (nrow(thresh) == 0L) {
    # Default: 365 days for QC 600
    if (days_since <= 365) return(600L)
    if (days_since <= 545) return(500L)
    return(400L)
  }

  .grade_value(days_since, thresh)
}


#' Grade field verification against NEMS tolerance
#'
#' Computes the difference between a reference measurement and the
#' logged in-situ value, then grades against NEMS verification thresholds.
#'
#' For temperature: simple absolute difference.
#' For DO: compound tolerance = base + (pct_of_ref / 100) * ref_value.
#'
#' @param var Character. Variable code (e.g. \code{"t_wtr"}, \code{"c_dosat"}).
#' @param ref_value Numeric. Reference sensor reading.
#' @param logged_value Numeric. Logged in-situ sensor reading.
#' @param thresholds Data frame. If \code{NULL}, read via
#'   \code{ltvocab::vocab_read("nems_thresholds")}.
#' @return Integer QC code (600, 500, 400, or 200).
#' @export
qc_grade_verification <- function(var, ref_value, logged_value,
                                  thresholds = NULL) {

  if (is.null(thresholds)) thresholds <- .read_thresholds()

  diff_abs <- abs(ref_value - logged_value)

  # Check for compound tolerance (DO has verification_base + verification_pct_of_ref)
  base_row <- thresholds[thresholds$variable == var &
                         thresholds$criterion == "verification_base", , drop = FALSE]
  pct_row  <- thresholds[thresholds$variable == var &
                         thresholds$criterion == "verification_pct_of_ref", , drop = FALSE]

  if (nrow(base_row) > 0L && nrow(pct_row) > 0L) {
    # Compound tolerance: base + pct of reference value
    tol_600 <- as.numeric(base_row$qc600_max) +
               (as.numeric(pct_row$qc600_max) / 100) * abs(ref_value)
    tol_400 <- as.numeric(base_row$qc400_max) +
               (as.numeric(pct_row$qc400_max) / 100) * abs(ref_value)
    tol_500 <- if (!is.na(base_row$qc500_max))
      as.numeric(base_row$qc500_max) +
      (as.numeric(pct_row$qc500_max) / 100) * abs(ref_value)
    else
      (tol_600 + tol_400) / 2  # midpoint if QC 500 not explicit

    if (diff_abs <= tol_600) return(600L)
    if (diff_abs <= tol_500) return(500L)
    if (diff_abs <= tol_400) return(400L)
    return(200L)
  }

  # Simple absolute threshold (e.g. temperature)
  diff_row <- thresholds[thresholds$variable == var &
                         thresholds$criterion == "verification_diff", , drop = FALSE]
  if (nrow(diff_row) == 0L) return(200L)

  .grade_value(diff_abs, diff_row)
}


#' Grade recording interval
#'
#' @param var Character. Variable code.
#' @param interval_s Numeric. Recording interval in seconds.
#' @param thresholds Data frame. If \code{NULL}, read via
#'   \code{ltvocab::vocab_read("nems_thresholds")}.
#' @return Integer QC code.
#' @export
qc_grade_interval <- function(var, interval_s, thresholds = NULL) {

  if (is.null(thresholds)) thresholds <- .read_thresholds()

  thresh <- thresholds[thresholds$variable == var &
                       thresholds$criterion == "recording_interval_s",
                       , drop = FALSE]
  if (nrow(thresh) == 0L) return(600L)

  .grade_value(interval_s, thresh)
}


#' Combine grades into initial quality code
#'
#' Returns the minimum (worst) of all supplied QC codes, per NEMS
#' hierarchy rule: "the initial quality code is the lowest of all boxes
#' ticked."
#'
#' @param ... Integer QC codes. \code{NULL} values are ignored.
#' @return Integer QC code.
#' @export
qc_initial <- function(...) {
  codes <- unlist(list(...))
  codes <- codes[!is.na(codes)]
  if (length(codes) == 0L) return(200L)
  min(codes)
}


# =========================================================================
# Internal helpers
# =========================================================================

#' Grade a value against a threshold row
#'
#' Thresholds are "lower is better" (value <= threshold = pass).
#' @noRd
.grade_value <- function(value, thresh_row) {
  qc600 <- as.numeric(thresh_row$qc600_max[1L])
  qc500 <- as.numeric(thresh_row$qc500_max[1L])
  qc400 <- as.numeric(thresh_row$qc400_max[1L])

  if (!is.na(qc600) && value <= qc600) return(600L)
  if (!is.na(qc500) && value <= qc500) return(500L)
  if (!is.na(qc400) && value <= qc400) return(400L)
  200L
}

#' Read model_specs via ltvocab
#' @noRd
.read_model_specs <- function() {
  if (!requireNamespace("ltvocab", quietly = TRUE))
    stop("ltvocab is required to read model_specs")
  ltvocab::vocab_metadata("model_specs")
}

#' Read models via ltvocab
#' @noRd
.read_models <- function() {
  if (!requireNamespace("ltvocab", quietly = TRUE))
    stop("ltvocab is required to read models")
  ltvocab::vocab_metadata("models")
}

#' Read device_calibrations via ltvocab
#' @noRd
.read_calibrations <- function() {
  if (!requireNamespace("ltvocab", quietly = TRUE))
    stop("ltvocab is required to read device_calibrations")
  ltvocab::vocab_metadata("device_calibrations")
}

#' Read nems_thresholds via ltvocab
#' @noRd
.read_thresholds <- function() {
  if (!requireNamespace("ltvocab", quietly = TRUE))
    stop("ltvocab is required to read nems_thresholds")
  as.data.frame(ltvocab::vocab_read("nems_thresholds"))
}

#' Fallback QC code table when ltvocab is unavailable
#' @noRd
.qc_nems_fallback <- function() {
  data.frame(
    code          = c(0L, 100L, 200L, 300L, 400L, 500L, 600L),
    label         = c("raw", "missing", "unverified", "estimated",
                      "compromised", "fair", "good"),
    description   = c("Unverified raw or auto-processed data",
                      "Data deleted or gap",
                      "Quality not yet determined",
                      "Synthetic or gap-filled data",
                      "Measured but has known issues",
                      "Measured and nearly meets QC 600",
                      "Meets all requirements"),
    colour_nems   = c("#9900FF", "#FF0000", "#8B5A00", "#D3D3D3",
                      "#FFA500", "#00BFFF", "#006400"),
    colour_nonnems = c("#888888", "#CC0066", "#996699", "#CCAA00",
                       "#D45500", "#0099CC", "#006666"),
    qartod_map    = c(2L, 9L, 2L, 3L, 3L, 1L, 1L),
    stringsAsFactors = FALSE
  )
}
