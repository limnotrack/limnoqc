#' Correct radiation data for sensor drift
#'
#' Applies a correction to `RadSWD` based on the analysis from
#' [qc_rad_analyse()]. If the sensor is stable and near 1.0 (within
#' `offset_threshold`), no correction is applied. If stable but offset
#' beyond the threshold, a static (mean ratio) correction is applied.
#' If drifting, a time-varying linear correction is applied.
#'
#' @param met Data frame with `DateTime` and `RadSWD` columns.
#' @param analysis Output from [qc_rad_analyse()].
#' @param method Correction method:
#'   \describe{
#'     \item{`"auto"`}{Choose automatically based on stability and offset.}
#'     \item{`"none"`}{No correction.}
#'     \item{`"static"`}{Constant offset correction using mean ratio.}
#'     \item{`"drift"`}{Time-varying linear drift correction.}
#'   }
#' @param damping Damping factor (0-1) controlling how much of the correction
#'   is applied. Default 0.9 for conservative correction.
#' @param offset_threshold Minimum deviation from 1.0 to trigger static
#'   correction on stable sensors. Default 0.05.
#'
#' @return `met` with `RadSWD` corrected and a new `RadSWD_raw` column
#'   containing the original values.
#' @export
qc_rad_correct <- function(met, analysis,
                           method = c("auto", "none", "static", "drift"),
                           damping = 0.9, offset_threshold = 0.05) {
  method <- match.arg(method)

  if (method == "auto") {
    if (!analysis$stable) {
      method <- "drift"
    } else if (abs(analysis$mean_ratio - 1) > offset_threshold) {
      method <- "static"
    } else {
      method <- "none"
    }
  }

  met$RadSWD_raw <- met$RadSWD

  if (method == "none") {
    # Stable sensor near unity — no correction needed

  } else if (method == "static") {
    # Stable but offset — constant correction
    ratio <- 1 + damping * (analysis$mean_ratio - 1)
    met$RadSWD <- met$RadSWD / ratio

  } else if (method == "drift") {
    day0 <- min(analysis$ratios$Date)
    local_date <- as.Date(format(met$DateTime, "%Y-%m-%d"))
    date_num <- as.numeric(local_date - day0)

    # Linear trend
    predicted_ratio <- predict(analysis$trend,
                               newdata = data.frame(date_num = date_num))
    predicted_ratio <- 1 + damping * (predicted_ratio - 1)
    predicted_ratio <- pmax(predicted_ratio, 0.1)
    met$RadSWD <- met$RadSWD / predicted_ratio
  }

  met
}
