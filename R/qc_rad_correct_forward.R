#' Apply radiation correction to new/incoming data
#'
#' Applies a forward correction based on a previous [qc_rad_analyse()]
#' result. With a linear trend model, this simply extrapolates the line
#' forward. This is intended as an interim correction until enough new
#' cloudless days accumulate to re-run the full analysis.
#'
#' @param met Data frame with at least `DateTime` (POSIXct) and `RadSWD`.
#' @param analysis Output from [qc_rad_analyse()].
#' @param damping Damping factor (0-1) controlling how much of the correction
#'   is applied. Default 0.9 for conservative correction.
#'
#' @return `met` with `RadSWD` corrected and a new `RadSWD_raw` column
#'   containing the original values.
#' @export
qc_rad_correct_forward <- function(met, analysis, damping = 0.9) {
  day0 <- min(analysis$ratios$Date)
  local_date <- as.Date(format(met$DateTime, "%Y-%m-%d"))
  date_num <- as.numeric(local_date - day0)

  # Linear trend extrapolates naturally
  predicted_ratio <- predict(analysis$trend,
                             newdata = data.frame(date_num = date_num))
  predicted_ratio <- 1 + damping * (predicted_ratio - 1)
  predicted_ratio <- pmax(predicted_ratio, 0.1)

  met$RadSWD_raw <- met$RadSWD
  met$RadSWD <- met$RadSWD / predicted_ratio
  met
}
