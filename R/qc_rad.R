#' Solar radiation drift detection and correction workflow
#'
#' Wrapper that runs the full radiation QC pipeline:
#' 1. Calculate clear-sky radiation ([qc_rad_clearsky()])
#' 2. Identify cloudless days ([qc_rad_clear_days()])
#' 3. Optionally plot diagnostics ([qc_rad_plot()])
#' 4. Analyse offset and drift ([qc_rad_analyse()])
#' 5. Optionally correct the data ([qc_rad_correct()])
#'
#' @param met Data frame with at least `DateTime` (POSIXct) and `RadSWD`
#'   columns.
#' @param lat Latitude in decimal degrees (negative for southern hemisphere).
#' @param lng Longitude in decimal degrees.
#' @param elev Elevation in metres above sea level.
#' @param Lz Longitude of the centre of the local time zone (degrees west of
#'   Greenwich). Default 180 for NZST.
#' @param min_rad Minimum radiation for daytime classification. Default 20.
#' @param max_roc_sd Maximum normalised rate-of-change SD for cloudless day
#'   detection. Default 0.02.
#' @param plot Logical; produce diagnostic plots? Default `TRUE`.
#' @param correct Logical; apply correction to `RadSWD`? Default `FALSE`.
#'   Set to `TRUE` to auto-correct based on the analysis.
#' @param correct_method Correction method passed to [qc_rad_correct()].
#'   Default `"auto"`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{met}{The (possibly corrected) met data frame with added columns
#'       `RadClr`, `Date`, `ClearDay`, and (if corrected) `RadSWD_raw`.}
#'     \item{analysis}{Output from [qc_rad_analyse()].}
#'   }
#' @export
qc_rad <- function(met, lat, lng, elev, Lz = 180,
                   min_rad = 20, max_roc_sd = 0.02,
                   plot = TRUE, correct = FALSE,
                   correct_method = "auto") {
  # 1. Clear-sky radiation
  met <- qc_rad_clearsky(met, lat = lat, lng = lng, elev = elev, Lz = Lz)

  # 2. Identify cloudless days
  met <- qc_rad_clear_days(met, min_rad = min_rad, max_roc_sd = max_roc_sd)

  n_clear <- length(unique(met$Date[met$ClearDay]))
  message(sprintf("Identified %d cloudless days out of %d total days.",
                  n_clear, length(unique(met$Date))))

  # 3. Plots
  if (plot) {
    qc_rad_plot(met, type = "overlay")
    qc_rad_plot(met, type = "daily")
    qc_rad_plot(met, type = "scatter")
  }

  # 4. Analysis
  analysis <- qc_rad_analyse(met, plot = plot)

  # 5. Correction
  if (correct) {
    met <- qc_rad_correct(met, analysis, method = correct_method)
    message(sprintf("Applied %s correction (mean ratio: %.3f, drift: %.4f/yr).",
                    if (analysis$stable) "offset" else "drift",
                    analysis$mean_ratio, analysis$drift_per_year))
  }

  list(met = met, analysis = analysis)
}
