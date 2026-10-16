#' Calculate clear-sky radiation for a met dataframe
#'
#' Adds a `RadClr` column with theoretical clear-sky shortwave radiation
#' (W/m2) calculated via [MeTo::Rso()].
#'
#' @param met Data frame with at least `DateTime` (POSIXct) and `RadSWD`
#'   columns.
#' @param lat Latitude in decimal degrees (negative for southern hemisphere).
#' @param lng Longitude in decimal degrees (east of Greenwich).
#' @param elev Elevation in metres above sea level.
#' @param Lz Longitude of the centre of the local time zone (degrees west of
#'   Greenwich). Default 180 for NZST (UTC+12). The MeTo convention is positive
#'   west of Greenwich.
#'
#' @return `met` with an additional `RadClr` column (W/m2).
#' @export
qc_rad_clearsky <- function(met, lat, lng, elev, Lz = 180) {
  # Detect modal timestep in seconds
  diffs <- as.numeric(diff(met$DateTime), units = "secs")
  dt_sec <- as.numeric(names(sort(table(diffs), decreasing = TRUE))[1])
  tl <- dt_sec / 3600

  # Build a regular sequence covering the full range, compute Rso on that,
  # then join back — avoids MeTo's requirement for consistent time steps
  dt_range <- range(met$DateTime)
  reg_seq <- seq(dt_range[1], dt_range[2], by = dt_sec)

  rso_mj <- MeTo::Rso(
    x = reg_seq,
    lat.deg = lat,
    long.deg = lng,
    elev = elev,
    tl = tl,
    control = list(Lz = Lz)
  )

  rso_wm2 <- rso_mj * 1e6 / dt_sec

  # Match observed timestamps to the regular grid
  reg_df <- data.frame(DateTime = reg_seq, RadClr = rso_wm2)
  met <- merge(met, reg_df, by = "DateTime", all.x = TRUE)
  met <- met[order(met$DateTime), ]
  rownames(met) <- NULL
  met
}
