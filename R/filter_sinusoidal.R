#' Sinusoidal envelope filter
#'
#' Flags values that fall above or below a user-defined sinusoidal
#' threshold curve. Useful for seasonal QC bounds — e.g. water
#' temperature should follow a seasonal cycle, and values far outside
#' the expected envelope are suspect.
#'
#' The threshold curve is:
#' \deqn{y(t) = y_{centre} + amplitude \cdot \sin(2\pi \cdot (doy - x_{centre}) / 365.25)}
#'
#' where \code{doy} is fractional day of year and \code{x_centre} is the
#' day of year at which the sinusoid peaks.
#'
#' @param x Numeric vector of values to test.
#' @param datetime POSIXct vector of timestamps.
#' @param y_centre Numeric. Vertical centre of the sinusoid (e.g. mean
#'   annual temperature).
#' @param amplitude Numeric. Half-range of the sinusoid (e.g. seasonal
#'   temperature range / 2).
#' @param x_centre Numeric. Day of year (1–366) at which the sinusoid
#'   peaks. For southern hemisphere water temperature, peak is
#'   ~mid-February (x_centre ≈ 46).
#' @param sides One of \code{"above"} (flag values above the curve),
#'   \code{"below"} (flag values below), or \code{"outside"} (flag
#'   values outside ± offset from the curve).
#' @param offset Numeric. Vertical buffer from the sinusoid before
#'   flagging. Default 0. For \code{sides = "outside"}, this widens
#'   the envelope by \code{offset} on each side.
#' @return Logical vector, same length as \code{x}. \code{TRUE} = flagged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Flag water temperatures above a seasonal upper bound
#' # Peak DOY 46 (~Feb 15), mean 14°C, amplitude 8°C, +3°C buffer
#' flagged <- filter_sinusoidal(
#'   pro$t_wtr_degc, pro$datetime,
#'   y_centre = 14, amplitude = 8, x_centre = 46,
#'   sides = "above", offset = 3
#' )
#'
#' # Flag outside a seasonal envelope (both sides)
#' flagged <- filter_sinusoidal(
#'   pro$t_wtr_degc, pro$datetime,
#'   y_centre = 14, amplitude = 8, x_centre = 46,
#'   sides = "outside", offset = 3
#' )
#' }
filter_sinusoidal <- function(x, datetime,
                              y_centre, amplitude, x_centre,
                              sides = c("above", "below", "outside"),
                              offset = 0) {
  sides <- match.arg(sides)

  # Fractional day of year
  doy <- as.numeric(format(datetime, "%j")) +
         as.numeric(format(datetime, "%H")) / 24 +
         as.numeric(format(datetime, "%M")) / 1440

  # Sinusoidal curve (period = 1 year)
  curve <- y_centre + amplitude * sin(2 * pi * (doy - x_centre) / 365.25)

  valid <- !is.na(x) & !is.na(doy)

  flagged <- switch(sides,
    above   = valid & x > (curve + offset),
    below   = valid & x < (curve - offset),
    outside = valid & (x > (curve + offset) | x < (curve - offset))
  )

  flagged[is.na(flagged)] <- FALSE
  flagged
}


#' Generate sinusoidal curve values for plotting
#'
#' Returns the sinusoidal threshold curve as a data frame, useful for
#' overlaying on plots.
#'
#' @inheritParams filter_sinusoidal
#' @param datetime POSIXct vector — the curve is evaluated at these times.
#' @return Data frame with columns: \code{datetime}, \code{curve},
#'   \code{upper} (curve + offset), \code{lower} (curve - offset).
#'
#' @export
sinusoidal_curve <- function(datetime, y_centre, amplitude, x_centre,
                             offset = 0) {

  doy <- as.numeric(format(datetime, "%j")) +
         as.numeric(format(datetime, "%H")) / 24 +
         as.numeric(format(datetime, "%M")) / 1440

  curve <- y_centre + amplitude * sin(2 * pi * (doy - x_centre) / 365.25)

  data.frame(
    datetime = datetime,
    curve    = curve,
    upper    = curve + offset,
    lower    = curve - offset
  )
}
