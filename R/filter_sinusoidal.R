#' Sinusoidal envelope filter
#'
#' Flags values that fall above or below a user-defined sinusoidal
#' threshold curve. Useful for seasonal QC bounds — e.g. water
#' temperature should follow a seasonal cycle, and values far outside
#' the expected envelope are suspect.
#'
#' The threshold curve is:
#' \deqn{y(t) = y_{centre} + amplitude \cdot \sin(2\pi \cdot (t - x_{centre}) / period)}
#'
#' where \code{t} is time in the same units as \code{period} (default:
#' decimal year).
#'
#' @param x Numeric vector of values to test.
#' @param datetime POSIXct vector of timestamps.
#' @param y_centre Numeric. Vertical centre of the sinusoid (e.g. mean
#'   annual temperature).
#' @param amplitude Numeric. Half-range of the sinusoid (e.g. seasonal
#'   temperature range / 2).
#' @param x_centre Numeric. Phase offset — the decimal year (or day of
#'   year if \code{period = 365}) at which the sinusoid peaks. For
#'   southern hemisphere water temperature, peak is ~mid-February
#'   (x_centre ≈ 2026.12 or DOY ≈ 45).
#' @param period Numeric. Period in the same units as \code{x_centre}.
#'   Default 365.25 (days, matching DOY-based x_centre). Set to 1 for
#'   decimal-year x_centre.
#' @param sides One of \code{"above"} (flag values above the curve),
#'   \code{"below"} (flag values below), or \code{"outside"} (flag
#'   values outside ± the curve, treating amplitude as a band).
#' @param offset Numeric. Additional vertical offset from the sinusoid
#'   before flagging. Default 0. For \code{sides = "outside"}, this
#'   widens the envelope by \code{offset} on each side.
#' @param time_unit How to convert datetime to the sinusoidal x-axis.
#'   \code{"doy"} (default) = day of year (1–366), matching
#'   \code{period = 365.25}. \code{"decimal_year"} = fractional year,
#'   matching \code{period = 1}.
#' @return Logical vector, same length as \code{x}. \code{TRUE} = flagged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Flag water temperatures above a seasonal upper bound
#' # Peak ~Feb 15 (DOY 46), mean 14°C, amplitude 8°C, + 3°C offset
#' flagged <- filter_sinusoidal(
#'   pro$t_wtr_degc, pro$datetime,
#'   y_centre = 14, amplitude = 8, x_centre = 46,
#'   period = 365.25, sides = "above", offset = 3
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
                              period = 365.25,
                              sides = c("above", "below", "outside"),
                              offset = 0,
                              time_unit = c("doy", "decimal_year")) {
  sides     <- match.arg(sides)
  time_unit <- match.arg(time_unit)

  n <- length(x)
  flagged <- rep(FALSE, n)

  # Convert datetime to sinusoidal x-axis
  t_val <- switch(time_unit,
    doy = {
      as.numeric(format(datetime, "%j")) +
        as.numeric(format(datetime, "%H")) / 24 +
        as.numeric(format(datetime, "%M")) / 1440
    },
    decimal_year = {
      yr  <- as.numeric(format(datetime, "%Y"))
      doy <- as.numeric(format(datetime, "%j"))
      yr + (doy - 1) / 365.25
    }
  )

  # Sinusoidal curve
  curve <- y_centre + amplitude * sin(2 * pi * (t_val - x_centre) / period)

  valid <- !is.na(x) & !is.na(t_val)

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
                             period = 365.25, offset = 0,
                             time_unit = c("doy", "decimal_year")) {
  time_unit <- match.arg(time_unit)

  t_val <- switch(time_unit,
    doy = {
      as.numeric(format(datetime, "%j")) +
        as.numeric(format(datetime, "%H")) / 24 +
        as.numeric(format(datetime, "%M")) / 1440
    },
    decimal_year = {
      yr  <- as.numeric(format(datetime, "%Y"))
      doy <- as.numeric(format(datetime, "%j"))
      yr + (doy - 1) / 365.25
    }
  )

  curve <- y_centre + amplitude * sin(2 * pi * (t_val - x_centre) / period)

  data.frame(
    datetime = datetime,
    curve    = curve,
    upper    = curve + offset,
    lower    = curve - offset
  )
}
