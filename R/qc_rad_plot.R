#' Plot observed vs clear-sky radiation on cloudless days
#'
#' Produces diagnostic plots comparing observed shortwave radiation against
#' the modelled clear-sky envelope on identified cloudless days.
#'
#' @param met Data frame with `DateTime`, `RadSWD`, `RadClr`, `Date`, and
#'   `ClearDay` columns. Run [qc_rad_clear_days()] first.
#' @param type Type of plot to produce:
#'   \describe{
#'     \item{`"overlay"`}{Time series of observed vs clear-sky, highlighting
#'       cloudless days.}
#'     \item{`"daily"`}{Individual day panels for each cloudless day.}
#'     \item{`"scatter"`}{Scatter of daily max observed vs max clear-sky for
#'       cloudless days.}
#'   }
#' @param max_panels Maximum number of daily panels to show (for
#'   `type = "daily"`). Default 16.
#'
#' @return Invisible `NULL`. Called for side effect (plotting).
#' @export
qc_rad_plot <- function(met, type = c("overlay", "daily", "scatter"),
                        max_panels = 16) {
  type <- match.arg(type)

  if (type == "overlay") {
    plot(met$DateTime, met$RadSWD, type = "l", col = "grey70",
         xlab = "Date", ylab = expression(W/m^2),
         main = "Observed vs Clear-Sky Radiation")
    lines(met$DateTime, met$RadClr, col = "orange", lwd = 0.5)
    clear_idx <- met$ClearDay
    points(met$DateTime[clear_idx], met$RadSWD[clear_idx],
           col = "blue", pch = ".", cex = 2)
    legend("topright",
           legend = c("Observed", "Clear-sky", "Cloudless days"),
           col = c("grey70", "orange", "blue"), lty = c(1, 1, NA),
           pch = c(NA, NA, 16), cex = 0.8)

  } else if (type == "daily") {
    clear_dates <- unique(met$Date[met$ClearDay])
    n <- min(length(clear_dates), max_panels)
    if (n == 0) {
      message("No cloudless days identified.")
      return(invisible(NULL))
    }
    # Sample evenly across the time range
    sel <- clear_dates[seq(1, length(clear_dates),
                           length.out = n)]
    dims <- grDevices::n2mfrow(n)
    old_par <- par(mfrow = dims, mar = c(2, 3, 2, 1))
    on.exit(par(old_par))
    for (d in as.character(sel)) {
      idx <- met$Date == as.Date(d)
      day <- met[idx, ]
      hr <- as.numeric(format(day$DateTime, "%H")) +
        as.numeric(format(day$DateTime, "%M")) / 60
      plot(hr, day$RadSWD, type = "l", col = "blue",
           xlab = "", ylab = "", main = d,
           xlim = c(0, 24), ylim = c(0, max(day$RadClr, na.rm = TRUE) * 1.1))
      lines(hr, day$RadClr, col = "orange", lwd = 2)
    }

  } else if (type == "scatter") {
    clear_dates <- unique(met$Date[met$ClearDay])
    if (length(clear_dates) == 0) {
      message("No cloudless days identified.")
      return(invisible(NULL))
    }
    obs_max <- tapply(met$RadSWD[met$ClearDay], met$Date[met$ClearDay], max,
                      na.rm = TRUE)
    clr_max <- tapply(met$RadClr[met$ClearDay], met$Date[met$ClearDay], max,
                      na.rm = TRUE)
    lims <- range(c(obs_max, clr_max), na.rm = TRUE)
    plot(clr_max, obs_max, pch = 16, col = "blue",
         xlab = expression("Clear-sky max " * (W/m^2)),
         ylab = expression("Observed max " * (W/m^2)),
         main = "Cloudless Day Peak Radiation",
         xlim = lims, ylim = lims)
    abline(0, 1, col = "grey50", lty = 2)
    abline(lm(obs_max ~ clr_max), col = "red")
    legend("topleft",
           legend = c("1:1", "Linear fit"),
           col = c("grey50", "red"), lty = c(2, 1), cex = 0.8)
  }

  invisible(NULL)
}
