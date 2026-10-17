#' Detect outliers using rolling median + scaled MAD
#'
#' Flags values that fall outside a rolling window centred on the local
#' median. The window width adapts to the measurement interval. Uses
#' median and MAD (median absolute deviation) rather than mean and SD
#' for robustness to the very outliers being detected.
#'
#' Two modes:
#' \describe{
#'   \item{\code{sides = "both"}}{Flag values more than \code{n_dev}
#'     MADs from the rolling median in either direction (default).}
#'   \item{\code{sides = "upper"}}{Flag only values above the median +
#'     \code{n_dev} * MAD.}
#'   \item{\code{sides = "lower"}}{Flag only values below the median -
#'     \code{n_dev} * MAD.}
#' }
#'
#' The MAD is scaled by 1.4826 to be consistent with SD for normal data,
#' so \code{n_dev = 3} is roughly equivalent to a 3-sigma test but
#' resistant to contamination.
#'
#' @param x Numeric vector.
#' @param window Integer. Number of observations in the rolling window.
#'   Should be odd (centred); even values are rounded up. Default 101.
#' @param n_dev Numeric. Number of scaled MADs from the median to flag.
#'   Default 3.
#' @param sides One of \code{"both"} (default), \code{"upper"}, or
#'   \code{"lower"}.
#' @param min_mad Numeric. Floor for the scaled MAD, preventing
#'   near-zero MAD (e.g. in constant segments) from flagging everything.
#'   Default 0 (no floor). Set to the sensor resolution for noisy data.
#' @return Logical vector, same length as \code{x}. \code{TRUE} = outlier.
#'
#' @details
#' Uses \code{stats::runmed()} for the rolling median (fast C
#' implementation, Turlach 1995). MAD is computed per window position
#' using a vectorised approach: \code{abs(x - rolling_median)} is
#' computed once, then its own rolling median gives the rolling MAD.
#'
#' Edge handling: the first and last \code{floor(window/2)} values use
#' shrinking windows (runmed's default endpoint rule). These edge
#' estimates are less stable but avoid leaving edges unflagged.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- sin(seq(0, 4*pi, length.out = 500)) + rnorm(500, sd = 0.1)
#' x[c(100, 250, 400)] <- c(10, -8, 12)  # inject spikes
#' outliers <- detect_outliers(x, window = 51, n_dev = 3)
#' sum(outliers)  # should catch the 3 spikes
#'
#' # One-sided: flag only high values (e.g. turbidity spikes)
#' detect_outliers(x, window = 51, n_dev = 3, sides = "upper")
detect_outliers <- function(x, window = 101L, n_dev = 3,
                            sides = c("both", "upper", "lower"),
                            min_mad = 0) {
  sides <- match.arg(sides)
  n <- length(x)
  result <- rep(FALSE, n)
  if (n < 3L) return(result)

  # Ensure odd window, clamp to data length
  window <- as.integer(window)
  if (window %% 2L == 0L) window <- window + 1L
  window <- min(window, n)
  if (window < 3L) window <- 3L

  # Rolling median (fast C implementation)
  roll_med <- stats::runmed(x, k = window, endrule = "median")

  # Rolling MAD: median of |x - rolling_median|
  abs_dev <- abs(x - roll_med)
  roll_mad <- stats::runmed(abs_dev, k = window, endrule = "median") * 1.4826

  # Apply floor
  roll_mad <- pmax(roll_mad, min_mad)

  # Flag outliers
  dev <- x - roll_med
  result <- switch(sides,
    both  = !is.na(dev) & abs(dev) > n_dev * roll_mad,
    upper = !is.na(dev) & dev > n_dev * roll_mad,
    lower = !is.na(dev) & dev < -n_dev * roll_mad
  )

  # NA values are not outliers
  result[is.na(x)] <- FALSE
  result
}
