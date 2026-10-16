#' Set columns to NA within a date range
#'
#' @param df Data frame with a datetime column.
#' @param vars Character vector of column names, or \code{"all"} to target
#'   every column except the datetime column.
#' @param date_from,date_to Start and end of the excision window. Character
#'   strings (parsed as NZST) or POSIXct. \code{NA} means the start/end of
#'   the data.
#' @param col_datetime Name of the datetime column (default \code{"DateTime"}).
#' @param tz Timezone for parsing character dates (default
#'   \code{"Etc/GMT-12"}).
#' @return \code{df} with matched cells set to \code{NA}.
#' @export
excise <- function(df, vars = "all", date_from = NA, date_to = NA,
                   col_datetime = "DateTime", tz = "Etc/GMT-12") {

  dt <- df[[col_datetime]]
  if (is.null(dt)) stop("Column '", col_datetime, "' not found in df")

  from <- .parse_dt(date_from, tz)
  to   <- .parse_dt(date_to,   tz)
  if (is.na(from)) from <- min(dt, na.rm = TRUE)
  if (is.na(to))   to   <- max(dt, na.rm = TRUE)

  if (tolower(vars[1]) == "all") {
    vars <- setdiff(names(df), col_datetime)
  } else {
    vars <- intersect(vars, names(df))
  }

  idx <- !is.na(dt) & dt >= from & dt <= to
  if (any(idx) && length(vars) > 0L) {
    df[idx, vars] <- NA
  }
  df
}


#' Apply all excisions from metadata for a site
#'
#' Reads \code{data_excisions.csv} via \pkg{ltvocab} and applies every
#' matching row to \code{df}. The \code{var} field in the CSV supports:
#' \itemize{
#'   \item \code{"all"} — excise all variable columns
#'   \item A single column name — exact match (e.g. \code{"DOconc"})
#'   \item Space-separated names — multiple exact matches
#'     (e.g. \code{"FlChlr FlPhyc TurbRT"})
#'   \item A variable prefix — matches any column starting with that prefix
#'     (e.g. \code{"t_wtr"} matches \code{t_wtr_degc}, \code{t_wtr_degc_d1},
#'     etc.)
#' }
#'
#' @param df Data frame to excise from.
#' @param site Character scalar: site code to filter excisions.
#' @param col_datetime Name of the datetime column (default \code{"DateTime"}).
#' @param tz Timezone (default \code{"Etc/GMT-12"}).
#' @param excisions Optional data frame of excisions. If \code{NULL} (default),
#'   read from \code{ltvocab::vocab_metadata("data_excisions")}.
#' @return \code{df} with all matching excisions applied.
#' @export
excise_from_metadata <- function(df, site, col_datetime = "DateTime",
                                 tz = "Etc/GMT-12", excisions = NULL) {

  if (is.null(excisions)) {
    if (!requireNamespace("ltvocab", quietly = TRUE))
      stop("ltvocab is required to read excisions metadata")
    excisions <- ltvocab::vocab_metadata("data_excisions")
  }

  ex <- excisions[excisions$site == site, , drop = FALSE]
  if (nrow(ex) == 0L) return(df)

  all_cols <- setdiff(names(df), col_datetime)

  for (i in seq_len(nrow(ex))) {
    vars <- .resolve_var_pattern(ex$var[i], all_cols)
    df <- excise(df, vars = vars,
                 date_from = ex$date_from[i], date_to = ex$date_to[i],
                 col_datetime = col_datetime, tz = tz)
  }
  df
}


# -- helpers ----------------------------------------------------------

#' Parse a scalar date/time value to POSIXct
#' @noRd
.parse_dt <- function(x, tz = "Etc/GMT-12") {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return(as.POSIXct(NA))
  if (inherits(x, "POSIXct")) return(x)
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = tz))
  # try ymd HMS then ymd HM then ymd
  for (fmt in c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d")) {
    out <- as.POSIXct(x, format = fmt, tz = tz)
    if (!is.na(out)) return(out)
  }
  stop("Cannot parse datetime: ", x)
}


#' Resolve a var field from the excisions CSV into column names
#'
#' Supports "all", exact names (space-separated), and prefix matching.
#' @noRd
.resolve_var_pattern <- function(pattern, all_cols) {
  pattern <- trimws(pattern)
  if (tolower(pattern) == "all") return("all")

  # space-separated tokens

  tokens <- strsplit(pattern, "\\s+")[[1L]]

  # try exact match first
  exact <- intersect(tokens, all_cols)

  # for any token that didn't match exactly, try prefix match
  unmatched <- setdiff(tokens, exact)
  prefix_hits <- character(0L)
  for (tok in unmatched) {
    hits <- all_cols[startsWith(all_cols, tok)]
    if (length(hits) > 0L) {
      prefix_hits <- c(prefix_hits, hits)
    } else {
      # no match at all — keep the token so excise() can silently skip it
      exact <- c(exact, tok)
    }
  }

  unique(c(exact, prefix_hits))
}
