#' Export dataset with QC flags joined from metadata tables
#'
#' Materialises a QC'd dataset by joining corrected data with flag sources:
#' deployment states, data excisions, interactive edits, and computed range
#' filters. Raw values are never modified — flags are added as
#' \code{{var}_qc} columns alongside the data.
#'
#' Flag sources are evaluated in order and combined with worst-case-wins:
#' \enumerate{
#'   \item \strong{Deployment states} — site-wide spans from \code{deployment.csv}
#'   \item \strong{Data excisions} — var + date range from \code{data_excisions.csv}
#'   \item \strong{Interactive edits} — granular conditions from \code{qc_edits.csv}
#'   \item \strong{Range filters} — computed on the fly from \code{qc_filters.csv}
#' }
#'
#' @param df Data frame of corrected data (profile or met).
#' @param site Site code.
#' @param col_datetime Datetime column (default \code{"datetime"}).
#' @param tier Quality tier for filtering: \code{"all"} (default, no filtering),
#'   \code{"operational"} (flag <= 3), or \code{"research"} (flag == 1 only).
#'   Filtering sets rejected values to NA but preserves the flag column.
#' @param deployment,excisions,edits,filters Optional data frames. If
#'   \code{NULL}, reads from \code{ltvocab::vocab_metadata()}.
#' @param tz Timezone (default \code{"Etc/GMT-12"}).
#' @return \code{df} with \code{{var}_qc} columns added (integer QARTOD flags).
#'   If \code{tier} is not \code{"all"}, values failing the tier threshold
#'   are set to NA.
#'
#' @export
export_qc <- function(df, site, col_datetime = "datetime",
                      tier = c("all", "operational", "research"),
                      deployment = NULL, excisions = NULL,
                      edits = NULL, filters = NULL,
                      tz = "Etc/GMT-12") {
  tier <- match.arg(tier)

  # -- Identify data columns (numeric, not datetime/structural/existing qc) --
  skip_patterns <- paste0("^", col_datetime, "$|_qc$|^profile_id$|^z_target_m$",
                          "|^z_raw_m$|^z_barocomp_m$|^device_")
  all_cols <- names(df)
  data_cols <- all_cols[vapply(df, is.numeric, logical(1L))]
  data_cols <- data_cols[!grepl(skip_patterns, data_cols)]

  # -- Initialise _qc columns to 2 (not evaluated) --------------------------
  for (v in data_cols) {
    qc_col <- paste0(v, "_qc")
    if (!(qc_col %in% all_cols))
      df[[qc_col]] <- 2L
  }

  dt <- df[[col_datetime]]

  # -- 1. Deployment state flags --------------------------------------------
  if (is.null(deployment))
    deployment <- .read_meta_safe("deployment")
  if (!is.null(deployment)) {
    dep <- deployment[deployment$site == site, , drop = FALSE]
    for (i in seq_len(nrow(dep))) {
      from <- .parse_dt(dep$valid_from[i], tz)
      to   <- .parse_dt(dep$valid_to[i], tz)
      if (is.na(from)) from <- min(dt, na.rm = TRUE)
      if (is.na(to))   to   <- max(dt, na.rm = TRUE)

      idx <- which(!is.na(dt) & dt >= from & dt <= to)
      if (length(idx) == 0L) next

      qartod <- .state_to_qartod(dep$state[i],
                                  as.integer(dep$default_qc[i]))
      target <- .resolve_target_vars(dep$vars[i], data_cols)
      for (v in target)
        df <- .flag_worst(df, v, idx, qartod)
    }
  }

  # -- 2. Data excisions ----------------------------------------------------
  if (is.null(excisions))
    excisions <- .read_meta_safe("data_excisions")
  if (!is.null(excisions)) {
    ex <- excisions[excisions$site == site, , drop = FALSE]
    for (i in seq_len(nrow(ex))) {
      from <- .parse_dt(ex$date_from[i], tz)
      to   <- .parse_dt(ex$date_to[i], tz)
      if (is.na(from)) from <- min(dt, na.rm = TRUE)
      if (is.na(to))   to   <- max(dt, na.rm = TRUE)

      target <- .resolve_var_pattern(ex$var[i], data_cols)
      if (identical(target, "all")) target <- data_cols
      target <- intersect(target, data_cols)

      idx <- which(!is.na(dt) & dt >= from & dt <= to)
      for (v in target)
        df <- .flag_worst(df, v, idx, 4L)
    }
  }

  # -- 3. Interactive edits (qc_edits.csv) ----------------------------------
  if (is.null(edits))
    edits <- .read_meta_safe("qc_edits")
  if (!is.null(edits) && nrow(edits) > 0L) {
    ed <- edits[edits$site == site, , drop = FALSE]
    for (i in seq_len(nrow(ed))) {
      from <- .parse_dt(ed$date_from[i], tz)
      to   <- .parse_dt(ed$date_to[i], tz)
      if (is.na(from)) from <- min(dt, na.rm = TRUE)
      if (is.na(to))   to   <- max(dt, na.rm = TRUE)

      var  <- ed$var[i]
      if (is.na(var) || !(var %in% data_cols)) next

      lo   <- suppressWarnings(as.numeric(ed$min[i]))
      hi   <- suppressWarnings(as.numeric(ed$max[i]))
      flag <- suppressWarnings(as.integer(ed$flag[i]))
      if (is.na(flag)) flag <- 4L
      orient <- if (!is.null(ed$orientation)) ed$orientation[i] else "within"

      idx <- which(!is.na(dt) & dt >= from & dt <= to)
      if (length(idx) == 0L) next

      vals <- df[[var]][idx]
      if (identical(orient, "outside")) {
        # Flag values OUTSIDE the range
        bad <- !is.na(vals) & ((!is.na(lo) & vals < lo) | (!is.na(hi) & vals > hi))
      } else {
        # Flag values WITHIN the range (default — "flag these values")
        bad <- !is.na(vals)
        if (!is.na(lo)) bad <- bad & vals >= lo
        if (!is.na(hi)) bad <- bad & vals <= hi
      }
      df <- .flag_worst(df, var, idx[bad], flag)
    }
  }

  # -- 4. Range filters (computed on the fly) --------------------------------
  if (is.null(filters))
    filters <- .read_meta_safe("qc_filters")
  if (!is.null(filters)) {
    # Site-specific override then default
    site_f <- filters[filters$site == site, , drop = FALSE]
    default_f <- filters[filters$site == "default", , drop = FALSE]
    f <- rbind(site_f, default_f[!(default_f$var %in% site_f$var), ])

    for (i in seq_len(nrow(f))) {
      lo <- suppressWarnings(as.numeric(f$low[i]))
      hi <- suppressWarnings(as.numeric(f$high[i]))
      if (is.na(lo) && is.na(hi)) next

      matched <- data_cols[startsWith(data_cols, f$var[i])]
      for (v in matched) {
        vals <- df[[v]]
        valid <- which(!is.na(vals))
        in_range <- rep(TRUE, length(vals))
        if (!is.na(lo)) in_range <- in_range & (is.na(vals) | vals >= lo)
        if (!is.na(hi)) in_range <- in_range & (is.na(vals) | vals <= hi)

        pass <- valid[in_range[valid]]
        fail <- valid[!in_range[valid]]
        df <- .flag_worst(df, v, pass, 1L)
        df <- .flag_worst(df, v, fail, 4L)
      }
    }
  }

  # -- 5. Tier filtering (set rejected values to NA) -------------------------
  if (tier != "all") {
    threshold <- switch(tier, research = 1L, operational = 3L)
    for (v in data_cols) {
      qc_col <- paste0(v, "_qc")
      if (!(qc_col %in% names(df))) next
      reject <- !is.na(df[[qc_col]]) & df[[qc_col]] > threshold
      df[[v]][reject] <- NA
    }
  }

  df
}


# -- internal helpers --------------------------------------------------------

#' Read a metadata table, returning NULL on failure
#' @noRd
.read_meta_safe <- function(table_name) {
  tryCatch({
    if (requireNamespace("ltvocab", quietly = TRUE))
      ltvocab::vocab_metadata(table_name)
    else NULL
  }, error = function(e) NULL)
}

#' Map deployment state to QARTOD flag
#' @noRd
.state_to_qartod <- function(state, default_qc = 2L) {
  if (is.na(state)) return(default_qc)
  switch(state,
    deployed       = if (!is.na(default_qc)) default_qc else 2L,
    service        = 4L,
    settling       = 3L,
    removed        = 4L,
    fouled         = 3L,
    failed         = 4L,
    fault          = 4L,
    maintenance    = 4L,
    battery_failed = 4L,
    default_qc
  )
}

#' Resolve vars field to column names
#' @noRd
.resolve_target_vars <- function(vars_field, data_cols) {
  if (is.na(vars_field) || tolower(trimws(vars_field)) == "all")
    return(data_cols)
  tokens <- strsplit(trimws(vars_field), "\\s+")[[1L]]
  out <- character()
  for (tok in tokens) {
    hits <- data_cols[startsWith(data_cols, tok)]
    out <- c(out, if (length(hits)) hits else tok)
  }
  intersect(unique(out), data_cols)
}

#' Update a QC flag column — worst case wins, never improves
#' @noRd
.flag_worst <- function(df, var, idx, flag) {
  if (length(idx) == 0L) return(df)
  qc_col <- paste0(var, "_qc")
  if (!(qc_col %in% names(df))) df[[qc_col]] <- 2L
  current <- df[[qc_col]][idx]
  # 9 always wins; otherwise higher value = worse
  if (flag == 9L) {
    df[[qc_col]][idx] <- 9L
  } else {
    worse <- is.na(current) | (current != 9L & current < flag)
    update_idx <- idx[worse]
    if (length(update_idx) > 0L)
      df[[qc_col]][update_idx] <- flag
  }
  df
}
