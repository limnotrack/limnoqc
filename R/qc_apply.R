#' Initialise or update a QC flag column
#'
#' Creates \code{{var}_qc} if absent (default QARTOD 2 = not evaluated),
#' or worsens existing flags. Flags never improve.
#'
#' @param df Data frame.
#' @param var Column name.
#' @param idx Logical vector (rows to flag).
#' @param flag QARTOD flag to apply (1=pass, 2=not eval, 3=suspect, 4=fail, 9=missing).
#' @return Modified data frame.
#' @keywords internal
.qc_flag <- function(df, var, idx, flag) {
  qc_col <- paste0(var, "_qc")
  if (!(qc_col %in% names(df)))
    df[[qc_col]] <- 2L  # not evaluated
  # Severity order: 1 < 2 < 3 < 4; 9 is missing (separate)
  # Only worsen: set flag if it's more severe than current
  current <- df[[qc_col]][idx]
  # 9 (missing) always overwrites; otherwise higher = worse (except 9)
  if (flag == 9L) {
    df[[qc_col]][idx] <- 9L
  } else {
    worse <- is.na(current) | current < flag | current == 9L
    df[[qc_col]][idx & worse] <- flag
  }
  df
}


#' Apply range filters to a data frame
#'
#' For each variable, values outside \code{[low, high]} are flagged 4 (fail)
#' in the corresponding \code{{var}_qc} column. Values within range are
#' flagged 1 (pass) if not already flagged worse.
#'
#' @param df Data frame with datetime and variable columns.
#' @param site Site code. Used to look up site-specific filters; falls back
#'   to \code{"default"} rows.
#' @param filters Data frame of filter rules with columns: \code{site},
#'   \code{var}, \code{low}, \code{high}. If \code{NULL}, reads from
#'   \code{ltvocab::vocab_metadata("qc_filters")}.
#' @param col_datetime Name of the datetime column (default \code{"datetime"}).
#' @return \code{df} with \code{{var}_qc} columns added or updated.
#'   Values outside range are flagged QARTOD 4; values within flagged 1.
#' @export
apply_range_filters <- function(df, site = NULL, filters = NULL,
                                col_datetime = "datetime") {
  if (is.null(filters)) {
    if (!requireNamespace("ltvocab", quietly = TRUE))
      stop("ltvocab required to read qc_filters")
    filters <- ltvocab::vocab_metadata("qc_filters")
  }

  # Site-specific override, then default
  if (!is.null(site)) {
    site_f <- filters[filters$site == site, , drop = FALSE]
    default_f <- filters[filters$site == "default", , drop = FALSE]
    # Site-specific overrides default for same var
    f <- rbind(site_f, default_f[!(default_f$var %in% site_f$var), ])
  } else {
    f <- filters[filters$site == "default", , drop = FALSE]
  }

  data_cols <- setdiff(names(df), c(col_datetime, grep("_qc$", names(df), value = TRUE)))

  for (i in seq_len(nrow(f))) {
    var_pattern <- f$var[i]
    lo  <- suppressWarnings(as.numeric(f$low[i]))
    hi  <- suppressWarnings(as.numeric(f$high[i]))
    if (is.na(lo) && is.na(hi)) next

    # Match columns by var prefix (e.g. "t_wtr" matches "t_wtr_degc")
    matched <- data_cols[startsWith(data_cols, var_pattern)]
    if (length(matched) == 0) next

    for (v in matched) {
      if (!is.numeric(df[[v]])) next
      valid <- !is.na(df[[v]])

      # Flag pass
      in_range <- valid
      if (!is.na(lo)) in_range <- in_range & df[[v]] >= lo
      if (!is.na(hi)) in_range <- in_range & df[[v]] <= hi

      df <- .qc_flag(df, v, which(valid & in_range), 1L)
      df <- .qc_flag(df, v, which(valid & !in_range), 4L)
    }
  }
  df
}


#' Apply deployment state flags
#'
#' Joins deployment spans to data by datetime and applies the
#' \code{default_qc} QARTOD flag for each span. Data outside any deployment
#' span is flagged 2 (not evaluated).
#'
#' @param df Data frame with a datetime column.
#' @param site Site code.
#' @param deployment Data frame of deployment spans with columns: \code{site},
#'   \code{state}, \code{valid_from}, \code{valid_to}, \code{default_qc},
#'   \code{vars}. If \code{NULL}, reads from
#'   \code{ltvocab::vocab_metadata("deployment")}.
#' @param col_datetime Name of the datetime column (default \code{"datetime"}).
#' @param tz Timezone (default \code{"Etc/GMT-12"}).
#' @return \code{df} with \code{{var}_qc} columns updated per deployment state.
#' @export
apply_deployment_flags <- function(df, site, deployment = NULL,
                                   col_datetime = "datetime",
                                   tz = "Etc/GMT-12") {
  if (is.null(deployment)) {
    if (!requireNamespace("ltvocab", quietly = TRUE))
      stop("ltvocab required to read deployment")
    deployment <- ltvocab::vocab_metadata("deployment")
  }

  dep <- deployment[deployment$site == site, , drop = FALSE]
  if (nrow(dep) == 0) return(df)

  dt <- df[[col_datetime]]
  data_cols <- setdiff(names(df), c(col_datetime, grep("_qc$", names(df), value = TRUE)))
  # Only flag numeric data columns
  data_cols <- data_cols[vapply(df[data_cols], is.numeric, logical(1))]

  for (i in seq_len(nrow(dep))) {
    from <- .parse_dt(dep$valid_from[i], tz)
    to   <- .parse_dt(dep$valid_to[i], tz)
    if (is.na(from)) from <- min(dt, na.rm = TRUE)
    if (is.na(to))   to   <- max(dt, na.rm = TRUE)

    qc_val <- suppressWarnings(as.integer(dep$default_qc[i]))
    if (is.na(qc_val)) qc_val <- 2L

    idx <- which(!is.na(dt) & dt >= from & dt <= to)
    if (length(idx) == 0) next

    # Which vars does this span apply to?
    vars_field <- dep$vars[i]
    if (is.na(vars_field) || tolower(trimws(vars_field)) == "all") {
      target_vars <- data_cols
    } else {
      tokens <- strsplit(trimws(vars_field), "\\s+")[[1]]
      target_vars <- character()
      for (tok in tokens) {
        hits <- data_cols[startsWith(data_cols, tok)]
        if (length(hits) > 0) target_vars <- c(target_vars, hits)
        else target_vars <- c(target_vars, tok)
      }
      target_vars <- intersect(unique(target_vars), data_cols)
    }

    # Map deployment state to QARTOD flag
    qartod <- switch(dep$state[i],
      deployed          = qc_val,
      service           = 4L,
      settling          = 3L,
      removed           = 4L,
      fouled            = 3L,
      failed            = 4L,
      fault             = 4L,
      maintenance       = 4L,
      battery_failed    = 4L,
      qc_val  # fallback
    )

    for (v in target_vars) {
      df <- .qc_flag(df, v, idx, qartod)
    }
  }
  df
}


#' Apply data excisions with QC flags
#'
#' Like [excise_from_metadata()] but adds QARTOD 9 (missing) flags
#' rather than (or in addition to) setting values to NA.
#'
#' @param df Data frame.
#' @param site Site code.
#' @param excisions Data frame with columns: \code{site}, \code{var},
#'   \code{date_from}, \code{date_to}, \code{reason}. If \code{NULL}, reads
#'   from \code{ltvocab::vocab_metadata("data_excisions")}.
#' @param col_datetime Name of datetime column (default \code{"datetime"}).
#' @param tz Timezone (default \code{"Etc/GMT-12"}).
#' @param set_na If \code{TRUE} (default), also set excised values to NA.
#' @return \code{df} with excised values flagged 9 (and optionally set to NA).
#' @export
apply_excisions <- function(df, site, excisions = NULL,
                            col_datetime = "datetime",
                            tz = "Etc/GMT-12", set_na = TRUE) {
  if (is.null(excisions)) {
    if (!requireNamespace("ltvocab", quietly = TRUE))
      stop("ltvocab required to read data_excisions")
    excisions <- ltvocab::vocab_metadata("data_excisions")
  }

  ex <- excisions[excisions$site == site, , drop = FALSE]
  if (nrow(ex) == 0) return(df)

  dt <- df[[col_datetime]]
  all_cols <- setdiff(names(df), c(col_datetime, grep("_qc$", names(df), value = TRUE)))

  for (i in seq_len(nrow(ex))) {
    from <- .parse_dt(ex$date_from[i], tz)
    to   <- .parse_dt(ex$date_to[i], tz)
    if (is.na(from)) from <- min(dt, na.rm = TRUE)
    if (is.na(to))   to   <- max(dt, na.rm = TRUE)

    vars <- .resolve_var_pattern(ex$var[i], all_cols)
    if (identical(vars, "all")) vars <- all_cols
    vars <- intersect(vars, names(df))

    idx <- which(!is.na(dt) & dt >= from & dt <= to)
    if (length(idx) == 0 || length(vars) == 0) next

    for (v in vars) {
      df <- .qc_flag(df, v, idx, 9L)
      if (set_na) df[[v]][idx] <- NA
    }
  }
  df
}


#' Flag equilibration scans at the start of each hold
#'
#' After a profiler descends to a new target depth, sensors need time to
#' equilibrate. Scans within \code{seconds} of the first reading at each
#' hold are flagged QARTOD 3 (suspect).
#'
#' @param df Data frame with datetime, depth, and profile_id columns.
#' @param seconds Equilibration window in seconds.
#' @param col_datetime Datetime column (default \code{"datetime"}).
#' @param col_depth Depth column (default \code{"z_meas_m"}).
#' @param col_target Target depth column (default \code{"z_target_m"}).
#'   If present, holds are grouped by target depth within each profile.
#'   If absent, falls back to depth change detection.
#' @param col_profile Profile ID column (default \code{"profile_id"}).
#' @param vars Character vector of variable columns to flag, or \code{NULL}
#'   to flag all numeric columns.
#' @return \code{df} with \code{{var}_qc} columns updated.
#' @export
trim_equilibration <- function(df, seconds,
                               col_datetime = "datetime",
                               col_depth    = "z_meas_m",
                               col_target   = "z_target_m",
                               col_profile  = "profile_id",
                               vars = NULL) {
  if (is.null(vars)) {
    skip <- c(col_datetime, col_depth, col_target, col_profile,
              grep("_qc$", names(df), value = TRUE))
    vars <- setdiff(names(df), skip)
    vars <- vars[vapply(df[vars], is.numeric, logical(1))]
  }

  dt <- df[[col_datetime]]
  has_target <- col_target %in% names(df) && any(!is.na(df[[col_target]]))
  has_profile <- col_profile %in% names(df) && any(!is.na(df[[col_profile]]))

  if (has_target && has_profile) {
    # Group by profile + target depth
    grp <- paste(df[[col_profile]], df[[col_target]], sep = "_")
  } else if (has_profile) {
    # Group by profile only — flag first N seconds of each profile
    grp <- df[[col_profile]]
  } else {
    return(df)  # can't identify holds without profile info
  }

  # Find first timestamp per group
  ugrp <- unique(grp[!is.na(grp)])
  for (g in ugrp) {
    g_idx <- which(grp == g & !is.na(dt))
    if (length(g_idx) == 0) next
    t_start <- min(dt[g_idx])
    equil_idx <- g_idx[as.numeric(difftime(dt[g_idx], t_start, units = "secs")) < seconds]
    if (length(equil_idx) == 0) next
    for (v in vars) {
      df <- .qc_flag(df, v, equil_idx, 3L)
    }
  }
  df
}
