# ==============================================================================
# scar_input.R
# Shared input handling for the scar-table functions (internal):
#   - checking the columns of a scar table
#   - splitting it into cores and deriving each core's scar vectors
#   - carrying per-core columns (`keep`) into the output
#   - collecting the notes of one call into a single message
#   - setting a seed without disturbing the caller's random stream
# ==============================================================================

# A scar shorter than this has no direction; a start this close to the start
# centroid has no start vector.
.SCAR_LEN_MIN    <- 1e-10
.SCAR_RADIUS_MIN <- 1e-12

.START_COLS <- c("Start_X", "Start_Y", "Start_Z")
.END_COLS   <- c("End_X", "End_Y", "End_Z")

# Check that `data` is a scar table with the columns the caller names.
.check_scar_table <- function(data, id_col = "ID", start_cols = .START_COLS,
                              end_cols = .END_COLS, keep = NULL,
                              need_id = TRUE, min_length = 0) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame of scars, one row per scar.", call. = FALSE)
  }
  if (!is.character(start_cols) || length(start_cols) != 3L ||
      !is.character(end_cols) || length(end_cols) != 3L) {
    stop("`start_cols` and `end_cols` must each name three columns (X, Y, Z).",
         call. = FALSE)
  }
  need <- c(if (need_id) id_col, start_cols, end_cols, keep)
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    stop("`data` lacks column(s): ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  for (cc in c(start_cols, end_cols)) {
    if (!is.numeric(data[[cc]])) {
      stop("Coordinate column `", cc, "` must be numeric.", call. = FALSE)
    }
  }
  if (!is.numeric(min_length) || length(min_length) != 1L ||
      is.na(min_length) || min_length < 0) {
    stop("`min_length` must be a single non-negative number.", call. = FALSE)
  }
  invisible(TRUE)
}

# The usable scars of one set of rows: start points S, unit directions v_hat,
# lengths L, and the start vectors r_hat (each start minus the centroid of the
# usable starts, normalised; NA where undefined). Rows with a missing
# coordinate, zero length, or a length below `min_length` are left out and
# counted.
.scar_vectors <- function(S0, E0, rows = seq_len(nrow(S0)), min_length = 0) {
  storage.mode(S0) <- "double"
  storage.mode(E0) <- "double"
  v <- E0 - S0
  L <- sqrt(rowSums(v^2))
  na    <- !stats::complete.cases(S0, E0)
  zero  <- !na & !(L > .SCAR_LEN_MIN)
  short <- !na & !zero & L < min_length
  keep  <- !na & !zero & !short

  S <- S0[keep, , drop = FALSE]
  v <- v[keep, , drop = FALSE]
  L <- L[keep]
  v_hat <- if (nrow(v)) v / L else v

  cen <- if (nrow(S)) colMeans(S) else rep(NA_real_, 3L)
  r <- if (nrow(S)) sweep(S, 2, cen) else S
  rn <- sqrt(rowSums(r^2))
  good <- rn > .SCAR_RADIUS_MIN
  r_hat <- matrix(NA_real_, nrow(r), 3L)
  if (any(good)) r_hat[good, ] <- r[good, , drop = FALSE] / rn[good]

  list(rows = rows[keep], S = S, v_hat = v_hat, L = L, r_hat = r_hat,
       n = nrow(S), n_missing = sum(na), n_zero = sum(zero),
       n_short = sum(short))
}

# Split a scar table into cores, in order of first appearance, and derive the
# vectors of each. Returns a list with `ids` (the ID values) and `cores` (one
# .scar_vectors() result per core, named by the ID as character).
.scar_cores <- function(data, id_col, start_cols, end_cols, min_length = 0) {
  key <- as.character(data[[id_col]])
  if (anyNA(key)) {
    stop("`", id_col, "` has missing values; every scar needs a core ID.",
         call. = FALSE)
  }
  ids <- unique(data[[id_col]])
  ukey <- as.character(ids)
  S0 <- as.matrix(data[, start_cols, drop = FALSE])
  E0 <- as.matrix(data[, end_cols, drop = FALSE])
  idx <- split(seq_len(nrow(data)), factor(key, levels = ukey))
  cores <- lapply(ukey, function(k) {
    i <- idx[[k]]
    .scar_vectors(S0[i, , drop = FALSE], E0[i, , drop = FALSE], rows = i,
                  min_length = min_length)
  })
  names(cores) <- ukey
  list(ids = ids, cores = cores)
}

# The per-core values of the `keep` columns: one row per core, in the order of
# `ids`. Each column must be constant within a core.
.keep_values <- function(data, id_col, keep, ids) {
  if (is.null(keep) || !length(keep)) return(NULL)
  keep <- setdiff(keep, id_col)
  if (!length(keep)) return(NULL)
  key <- as.character(data[[id_col]])
  ukey <- as.character(ids)
  first <- match(ukey, key)
  out <- data.frame(row.names = seq_along(ukey))
  for (cc in keep) {
    x <- data[[cc]]
    n_distinct <- tapply(x, factor(key, levels = ukey),
                         function(v) length(unique(v[!is.na(v)])))
    bad <- names(n_distinct)[n_distinct > 1L]
    if (length(bad)) {
      stop("`keep` column `", cc, "` varies within core(s): ", .fmt_ids(bad),
           ". Only per-core columns can be kept.", call. = FALSE)
    }
    vals <- x[first]
    has <- tapply(!is.na(x), factor(key, levels = ukey), any)
    miss_first <- is.na(vals) & has
    if (any(miss_first)) {
      vals[miss_first] <- vapply(ukey[miss_first], function(k) {
        v <- x[key == k]
        v[!is.na(v)][1]
      }, x[1])
    }
    out[[cc]] <- vals
  }
  out
}

# The first columns of every per-core output: the ID, then the kept columns.
.core_frame <- function(data, id_col, keep, ids) {
  out <- data.frame(ids, stringsAsFactors = FALSE)
  names(out) <- id_col
  kv <- .keep_values(data, id_col, keep, ids)
  if (!is.null(kv)) out <- cbind(out, kv)
  out
}

# "A, B, C" for up to `max` IDs, then "and k more".
.fmt_ids <- function(ids, max = 8L) {
  ids <- as.character(ids)
  if (length(ids) <= max) return(paste(ids, collapse = ", "))
  paste0(paste(ids[seq_len(max)], collapse = ", "), " and ",
         length(ids) - max, " more")
}

# Emit the notes of one call as a single message.
.say <- function(notes, quiet = FALSE) {
  notes <- notes[nzchar(notes)]
  if (!quiet && length(notes)) message(paste(notes, collapse = "\n"))
  invisible(notes)
}

# Notes on scars dropped from a set of cores.
.dropped_notes <- function(cores) {
  cnt <- function(f) vapply(cores, function(x) x[[f]], 0L)
  out <- character(0)
  lab <- c(n_missing = "with missing coordinates",
           n_zero    = "of zero length",
           n_short   = "shorter than `min_length`")
  for (f in names(lab)) {
    k <- cnt(f)
    if (any(k > 0L)) {
      out <- c(out, sprintf("Left out %d scar(s) %s (%s).", sum(k), lab[[f]],
                            .fmt_ids(names(k)[k > 0L])))
    }
  }
  out
}

# Run `expr` under `seed`, then put the caller's random stream back. A NULL
# seed runs `expr` on the current stream.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

.check_seed <- function(seed) {
  if (!is.null(seed) && (length(seed) != 1L || !is.numeric(seed) || is.na(seed))) {
    stop("`seed` must be a single number or NULL.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}
