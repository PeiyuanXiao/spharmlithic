# ==============================================================================
# outline_iges.R
# Reading the peripheral outline of a core from an IGES file: the curve traced
# along the periphery in CAD software (for example with the Curves tool of
# Geomagic Wrap) and exported as rational B-spline segments (IGES entity 126),
# which are strung head-to-tail into one closed loop.
# ==============================================================================

#' Read core outlines from IGES files
#'
#' Reads the peripheral outline of each core from an IGES file and returns it
#' as a closed loop of points, ready for [fit_intersection_plane()],
#' [cpi_pra_batch()] and [plot_core()].
#'
#' @param files Paths to IGES files (`.igs` / `.iges`), one per core.
#' @param ids Core IDs, one per file, matching the ID column of the scar
#'   table. Default: the file names without extension.
#' @param n_per_seg Number of points sampled along each B-spline segment.
#'   Default 300.
#' @param gap_frac Largest gap between segment ends that is bridged with a
#'   straight chord, as a fraction of the outline's bounding-box diagonal.
#'   Default 0.01.
#'
#' @return A named list (names = `ids`) of numeric matrices with columns
#'   `X`, `Y`, `Z`, one point per row, ordered along the loop with the first
#'   point repeated as the last. Attributes record how the loop was
#'   assembled: `n_segments`, `worst_gap` (largest bridged gap), `n_bridged`,
#'   `n_discarded` (separate curves that were not part of the loop) and
#'   `file`. Files that cannot be read are left out with a warning.
#'
#' @details
#' The outline must be exported as rational B-spline curves (IGES entity
#' type 126), the form Geomagic Wrap and most CAD programs use for curves.
#' Segments may come in any order and orientation; they are chained into
#' connected loops. If an export carries stray curves away from the
#' periphery, the longest closed loop is taken as the outline. Gaps between
#' segment ends up to `gap_frac` of the outline's size are closed with a
#' straight chord; larger gaps are an error, because a missing stretch of
#' outline must not be invented. Files with other curve entities (arcs,
#' polylines, offset curves) are rejected: export the outline as B-splines.
#'
#' The declared IGES unit is ignored. The outline must be in the same
#' coordinate frame as the scar coordinates of its core;
#' [cpi_pra_batch()] checks this.
#'
#' @examples
#' \dontrun{
#' files    <- list.files("outlines", "\\.(igs|iges)$", full.names = TRUE)
#' outlines <- read_outline_iges(files)
#' names(outlines)
#'
#' # File names that differ from the core IDs
#' outlines <- read_outline_iges(files, ids = c("EXP01", "EXP02"))
#' }
#'
#' @seealso [fit_intersection_plane()], [cpi_pra_batch()], [plot_core()]
#' @export
read_outline_iges <- function(files, ids = NULL, n_per_seg = 300, gap_frac = 0.01) {
  if (!is.character(files) || !length(files)) {
    stop("`files` must be one or more paths to IGES files.", call. = FALSE)
  }
  if (is.null(ids)) {
    ids <- sub("\\.(igs|iges)$", "", basename(files), ignore.case = TRUE)
  }
  if (length(ids) != length(files)) {
    stop("`ids` must give one ID per file.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop("Duplicated outline IDs: ", .fmt_ids(unique(ids[duplicated(ids)])), ".",
         call. = FALSE)
  }
  out <- list()
  failed <- character(0)
  for (i in seq_along(files)) {
    loop <- tryCatch(.read_outline_iges(files[i], as.integer(n_per_seg), gap_frac),
                     error = function(e) e)
    if (inherits(loop, "error")) {
      failed <- c(failed, sprintf("%s (%s)", basename(files[i]), conditionMessage(loop)))
      next
    }
    out[[as.character(ids[i])]] <- loop
  }
  if (length(failed)) {
    warning(length(failed), " outline file(s) could not be read and are left out:\n  ",
            paste(failed, collapse = "\n  "), call. = FALSE)
  }
  out
}

# ---- reading the file --------------------------------------------------------------

# IGES is a fixed 80-column record format: columns 1-72 carry content, column
# 73 the section letter, columns 74-80 the sequence number.
.iges_sections <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- formatC(lines, width = -80)          # pad short lines
  split(lines, substr(lines, 73, 73))
}

# Split an IGES free-format parameter string. Hollerith strings (`12HSOME_TEXT`)
# may contain the delimiter, so they are consumed by declared length.
.iges_split <- function(text, pd = ",", rd = ";") {
  # no "H", no Hollerith string: a plain split, up to the first record delimiter
  if (!grepl("H", text, fixed = TRUE)) {
    cut <- regexpr(rd, text, fixed = TRUE)
    if (cut > 0L)
      return(strsplit(paste0(substr(text, 1L, cut - 1L), pd), pd, fixed = TRUE)[[1]])
    out <- strsplit(text, pd, fixed = TRUE)[[1]]
    if (!endsWith(text, pd) && length(out) && !nzchar(trimws(out[length(out)])))
      out <- out[-length(out)]                    # a blank remainder is not a field
    return(out)
  }
  ch <- strsplit(text, "", fixed = TRUE)[[1]]
  out <- character(0)
  buf <- ""
  i <- 1L
  n <- length(ch)
  while (i <= n) {
    c1 <- ch[i]
    if (c1 == pd || c1 == rd) {
      out <- c(out, buf)
      buf <- ""
      i <- i + 1L
      if (c1 == rd) break
      next
    }
    if (c1 == "H" && grepl("^[0-9]+$", trimws(buf))) {
      len <- as.integer(trimws(buf))
      out <- c(out, paste(ch[(i + 1L):(i + len)], collapse = ""))
      buf <- ""
      i <- i + 1L + len
      if (i <= n && (ch[i] == pd || ch[i] == rd)) i <- i + 1L
      next
    }
    buf <- paste0(buf, c1)
    i <- i + 1L
  }
  if (nzchar(trimws(buf))) out <- c(out, buf)
  out
}

# IGES writes doubles in Fortran D-exponent form (1.5D-2).
.iges_num <- function(s) as.numeric(sub("[dD]", "E", trimws(s)))

# Curve entities other than 126 (B-spline) and 110 (line) that an outline
# could have been exported as; the reader does not evaluate them.
.IGES_OTHER_CURVES <- c(`100` = "circular arc", `102` = "composite curve",
                        `104` = "conic arc", `106` = "copious data / polyline",
                        `112` = "parametric spline curve", `130` = "offset curve")

# One IGES file: its entities of type 110 (line) and 126 (rational B-spline).
.read_iges <- function(path) {
  sec <- .iges_sections(path)
  if (is.null(sec[["D"]]) || is.null(sec[["P"]])) {
    stop("not an IGES file (no directory or parameter section)", call. = FALSE)
  }

  # directory entries come in pairs of lines; DE numbers are odd
  D <- sec[["D"]]
  de_type <- integer(0); de_num <- integer(0)
  for (k in seq(1L, length(D), by = 2L)) {
    de_type <- c(de_type, as.integer(trimws(substr(D[k], 1, 8))))
    de_num  <- c(de_num, as.integer(trimws(substr(D[k], 74, 80))))
  }
  other <- intersect(names(.IGES_OTHER_CURVES), as.character(de_type))
  if (length(other)) {
    stop("it contains curve entities that are not read (",
         paste(sprintf("%s %s", other, .IGES_OTHER_CURVES[other]), collapse = ", "),
         "); export the outline as B-spline curves (entity 126)", call. = FALSE)
  }

  # parameter lines are keyed by the DE number in columns 65-72
  P <- sec[["P"]]
  p_key <- as.integer(trimws(substr(P, 65, 72)))
  p_txt <- substr(P, 1, 64)

  entities <- list()
  for (j in seq_along(de_num)) {
    if (!de_type[j] %in% c(110L, 126L)) next      # nothing else is read
    par <- .iges_split(paste(p_txt[p_key == de_num[j]], collapse = ""))
    ent <- list(de = de_num[j], type = de_type[j])
    if (de_type[j] == 110L) {
      v <- .iges_num(par[2:7])
      ent$A <- v[1:3]
      ent$B <- v[4:6]
    } else {
      K <- as.integer(.iges_num(par[2]))          # upper index of the sum
      M <- as.integer(.iges_num(par[3]))          # degree
      i <- 8L
      nk <- K + M + 2L
      ent$degree  <- M
      ent$knots   <- .iges_num(par[i:(i + nk - 1L)]);            i <- i + nk
      ent$weights <- .iges_num(par[i:(i + K)]);                  i <- i + K + 1L
      ent$cps     <- matrix(.iges_num(par[i:(i + 3L * (K + 1L) - 1L)]),
                            ncol = 3, byrow = TRUE);             i <- i + 3L * (K + 1L)
      ent$t0 <- .iges_num(par[i])
      ent$t1 <- .iges_num(par[i + 1L])
    }
    entities[[length(entities) + 1L]] <- ent
  }
  entities
}

# ---- curve evaluation --------------------------------------------------------------

# Cox-de Boor basis values N_{k,degree}(t), k = 1..n, by the iterative
# triangular scheme, one column per value of t. `degree` is read per segment:
# exports may mix cubic and quadratic segments in the same file.
.bspline_basis <- function(t, degree, knots, n) {
  m <- length(knots) - 1L
  N <- matrix(0, m, length(t))
  for (j in seq_along(t)) {
    hit <- which(knots[-length(knots)] <= t[j] & t[j] < knots[-1L])
    if (length(hit)) {
      N[hit[1L], j] <- 1
    } else if (t[j] >= knots[length(knots)]) {    # closed right end
      last <- max(which(knots[-length(knots)] < knots[-1L]))
      N[last, j] <- 1
    }
  }
  for (d in seq_len(degree)) {
    Nn <- matrix(0, m, length(t))
    for (k in seq_len(m - d)) {
      a <- 0
      d1 <- knots[k + d] - knots[k]
      if (d1 > 0) a <- a + (t - knots[k]) / d1 * N[k, ]
      d2 <- knots[k + d + 1L] - knots[k + 1L]
      if (d2 > 0) a <- a + (knots[k + d + 1L] - t) / d2 * N[k + 1L, ]
      Nn[k, ] <- a
    }
    N <- Nn
  }
  N[seq_len(n), , drop = FALSE]
}

# One entity sampled as a polyline (lines as their two endpoints).
.iges_polyline <- function(ent, n = 300L) {
  if (ent$type == 110L) return(rbind(ent$A, ent$B))
  ts <- seq(ent$t0, ent$t1, length.out = n)
  ncp <- nrow(ent$cps)
  B <- .bspline_basis(ts, ent$degree, ent$knots, ncp) * ent$weights
  out <- matrix(NA_real_, n, 3)
  for (j in seq_len(n)) {
    b <- B[, j]
    s <- sum(b)
    out[j, ] <- if (s > 0) as.vector(b %*% ent$cps) / s else ent$cps[ncp, ]
  }
  out
}

# ---- assembling the closed loop ----------------------------------------------------

# String open segments head-to-tail into connected components. Segments are
# stored in no particular order and in either orientation, so a component is
# grown from both ends: at each step the nearest free endpoint is attached to
# whichever end of the chain it is closer to. Each component carries the
# attributes `n_segments`, `closing` (distance from its last point back to its
# first) and `bridged` (the non-zero gaps that were closed to build it).
.chain_components <- function(segs, tol = 0) {
  rest <- seq_along(segs)
  out <- list()
  while (length(rest)) {
    cur <- rest[1L]
    rest <- rest[-1L]
    loop <- segs[[cur]]
    n_seg <- 1L
    bridged <- numeric(0)
    while (length(rest)) {
      head <- loop[1L, ]
      tail <- loop[nrow(loop), ]
      cand <- lapply(rest, function(i) {
        s <- segs[[i]]
        c(sqrt(sum((s[1L, ] - tail)^2)),            # append as-is
          sqrt(sum((s[nrow(s), ] - tail)^2)),       # append reversed
          sqrt(sum((s[nrow(s), ] - head)^2)),       # prepend as-is
          sqrt(sum((s[1L, ] - head)^2)))            # prepend reversed
      })
      D <- do.call(rbind, cand)
      k <- which.min(apply(D, 1, min))
      mode <- which.min(D[k, ])
      gap <- D[k, mode]
      if (gap > tol) break
      s <- segs[[rest[k]]]
      if (mode == 2L || mode == 4L) s <- s[nrow(s):1L, , drop = FALSE]
      loop <- if (mode <= 2L) rbind(loop, s[-1L, , drop = FALSE])
              else rbind(s[-nrow(s), , drop = FALSE], loop)
      if (gap > 0) bridged <- c(bridged, gap)
      n_seg <- n_seg + 1L
      rest <- rest[-k]
    }
    attr(loop, "n_segments") <- n_seg
    attr(loop, "closing") <- sqrt(sum((loop[nrow(loop), ] - loop[1L, ])^2))
    attr(loop, "bridged") <- bridged
    out[[length(out) + 1L]] <- loop
  }
  out
}

# The outline of one IGES file as a closed polyline, first point repeated last.
.read_outline_iges <- function(path, n_per_seg, gap_frac) {
  if (!file.exists(path)) stop("file not found", call. = FALSE)
  ents <- .read_iges(path)
  types <- vapply(ents, function(e) e$type, 0L)
  curves <- ents[types == 126L]
  if (!length(curves)) stop("no B-spline curve (entity 126) in the file", call. = FALSE)
  segs <- lapply(curves, .iges_polyline, n = n_per_seg)

  allP <- do.call(rbind, segs)
  tol <- gap_frac * sqrt(sum((apply(allP, 2, max) - apply(allP, 2, min))^2))

  # .chain_components() always takes the globally nearest free endpoint, so an
  # exact join wins over a bridge wherever one exists; `tol` only ever comes
  # into play where nothing fits exactly
  comps <- .chain_components(segs, tol = tol)

  closed <- vapply(comps, function(p) attr(p, "closing") <= tol, TRUE)
  if (!any(closed)) {
    gaps <- vapply(comps, function(p) attr(p, "closing"), 0)
    stop(sprintf("no closed outline: %d component(s), smallest closing gap %.4g, tolerance %.4g",
                 length(comps), min(gaps), tol), call. = FALSE)
  }
  lens <- vapply(comps, .polyline_length, 0)
  lens[!closed] <- -Inf
  k <- which.max(lens)
  loop <- comps[[k]]

  bridged <- attr(loop, "bridged")
  closing <- attr(loop, "closing")
  n_seg <- attr(loop, "n_segments")
  if (closing > 0) {
    bridged <- c(bridged, closing)
    loop <- rbind(loop, loop[1L, , drop = FALSE])  # close it with a chord
  }
  attributes(loop) <- list(dim = dim(loop))
  colnames(loop) <- c("X", "Y", "Z")
  attr(loop, "n_segments") <- n_seg
  attr(loop, "worst_gap") <- if (length(bridged)) max(bridged) else 0
  attr(loop, "n_bridged") <- length(bridged)
  attr(loop, "n_discarded") <- length(comps) - 1L
  attr(loop, "file") <- basename(path)
  loop
}

# Total length of a polyline.
.polyline_length <- function(P) sum(sqrt(rowSums(diff(P)^2)))
