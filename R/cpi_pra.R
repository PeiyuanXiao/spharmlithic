# ==============================================================================
# cpi_pra.R
# The plane of intersection of a peripherally exploited core, fitted to its
# peripheral outline, and the two angles each scar makes against it:
#   CPI  Centripetality Index: the plan-view angle between the removal
#        direction and the direction from the centre to the scar's start
#        (180 = towards the centre, 90 = along the edge, 0 = away from it)
#   PRA  Plane-Relative Angle: the angle between the removal direction and
#        the plane (0 = along the plane, 90 = perpendicular to it)
# ==============================================================================

# A start cloud sits this close to its own outline when the two share a
# coordinate frame, as a fraction of the outline's bounding-box diagonal.
# Measured over 26 experimental peripheral cores (Xiao et al., in prep.):
# in frame 0.2-0.6 %, out of frame 16.7-56.6 %.
.FRAME_OK_FRAC <- 0.03
# The fraction of starts the frame check keeps: a start whose true
# initiation point is covered by later scars is digitised inward of the
# outline and says nothing about the frame.
.FRAME_TRIM <- 0.85

#' Fit the plane of intersection of a core to its outline
#'
#' Fits a plane to the peripheral outline of a core by total least squares
#' and finds the centre of the region the outline encloses. The plane and
#' the centre are the reference for the Centripetality Index (CPI) and the
#' Plane-Relative Angle (PRA) of [cpi_pra_batch()].
#'
#' @param outline The outline of one core: a numeric matrix or data frame
#'   with three columns (X, Y, Z), one point per row, **ordered along the
#'   loop**, as returned by [read_outline_iges()]. The loop is closed
#'   automatically if its last point differs from its first. A fitted
#'   `intersection_plane` object is returned unchanged.
#' @param K Number of equally spaced points (by arc length) the plane is
#'   fitted to. Default 64.
#'
#' @return An object of class `intersection_plane`: a list with
#'   \describe{
#'     \item{centre}{The centroid of the area enclosed by the outline (the
#'       centre against which CPI is measured).}
#'     \item{normal}{Unit normal of the fitted plane.}
#'     \item{plane_centroid}{The centroid of the `K` resampled points (a point
#'       on the plane).}
#'     \item{basis}{3 x 3 matrix whose rows are two in-plane axes and the
#'       normal (right-handed).}
#'     \item{flatness}{The smallest singular value of the resampled outline
#'       over the middle one: 0 for a perfectly flat outline; the larger it
#'       is, the more the single plane is an approximation.}
#'     \item{rms_out_of_plane}{Root mean square distance of the outline
#'       points from the plane.}
#'     \item{diagonal}{Length of the outline's bounding-box diagonal.}
#'     \item{loop}{The closed outline (first point repeated as the last).}
#'   }
#'
#' @details
#' The plane is fitted to `K` points spaced equally along the outline, so
#' that densely sampled stretches do not dominate the fit. The centre is the
#' centroid of the enclosed **area**, not of the outline points: the loop is
#' projected onto the fitted plane, the polygon centroid is taken there, and
#' the result is lifted back.
#'
#' The points must follow the outline in order. If the projected outline
#' crosses itself, a warning is given: the points may be out of order, and
#' the centre would then be wrong.
#'
#' @examples
#' # An ellipse tilted out of the XY plane
#' t <- seq(0, 2 * pi, length.out = 200)
#' loop <- cbind(30 * cos(t), 20 * sin(t), 3 * cos(t))
#' fit_intersection_plane(loop)
#'
#' @seealso [read_outline_iges()], [cpi_pra_batch()], [plot_core()]
#' @export
fit_intersection_plane <- function(outline, K = 64) {
  if (inherits(outline, "intersection_plane")) return(outline)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) || K < 8 || K != round(K)) {
    stop("`K` must be a whole number of at least 8.", call. = FALSE)
  }
  P <- .as_loop(outline)
  L <- .resample_arclength(P, as.integer(K))
  pl <- .fit_plane(L)
  centre <- .polygon_centroid(P, pl)
  off <- as.vector(sweep(P, 2, pl$centroid) %*% pl$normal)
  diag_len <- sqrt(sum((apply(P, 2, max) - apply(P, 2, min))^2))

  uv <- sweep(L, 2, pl$centroid) %*% t(pl$basis[1:2, , drop = FALSE])
  k <- .count_crossings(uv)
  if (k > 0L) {
    warning("The outline crosses itself ", k, " time(s) when projected onto ",
            "its plane: its points may not be in order along the loop, and ",
            "the centre may be wrong.", call. = FALSE)
  }
  structure(list(centre = centre, normal = pl$normal,
                 plane_centroid = pl$centroid, basis = pl$basis,
                 flatness = pl$sv[3] / pl$sv[2],
                 rms_out_of_plane = sqrt(mean(off^2)),
                 diagonal = diag_len, loop = P),
            class = "intersection_plane")
}

#' @export
print.intersection_plane <- function(x, digits = 4, ...) {
  f <- function(v) paste(formatC(v, digits = digits, format = "g"), collapse = "  ")
  cat("<intersection_plane>\n")
  cat("  centre (area centroid):", f(x$centre), "\n")
  cat("  normal:                ", f(x$normal), "\n")
  cat("  flatness:              ", f(x$flatness), "\n")
  cat("  rms out of plane:      ", f(x$rms_out_of_plane), "\n")
  cat("  outline:               ", nrow(x$loop) - 1L, "points, diagonal",
      f(x$diagonal), "\n")
  invisible(x)
}

#' CPI and PRA for every peripherally exploited core
#'
#' Computes the Centripetality Index (CPI) and the Plane-Relative Angle (PRA)
#' of every scar of each core that has an outline, and summarises them per
#' core by their mean and Gini's mean difference (GMD).
#'
#' @inheritParams spi_ei_batch
#' @param outlines A named list with one element per core, named by core ID:
#'   the outline points (a three-column matrix or data frame ordered along the
#'   loop, as returned by [read_outline_iges()]) or a plane already fitted
#'   with [fit_intersection_plane()]. Cores without an outline are left out.
#'   The outlines must be in the same coordinate frame as the scar columns
#'   named by `start_cols` and `end_cols`.
#'
#' @return A tibble with one row per core that has an outline: the ID
#'   column, the `keep` columns, `n_scar` (scars used), `n_cpi` (scars with a
#'   defined CPI), `CPI_mean`, `CPI_GMD`, `PRA_mean`, `PRA_GMD` (degrees),
#'   and three checks of the outline: `flatness` (see
#'   [fit_intersection_plane()]), `start_dist` (the typical distance of the
#'   scar starts from the outline, as a fraction of the outline's size) and
#'   `frame_ok` (`start_dist` at most 0.03).
#'
#' @details
#' **The angles.** The plane of intersection is fitted to the outline and its
#' centre is the centroid of the enclosed area (see
#' [fit_intersection_plane()]). CPI is the angle, in the plane, between a
#' scar's removal direction and the direction from the centre to its start:
#' 180 degrees for a scar running towards the centre, 90 for one running
#' along the edge, 0 for one running away from the centre. PRA is the angle
#' between the removal direction and the plane: 0 for a scar running along
#' the plane, 90 for one driven into the core perpendicular to it. PRA does
#' not distinguish the two sides of the plane. Neither angle wraps around,
#' so ordinary means apply.
#'
#' **The core summaries.** The mean summarises the overall centripetality
#' or inclination of a core; GMD, the mean absolute difference between the
#' angles of all pairs of its scars, summarises how much they vary. Both are
#' averages over single scars or pairs of scars, so their expectations do
#' not change with the scar count, and they can be compared directly between
#' cores and groups (for example with [compare_indices()]); the median,
#' interquartile range and standard deviation do not have this property
#' (Xiao et al., in prep.).
#'
#' **When to use.** CPI and PRA describe cores flaked from their periphery
#' whose plane of intersection can be identified, such as discoidal and
#' Levallois cores. Supply outlines only for such cores.
#'
#' **The frame check.** Scar starts are recorded on the periphery, so they
#' lie on or near the outline when scars and outline share a coordinate
#' frame. `start_dist` is the median distance of the starts from the outline
#' (over the 85% of starts closest to it) divided by the outline's
#' bounding-box diagonal. Cores whose `start_dist` exceeds 0.03 are flagged
#' with a warning: their scars and outline are probably in different frames
#' or units, and their angles are not meaningful. Use [plot_core()] to look.
#'
#' Rows with missing coordinates or zero length, and scars shorter than
#' `min_length`, are left out. Cores with fewer than three usable scars get
#' `NA`. A scar whose start lies on the centre, or which runs exactly along
#' the normal, has no CPI; `n_cpi` counts the scars that have one.
#'
#' @references
#' Gini, C. (1912). \emph{Variabilità e mutabilità}. Bologna: C. Cuppini.
#'
#' Yitzhaki, S. (2003). Gini's mean difference: a superior measure of
#' variability for non-normal distributions. \emph{Metron}, \strong{61}(2),
#' 285--316.
#'
#' @examples
#' # A disc-shaped core: twelve scars from the edge towards the centre
#' t <- seq(0, 2 * pi, length.out = 13)[-13]
#' start <- cbind(30 * cos(t), 30 * sin(t), 0)
#' end <- start * 0.4
#' end[, 3] <- -5
#' scars <- data.frame(ID = "disc", Start_X = start[, 1], Start_Y = start[, 2],
#'                     Start_Z = start[, 3], End_X = end[, 1], End_Y = end[, 2],
#'                     End_Z = end[, 3])
#' u <- seq(0, 2 * pi, length.out = 200)
#' outlines <- list(disc = cbind(30 * cos(u), 30 * sin(u), 0))
#' cpi_pra_batch(scars, outlines)
#'
#' @seealso [read_outline_iges()], [fit_intersection_plane()],
#'   [cpi_pra_per_scar()], [plot_core()], [spi_ei_batch()]
#' @export
cpi_pra_batch <- function(data, outlines, id_col = "ID",
                          start_cols = c("Start_X", "Start_Y", "Start_Z"),
                          end_cols = c("End_X", "End_Y", "End_Z"), min_length = 0,
                          keep = NULL, quiet = FALSE) {
  .check_scar_table(data, id_col, start_cols, end_cols, keep, min_length = min_length)
  .check_flag(quiet, "quiet")
  outlines <- .check_outlines(outlines)
  sc <- .scar_cores(data, id_col, start_cols, end_cols, min_length)
  ukey <- names(sc$cores)
  has <- ukey %in% names(outlines)
  unused <- setdiff(names(outlines), ukey)
  if (!any(has)) {
    stop("No outline matches a core ID. Outline names: ", .fmt_ids(names(outlines)),
         "; core IDs: ", .fmt_ids(ukey), ".", call. = FALSE)
  }

  res <- lapply(ukey[has], function(k) {
    plane <- .plane_for(outlines[[k]], k)
    x <- sc$cores[[k]]
    if (x$n < 3L) {
      a <- c(CPI_mean = NA_real_, CPI_GMD = NA_real_, PRA_mean = NA_real_,
             PRA_GMD = NA_real_, n_cpi = NA_real_)
    } else {
      ang <- .ring_angles(x$S, x$v_hat, plane$centre, plane$normal)
      a <- c(CPI_mean = mean(ang$cpi, na.rm = TRUE), CPI_GMD = .gini_md(ang$cpi),
             PRA_mean = mean(ang$pra), PRA_GMD = .gini_md(ang$pra),
             n_cpi = sum(!is.na(ang$cpi)))
    }
    sd <- if (x$n >= 1L) .start_dist(x$S, plane) else NA_real_
    c(a, flatness = plane$flatness, start_dist = sd)
  })
  res <- do.call(rbind, res)

  ids <- sc$ids[has]
  out <- .core_frame(data[as.character(data[[id_col]]) %in% ukey[has], , drop = FALSE],
                     id_col, keep, ids)
  n <- vapply(sc$cores[has], function(x) x$n, 0L)
  out$n_scar <- unname(n)
  out$n_cpi <- as.integer(res[, "n_cpi"])
  for (cc in c("CPI_mean", "CPI_GMD", "PRA_mean", "PRA_GMD", "flatness", "start_dist")) {
    out[[cc]] <- unname(res[, cc])
  }
  out$frame_ok <- out$start_dist <= .FRAME_OK_FRAC

  few <- ukey[has][n < 3L]
  nout <- ukey[!has]
  .say(c(.dropped_notes(sc$cores[has]),
         sprintf("CPI and PRA computed for %d core(s) with an outline.", sum(has)),
         if (length(nout)) sprintf("%d core(s) have no outline and are left out: %s.",
                                   length(nout), .fmt_ids(nout)),
         if (length(unused)) sprintf("%d outline(s) match no core ID: %s.",
                                     length(unused), .fmt_ids(unused)),
         if (length(few)) sprintf("%d core(s) with fewer than 3 usable scars are NA: %s.",
                                  length(few), .fmt_ids(few))),
       quiet)
  bad <- as.character(out[[id_col]])[!is.na(out$frame_ok) & !out$frame_ok]
  if (length(bad)) {
    warning(sprintf(paste("%d core(s) fail the frame check: their scar starts lie far",
                          "from the outline (start_dist > %.2f), so scars and outline",
                          "are probably in different coordinate frames or units: %s.",
                          "Check them with plot_core()."),
                    length(bad), .FRAME_OK_FRAC, .fmt_ids(bad)), call. = FALSE)
  }
  tibble::as_tibble(out)
}

#' CPI and PRA of every scar of one core
#'
#' Returns the Centripetality Index (CPI) and the Plane-Relative Angle (PRA)
#' of each scar of one core, for example to look at their distribution.
#' See [cpi_pra_batch()] for the definitions.
#'
#' @param data The scars of **one** core (a data frame, one row per scar). If
#'   it contains an `id_col` column, it must hold a single ID.
#' @param outline The outline of the core (see [fit_intersection_plane()]),
#'   or a plane already fitted with [fit_intersection_plane()].
#' @inheritParams spi_ei_batch
#'
#' @return A tibble with one row per usable scar: `scar` (the `Scar_ID`
#'   column if present, otherwise the row number in `data`), `CPI` and `PRA`
#'   (degrees), and `dist_to_outline` (the distance of the scar's start from
#'   the outline, in the units of the coordinates).
#'
#' @examples
#' t <- seq(0, 2 * pi, length.out = 13)[-13]
#' start <- cbind(30 * cos(t), 30 * sin(t), 0)
#' end <- start * 0.4
#' end[, 3] <- -5
#' one <- data.frame(Start_X = start[, 1], Start_Y = start[, 2], Start_Z = start[, 3],
#'                   End_X = end[, 1], End_Y = end[, 2], End_Z = end[, 3])
#' u <- seq(0, 2 * pi, length.out = 200)
#' cpi_pra_per_scar(one, cbind(30 * cos(u), 30 * sin(u), 0))
#'
#' @seealso [cpi_pra_batch()], [plot_core()]
#' @export
cpi_pra_per_scar <- function(data, outline, id_col = "ID",
                             start_cols = c("Start_X", "Start_Y", "Start_Z"),
                             end_cols = c("End_X", "End_Y", "End_Z"), min_length = 0) {
  .check_scar_table(data, id_col, start_cols, end_cols, need_id = FALSE,
                    min_length = min_length)
  if (id_col %in% names(data) && length(unique(data[[id_col]])) > 1L) {
    stop("`data` holds more than one core (", .fmt_ids(unique(data[[id_col]])),
         "). Subset it to one core, or use cpi_pra_batch().", call. = FALSE)
  }
  plane <- fit_intersection_plane(outline)
  x <- .scar_vectors(as.matrix(data[, start_cols, drop = FALSE]),
                     as.matrix(data[, end_cols, drop = FALSE]),
                     min_length = min_length)
  ang <- .ring_angles(x$S, x$v_hat, plane$centre, plane$normal)
  scar <- if ("Scar_ID" %in% names(data)) data$Scar_ID[x$rows] else x$rows
  tibble::tibble(scar = scar, CPI = ang$cpi, PRA = ang$pra,
                 dist_to_outline = .point_to_polyline(x$S, plane$loop)$dist)
}

# ---- outlines --------------------------------------------------------------------

.check_outlines <- function(outlines) {
  if (is.matrix(outlines) || is.data.frame(outlines) ||
      inherits(outlines, "intersection_plane")) {
    stop("`outlines` must be a list of outlines named by core ID, ",
         "e.g. list(EXP01 = outline). For one core, use cpi_pra_per_scar().",
         call. = FALSE)
  }
  if (!is.list(outlines) || !length(outlines) || is.null(names(outlines)) ||
      any(!nzchar(names(outlines)))) {
    stop("`outlines` must be a non-empty list named by core ID ",
         "(as returned by read_outline_iges()).", call. = FALSE)
  }
  if (anyDuplicated(names(outlines))) {
    stop("Duplicated outline names: ",
         .fmt_ids(unique(names(outlines)[duplicated(names(outlines))])), ".",
         call. = FALSE)
  }
  outlines
}

# The fitted plane of one outline; warnings name the core.
.plane_for <- function(outline, id) {
  withCallingHandlers(
    tryCatch(fit_intersection_plane(outline),
             error = function(e) stop("Outline of core ", id, ": ",
                                      conditionMessage(e), call. = FALSE)),
    warning = function(w) {
      warning("Outline of core ", id, ": ", conditionMessage(w), call. = FALSE)
      invokeRestart("muffleWarning")
    })
}

# A closed loop of points (first point repeated as the last) from a matrix or
# data frame with three numeric columns.
.as_loop <- function(outline) {
  if (is.data.frame(outline)) {
    num <- vapply(outline, is.numeric, TRUE)
    if (sum(num) != 3L) {
      stop("An outline needs exactly three numeric columns (X, Y, Z).", call. = FALSE)
    }
    outline <- as.matrix(outline[, num, drop = FALSE])
  }
  if (!is.matrix(outline) || !is.numeric(outline) || ncol(outline) != 3L) {
    stop("An outline must be a numeric matrix or data frame with three ",
         "columns (X, Y, Z), one point per row.", call. = FALSE)
  }
  P <- unname(outline)
  storage.mode(P) <- "double"
  if (anyNA(P)) stop("The outline has missing coordinates.", call. = FALSE)
  if (nrow(unique(P)) < 4L) {
    stop("An outline needs at least four distinct points.", call. = FALSE)
  }
  n <- nrow(P)
  if (any(P[n, ] != P[1L, ])) {
    step <- sqrt(rowSums(diff(P)^2))
    gap <- sqrt(sum((P[n, ] - P[1L, ])^2))
    dg <- sqrt(sum((apply(P, 2, max) - apply(P, 2, min))^2))
    if (gap > 2 * max(step) && gap > 0.05 * dg) {
      warning(sprintf(paste("The outline is not closed: its ends are %.3g apart",
                            "(%.0f%% of its size); they are joined with a",
                            "straight segment."), gap, 100 * gap / dg),
              call. = FALSE)
    }
    P <- rbind(P, P[1L, , drop = FALSE])
  }
  colnames(P) <- c("X", "Y", "Z")
  P
}

# K landmarks at equal arc length along a closed polyline (first point
# repeated as the last); the duplicated closing point is dropped, so the K
# landmarks are distinct. Coincident consecutive samples are dropped before
# interpolating.
.resample_arclength <- function(P, K) {
  d <- sqrt(rowSums(diff(P)^2))
  keep <- c(TRUE, d > 0)
  P <- P[keep, , drop = FALSE]
  s <- c(0, cumsum(sqrt(rowSums(diff(P)^2))))
  target <- seq(0, s[length(s)], length.out = K + 1L)[seq_len(K)]
  vapply(1:3, function(j) stats::approx(s, P[, j], xout = target)$y, numeric(K))
}

# Total-least-squares plane through a point set, with a right-handed basis
# (rows e1, e2, normal).
.fit_plane <- function(X) {
  c0 <- colMeans(X)
  sv <- svd(sweep(X, 2, c0))
  B <- t(sv$v)
  n <- B[3, ]
  cr <- .cross3(B[1, ], B[2, ])
  if (sum(cr * n) < 0) B[2, ] <- -B[2, ]
  list(centroid = c0, normal = n, sv = sv$d, basis = B)
}

.cross3 <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

# Centroid of the region a closed loop encloses: the loop is projected onto
# its plane, the polygon centroid is taken there, and the result is lifted
# back.
.polygon_centroid <- function(P, pl) {
  uv <- sweep(P, 2, pl$centroid) %*% t(pl$basis[1:2, , drop = FALSE])
  n <- nrow(uv)
  if (all(abs(uv[n, ] - uv[1, ]) < 1e-12)) uv <- uv[-n, , drop = FALSE]
  x <- uv[, 1]; y <- uv[, 2]
  x1 <- c(x[-1], x[1]); y1 <- c(y[-1], y[1])
  cr <- x * y1 - x1 * y
  A <- sum(cr) / 2
  cx <- sum((x + x1) * cr) / (6 * A)
  cy <- sum((y + y1) * cr) / (6 * A)
  as.vector(pl$centroid + cx * pl$basis[1, ] + cy * pl$basis[2, ])
}

# Number of proper crossings between non-adjacent edges of a closed 2-D
# polygon (vertices as rows, not repeated).
.count_crossings <- function(uv) {
  n <- nrow(uv)
  a <- uv
  b <- uv[c(2:n, 1L), , drop = FALSE]
  cross2 <- function(o, p, q) (p[, 1] - o[, 1]) * (q[, 2] - o[, 2]) -
    (p[, 2] - o[, 2]) * (q[, 1] - o[, 1])
  pr <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i <- pr[, 1]; j <- pr[, 2]
  adj <- (j - i == 1L) | (i == 1L & j == n)
  i <- i[!adj]; j <- j[!adj]
  if (!length(i)) return(0L)
  d1 <- cross2(a[j, , drop = FALSE], b[j, , drop = FALSE], a[i, , drop = FALSE])
  d2 <- cross2(a[j, , drop = FALSE], b[j, , drop = FALSE], b[i, , drop = FALSE])
  d3 <- cross2(a[i, , drop = FALSE], b[i, , drop = FALSE], a[j, , drop = FALSE])
  d4 <- cross2(a[i, , drop = FALSE], b[i, , drop = FALSE], b[j, , drop = FALSE])
  sum(d1 * d2 < 0 & d3 * d4 < 0)
}

# ---- the angles -------------------------------------------------------------------

#   r_hat_i = (S_i - centre) / |S_i - centre|
#   CPI_i   = degrees( arccos( unit(v_par_i) . unit(r_par_i) ) )   in-plane
#   PRA_i   = degrees( arcsin( | v_hat_i . n_hat | ) )             out-of-plane
# Scars whose start sits on the centre have no defined radius, and scars driven
# exactly along the normal have no in-plane direction; both come back as NA in
# `cpi`.
.ring_angles <- function(S, v_hat, centre, normal) {
  r <- sweep(S, 2, centre)
  rn <- sqrt(rowSums(r^2))
  ok <- rn > .SCAR_RADIUS_MIN
  r_hat <- matrix(NA_real_, nrow(r), 3)
  if (any(ok)) r_hat[ok, ] <- r[ok, , drop = FALSE] / rn[ok]
  cosa <- rowSums(.project_to_plane(v_hat, normal) *
                  .project_to_plane(r_hat, normal))
  list(cpi = acos(pmin(pmax(cosa, -1), 1)) * 180 / pi,
       pra = asin(pmin(pmax(abs(as.vector(v_hat %*% normal)), 0), 1)) * 180 / pi,
       r_defined = ok)
}

# Project unit row vectors onto the plane with unit normal n and re-normalise.
# Rows whose in-plane part vanishes come back as NA.
.project_to_plane <- function(M, n) {
  P <- M - outer(as.vector(M %*% n), n)
  len <- sqrt(rowSums(P^2))
  out <- matrix(NA_real_, nrow(P), 3)
  ok <- !is.na(len) & len > 1e-9
  out[ok, ] <- P[ok, , drop = FALSE] / len[ok]
  out
}

# Gini mean difference: the mean of |x_i - x_j| over all pairs of the values
# (missing values dropped), NA below two values. Equals mean(stats::dist(x)).
.gini_md <- function(x) {
  x <- sort(x[!is.na(x)])
  n <- length(x)
  if (n < 2L) return(NA_real_)
  2 * sum((2 * seq_len(n) - n - 1) * x) / (n * (n - 1))
}

# ---- the frame check ---------------------------------------------------------------

# Shortest distance from each row of P to a polyline through Q (exact
# point-to-segment projection).
.point_to_polyline <- function(P, Q) {
  A <- Q[-nrow(Q), , drop = FALSE]
  B <- Q[-1, , drop = FALSE]
  AB <- B - A
  L2 <- pmax(rowSums(AB^2), 1e-300)
  dist <- numeric(nrow(P))
  for (i in seq_len(nrow(P))) {
    w <- sweep(-A, 2, -P[i, ])                    # P[i,] - A, row-wise
    tt <- pmin(pmax(rowSums(w * AB) / L2, 0), 1)
    proj <- A + tt * AB
    dist[i] <- min(sqrt(rowSums(sweep(proj, 2, P[i, ])^2)))
  }
  list(dist = dist)
}

# Median distance of the starts from the outline over the closest 85% of
# starts, as a fraction of the outline's bounding-box diagonal.
.start_dist <- function(S, plane) {
  d <- sort(.point_to_polyline(S, plane$loop)$dist)
  k <- max(min(3L, length(d)), round(.FRAME_TRIM * length(d)))
  stats::median(d[seq_len(k)]) / plane$diagonal
}
