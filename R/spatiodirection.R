# ==============================================================================
# spatiodirection.R
# Spatiodirectional SP-SPHARM: a holistic description of the scar pattern of
# a core from the removal directions, the start positions and their
# correspondence, computed as pair means of Legendre polynomials (base R; no
# Python needed).
# ==============================================================================

.SD_LMAX  <- 6L
.SD_BLOCKS <- list(dir = paste0("dir_", 1:6), start = paste0("start_", 2:6),
                   same = paste0("same_", 1:6), diff = paste0("diff_", 2:6))
.SD_FEATURES <- unlist(.SD_BLOCKS, use.names = FALSE)
.SD_UNIT_TOL <- 1e-6

#' Spatiodirectional SP-SPHARM features of every core
#'
#' Describes the scar pattern of each core as a whole, with no feature
#' specified in advance, by 22 rotation-invariant features computed from the
#' removal directions of its scars, the positions of their starts, and the
#' correspondence between the two (Xiao et al., in prep.; extending
#' SP-SPHARM of Xiao et al., 2027). Use the features to compare the overall
#' scar patterns of groups of cores with [compare_groups()].
#'
#' @inheritParams spi_ei_batch
#' @param features Which blocks of features to return; any of `"dir"`,
#'   `"start"`, `"same"`, `"diff"` (default: all four, 22 features).
#'
#' @return A tibble with one row per core: the ID column, the `keep`
#'   columns, `n_scar` (scars used) and the features. The names of the
#'   feature columns are stored in the attribute `"features"`. Cores with
#'   fewer than three usable scars get `NA`.
#'
#' @details
#' Each scar \eqn{i} is represented by two unit vectors: its removal
#' direction \eqn{\hat v_i}, and its start vector \eqn{\hat u_i}, the
#' direction from the centroid of all the core's starts to its start. For a
#' degree \eqn{l} the features average the Legendre polynomial
#' \eqn{P_l} of the cosine between two vectors:
#'
#' * `dir_l` (l = 1..6), direction-direction pairs: the mean of
#'   \eqn{P_l(\hat v_i \cdot \hat v_j)} over the distinct pairs of scars
#'   \eqn{i \ne j}. These describe the directional relationships among the
#'   scars; `dir_1` is the mean cosine between the directions of two
#'   different scars.
#' * `start_l` (l = 2..6), start-start pairs: the same for the start
#'   vectors, describing how the starts are spread around the core.
#' * `same_l` (l = 1..6), within-scar start-direction pairs: the mean of
#'   \eqn{P_l(\hat v_i \cdot \hat u_i)} over the scars, describing how each
#'   scar's direction relates to its own start position (for example
#'   whether scars run inward from where they start).
#' * `diff_l` (l = 2..6), between-scar direction-start pairs: the mean of
#'   \eqn{P_l(\hat v_i \cdot \hat u_j)} over the ordered pairs \eqn{i \ne j},
#'   describing how a scar's direction relates to the starts of the other
#'   scars.
#'
#' Each feature is multiplied by \eqn{\sqrt{(2l+1)/2}}, which gives every
#' degree the same variance under random directions, so the features can be
#' compared with Euclidean distances without further scaling (do **not**
#' standardise them again, for example in [stats::prcomp()]).
#'
#' **Scar count.** Pairing each scar with itself contributes a term that
#' depends only on the scar count, so those self-pairs are left out: every
#' feature is a mean over single scars or distinct pairs of scars (a
#' U-statistic), and its expectation does not change with the scar count.
#' For the same reason `start_1` and `diff_1`, which contain a term of
#' order \eqn{1/(n-1)}, are not returned. Compare cores with different scar
#' counts directly; do not add the scar count as a covariate.
#'
#' The features use only angles between the vectors of a core, so they do
#' not change when the core is rotated or moved: raw and aligned coordinates
#' give the same features.
#'
#' @references
#' Steinhardt, P. J., Nelson, D. R., & Ronchetti, M. (1983).
#' Bond-orientational order in liquids and glasses. \emph{Physical Review B},
#' \strong{28}(2), 784--805.
#'
#' Hoeffding, W. (1948). A class of statistics with asymptotically normal
#' distribution. \emph{Annals of Mathematical Statistics}, \strong{19}(3),
#' 293--325.
#'
#' @examples
#' set.seed(1)
#' scars <- do.call(rbind, lapply(1:3, function(i) {
#'   s <- matrix(rnorm(30, sd = 10), 10, 3)
#'   e <- s + matrix(rnorm(30, sd = 5), 10, 3)
#'   data.frame(ID = paste0("core", i), Start_X = s[, 1], Start_Y = s[, 2],
#'              Start_Z = s[, 3], End_X = e[, 1], End_Y = e[, 2], End_Z = e[, 3])
#' }))
#' sp <- spharm_spatiodirection(scars)
#' sp
#' attr(sp, "features")
#'
#' # Direction features only
#' spharm_spatiodirection(scars, features = "dir")
#'
#' @seealso [compare_groups()], [spi_ei_batch()], [cpi_pra_batch()]
#' @export
spharm_spatiodirection <- function(data, id_col = "ID",
                                   start_cols = c("Start_X", "Start_Y", "Start_Z"),
                                   end_cols = c("End_X", "End_Y", "End_Z"),
                                   min_length = 0, keep = NULL,
                                   features = c("dir", "start", "same", "diff"),
                                   quiet = FALSE) {
  .check_scar_table(data, id_col, start_cols, end_cols, keep, min_length = min_length)
  .check_flag(quiet, "quiet")
  if (!is.character(features) || !length(features) ||
      !all(features %in% names(.SD_BLOCKS))) {
    stop("`features` must be one or more of \"dir\", \"start\", \"same\", \"diff\".",
         call. = FALSE)
  }
  cols <- unlist(.SD_BLOCKS[names(.SD_BLOCKS) %in% features], use.names = FALSE)
  sc <- .scar_cores(data, id_col, start_cols, end_cols, min_length)
  n <- vapply(sc$cores, function(x) x$n, 0L)
  X <- t(vapply(sc$cores, function(x) {
    if (x$n < 3L) return(stats::setNames(rep(NA_real_, length(.SD_FEATURES)), .SD_FEATURES))
    .sd_features(x$v_hat, x$r_hat)
  }, numeric(length(.SD_FEATURES))))

  out <- .core_frame(data, id_col, keep, sc$ids)
  out$n_scar <- unname(n)
  for (cc in cols) out[[cc]] <- unname(X[, cc])
  out <- tibble::as_tibble(out)
  attr(out, "features") <- cols

  few <- names(n)[n < 3L]
  .say(c(.dropped_notes(sc$cores),
         if (length(few)) sprintf("%d core(s) with fewer than 3 usable scars are NA: %s.",
                                  length(few), .fmt_ids(few))),
       quiet)
  out
}

# The 22 features of one core from its unit directions V and start vectors R
# (rows aligned; R is NA where the start vector is undefined).
.sd_features <- function(V, R) {
  Rc <- R[stats::complete.cases(R), , drop = FALSE]
  cr <- .pair_mean_cross(V, R)
  out <- c(.pair_mean(V), .pair_mean(Rc)[-1], cr$same, cr$diff[-1])
  names(out) <- .SD_FEATURES
  out
}

# For every degree l, sqrt((2l + 1) / 2) times the mean of P_l(u_i . u_j) over
# the distinct pairs of the unit rows of U; NA when n < 3.
.pair_mean <- function(U, lmax = .SD_LMAX) {
  out <- rep(NA_real_, lmax)
  n <- nrow(U)
  if (n < 3L) return(out)
  if (any(abs(sqrt(rowSums(U^2)) - 1) > .SD_UNIT_TOL))
    stop("internal: the rows of U must be unit vectors", call. = FALSE)
  G <- tcrossprod(U)
  G[G > 1] <- 1
  G[G < -1] <- -1
  diag(G) <- 1                 # the self-pairs: P_l(1) = 1, removed below as n
  s <- numeric(lmax)
  p0 <- 1
  p1 <- G
  s[1] <- sum(p1)
  if (lmax > 1L) for (l in 1:(lmax - 1L)) {
    p2 <- ((2 * l + 1) * G * p1 - l * p0) / (l + 1)
    s[l + 1L] <- sum(p2)
    p0 <- p1
    p1 <- p2
  }
  (s - n) / (n * (n - 1)) * sqrt((2 * seq_len(lmax) + 1) / 2)
}

# The cross spectrum of directions V and start vectors U (rows aligned; rows
# with an NA are left out): same_l over single scars, diff_l over the
# n (n - 1) ordered pairs of different scars, both weighted by
# sqrt((2l + 1) / 2); NA when n < 3.
.pair_mean_cross <- function(V, U, lmax = .SD_LMAX) {
  ok <- stats::complete.cases(V) & stats::complete.cases(U)
  V <- V[ok, , drop = FALSE]
  U <- U[ok, , drop = FALSE]
  n <- nrow(V)
  same <- rep(NA_real_, lmax)
  diff <- rep(NA_real_, lmax)
  if (n < 3L) return(list(same = same, diff = diff, n = n))
  C <- tcrossprod(V, U)        # C[i, j] = v_i . u_j
  C[C > 1] <- 1
  C[C < -1] <- -1
  p0 <- 1
  p1 <- C
  same[1] <- mean(diag(p1))
  diff[1] <- (sum(p1) - sum(diag(p1))) / (n * (n - 1))
  if (lmax > 1L) for (l in 1:(lmax - 1L)) {
    p2 <- ((2 * l + 1) * C * p1 - l * p0) / (l + 1)
    same[l + 1L] <- mean(diag(p2))
    diff[l + 1L] <- (sum(p2) - sum(diag(p2))) / (n * (n - 1))
    p0 <- p1
    p1 <- p2
  }
  w <- sqrt((2 * seq_len(lmax) + 1) / 2)
  list(same = same * w, diff = diff * w, n = n)
}
