# ==============================================================================
# stats.R
# Directional statistics for 3-D scar orientation data, one core at a time:
#   - compute_spi       : Scar Pattern Index (Clarkson et al. 2006)
#   - compute_spi_angle : SPI converted to an angle (Clarkson et al. 2006)
#   - compute_ei        : Elongation / Isotropy from the orientation tensor
# Each gives the raw value, the rarefied value and the random interval at the
# core's scar count. For many cores at once see spi_ei_batch().
# ==============================================================================

#' Scar Pattern Index (Clarkson et al. 2006)
#'
#' Computes the Scar Pattern Index of one core: the length of the resultant
#' of its scar directions divided by the number of scars, from 0 (directions
#' cancel out) to 1 (all scars run the same way). By default the function
#' returns three readings:
#'
#' * **raw** SPI, computed from all scars together;
#' * **rarefied** SPI, the mean SPI over every combination of `m` scars
#'   (Xiao et al., in prep.). Its expectation does not change with the scar
#'   count, so it is the value to use when comparing cores or groups of
#'   cores;
#' * the **random interval** at the core's scar count: the central 95% of
#'   the SPI of cores whose scar directions are uniform on the sphere, and
#'   whether the raw SPI lies `"within"` or `"outside"` it.
#'
#' @param dx,dy,dz Numeric vectors of equal length: the X, Y and Z components
#'   of the scar direction vectors of **one** core, for example
#'   `End_X - Start_X` etc., or the `d_x`, `d_y`, `d_z` columns returned by
#'   [align_scar_batch()]. The vectors are normalised to unit length inside the
#'   function; rows of zero length or with missing values are left out.
#' @param lengths Optional numeric vector of scar lengths, one per vector. If
#'   given, each unit direction is weighted by its length, reproducing the
#'   length-weighted SPI of Clarkson et al. (2006). The rarefied value and the
#'   random interval are defined for unweighted SPI only and are `NA` when
#'   `lengths` is given.
#' @param raw Logical. Return the raw SPI (`SPI_raw`)? Default `TRUE`.
#' @param rarefied Logical. Return the rarefied SPI (`SPI_rarefied`)? Default
#'   `TRUE`.
#' @param random_interval Logical. Return the random interval at the core's
#'   scar count (`SPI_lower`, `SPI_upper`) and the position of the raw SPI
#'   relative to it (`SPI_vs_random`)? Default `TRUE`.
#' @param m Number of scars per combination for the rarefied value. Default 3.
#' @param max_comb Largest number of combinations averaged. Above it,
#'   `max_comb` distinct combinations are drawn at random (a Monte Carlo
#'   approximation of the mean over all combinations). Default 2000.
#' @param seed Seed for that random draw; the caller's random stream is
#'   restored afterwards. `NULL` draws from the current stream. Default 2226.
#'
#' @return A one-row tibble with `n_scar` (the number of scars used) and the
#'   columns selected by `raw`, `rarefied` and `random_interval`:
#'   `SPI_raw`, `SPI_rarefied`, `SPI_lower`, `SPI_upper` and `SPI_vs_random`
#'   (`"within"` or `"outside"`). All values are `NA` for cores with fewer
#'   than three scars.
#'
#' @details
#' **Which value for which question.** Raw SPI shifts with the scar count:
#' the fewer the scars, the less their directions cancel out, even when they
#' are random. A single core is therefore judged by comparing its raw SPI
#' with the random interval at its own scar count; cores and groups of cores
#' are compared with the rarefied SPI, whose expectation does not depend on
#' the scar count (Xiao et al., in prep.).
#'
#' **Reading the random interval.** SPI above the upper bound indicates scar
#' directions more concentrated on one side than random directions produce
#' at this scar count. SPI below the lower bound is typical of opposed
#' (bidirectional) flaking, whose two sets of removals cancel out. A value
#' `"within"` the interval cannot be told apart from random directions at
#' this scar count. See [random_interval()].
#'
#' **Changes in 0.2.0.** The function returns a one-row tibble instead of a
#' single number, and normalises the input vectors itself: raw displacement
#' vectors (end minus start) now give the unweighted SPI, and the
#' length-weighted SPI is requested with `lengths`.
#'
#' @references
#' Clarkson, C., Vinicius, L., & Lahr, M. M. (2006). Quantifying flake scar
#' patterning on cores using 3D recording techniques. \emph{Journal of
#' Archaeological Science}, \strong{33}(1), 132--142.
#'
#' Bretzke, K., & Conard, N. J. (2012). Evaluating morphological variability
#' in lithic assemblages using 3D models of stone artifacts. \emph{Journal of
#' Archaeological Science}, \strong{39}(12), 3741--3749.
#'
#' Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique
#' and alternative parameters. \emph{Ecology}, \strong{52}(4), 577--586.
#'
#' @examples
#' # Twelve scars running broadly the same way
#' set.seed(1)
#' d <- matrix(rnorm(36, sd = 0.4), 12, 3)
#' d[, 1] <- d[, 1] + 1
#' compute_spi(d[, 1], d[, 2], d[, 3])
#'
#' # Raw SPI only
#' compute_spi(d[, 1], d[, 2], d[, 3], rarefied = FALSE, random_interval = FALSE)
#'
#' # Length-weighted SPI (Clarkson et al. 2006)
#' compute_spi(d[, 1], d[, 2], d[, 3], lengths = runif(12, 5, 30),
#'             rarefied = FALSE, random_interval = FALSE)
#'
#' @seealso [compute_ei()], [spi_ei_batch()], [random_interval()],
#'   [compute_spi_angle()]
#' @export
compute_spi <- function(dx, dy, dz, lengths = NULL, raw = TRUE, rarefied = TRUE,
                        random_interval = TRUE, m = 3, max_comb = 2000,
                        seed = 2226) {
  .check_index_args(raw, rarefied, random_interval, m, max_comb, seed)
  u <- .unit_input(dx, dy, dz, lengths)
  n <- nrow(u$V)
  weighted <- !is.null(u$w)
  enough <- n >= 3L
  value <- if (enough) .spi_raw(u$V, u$w) else NA_real_

  out <- list(n_scar = n)
  if (raw) out$SPI_raw <- value
  if (rarefied) {
    out$SPI_rarefied <- if (enough && !weighted) {
      .rarefied_fabric(u$V, m, max_comb, seed)$values[["SPI"]]
    } else NA_real_
  }
  if (random_interval) {
    b <- if (enough && !weighted) .ri_bounds("SPI", n) else
      list(lower = NA_real_, upper = NA_real_)
    out$SPI_lower <- b$lower
    out$SPI_upper <- b$upper
    out$SPI_vs_random <- .ri_position(value, b$lower, b$upper)
  }
  .say(c(.single_core_notes(n, "SPI", random_interval),
         if (weighted && (rarefied || random_interval))
           paste("Length-weighted SPI: the rarefied value and the random interval",
                 "are defined for unweighted SPI only and are NA.")))
  tibble::as_tibble(out)
}


#' Convert SPI to a scar-pattern angle (Clarkson et al. 2006)
#'
#' Converts the raw Scar Pattern Index of one core into an angle,
#' \eqn{\theta = \arccos(\mathrm{SPI})}, following Clarkson et al. (2006).
#' SPI = 1 maps to 0 degrees (parallel scars) and SPI = 0 to 90 degrees.
#'
#' @inheritParams compute_spi
#' @param unit Either `"degrees"` (default) or `"radians"`.
#'
#' @return A single numeric value in \eqn{[0, 90]} degrees (or
#'   \eqn{[0, \pi/2]} radians).
#'
#' @details
#' Clarkson et al. (2006) read this angle as the mean angle between any two
#' scars. That reading is only approximate: SPI includes each scar paired
#' with itself. Over the distinct pairs of n unit vectors the mean cosine is
#' \eqn{(n \, \mathrm{SPI}^2 - 1) / (n - 1)}, which is not \eqn{\mathrm{SPI}};
#' it is the degree-1 direction feature of [spharm_spatiodirection()]
#' (`dir_1`, up to the factor \eqn{\sqrt{3/2}}). The angle is computed from
#' the raw SPI and, like it, shifts with the scar count.
#'
#' Before applying `acos()`, the SPI value is clamped to \eqn{[-1, 1]} to
#' guard against floating-point overshoot that would otherwise produce `NaN`.
#'
#' @references
#' Clarkson, C., Vinicius, L., & Lahr, M. M. (2006). Quantifying flake scar
#' patterning on cores using 3D recording techniques. \emph{Journal of
#' Archaeological Science}, \strong{33}(1), 132--142.
#'
#' @examples
#' compute_spi_angle(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))   # 0 (parallel)
#' compute_spi_angle(c(1, 0),    c(0, 1),    c(0, 0))       # 45 (orthogonal pair)
#'
#' @seealso [compute_spi()]
#' @export
compute_spi_angle <- function(dx, dy, dz, lengths = NULL,
                              unit = c("degrees", "radians")) {
  unit  <- match.arg(unit)
  u     <- .unit_input(dx, dy, dz, lengths)
  spi   <- .spi_raw(u$V, u$w)
  # Clamp to [-1, 1] to guard against floating-point overshoot before acos.
  angle <- acos(pmin(pmax(spi, -1), 1))
  if (unit == "degrees") angle * 180 / pi else angle
}


#' Elongation (ER) and Isotropy (IR) ratios from the orientation tensor
#'
#' Constructs the 3 x 3 orientation (fabric) tensor of the scar directions of
#' one core, computes its eigenvalues
#' \eqn{\lambda_1 \ge \lambda_2 \ge \lambda_3}, and returns the elongation
#' and isotropy ratios
#' \deqn{ER = 1 - \lambda_2 / \lambda_1, \quad IR = \lambda_3 / \lambda_1.}
#' As for [compute_spi()], each is returned raw, rarefied, and with its
#' random interval at the core's scar count.
#'
#' @inheritParams compute_spi
#' @param raw Logical. Return the raw ratios (`ER_raw`, `IR_raw`) and the
#'   eigenvalues? Default `TRUE`.
#' @param rarefied Logical. Return the rarefied ratios (`ER_rarefied`,
#'   `IR_rarefied`)? Default `TRUE`.
#' @param random_interval Logical. Return the random intervals (`ER_lower`,
#'   `ER_upper`, `IR_lower`, `IR_upper`) and the position of the raw ratios
#'   relative to them (`ER_vs_random`, `IR_vs_random`)? Default `TRUE`.
#'
#' @return A one-row tibble with `n_scar`, then for ER and for IR the columns
#'   selected by `raw`, `rarefied` and `random_interval` (`*_raw`,
#'   `*_rarefied`, `*_lower`, `*_upper`, `*_vs_random`), and, with
#'   `raw = TRUE`, the eigenvalues `lambda1`, `lambda2`, `lambda3` of the
#'   orientation tensor of all scars. All values are `NA` for cores with
#'   fewer than three scars.
#'
#' @details
#' ER measures how far the scars are aligned along one axis, whether they run
#' the same way or opposite ways; IR measures how far the scar directions
#' spread into three dimensions (low IR: directions concentrated in a plane
#' or along an axis). ER above the upper bound and IR below the lower bound
#' of the random interval indicate directions more concentrated than random.
#'
#' Raw ER and IR shift with the scar count; compare cores and groups with the
#' rarefied values (see [compute_spi()]). Negative eigenvalues caused by
#' numerical error are clamped to zero, and both ratios are `NA` when
#' \eqn{\lambda_1 \approx 0}.
#'
#' **Changes in 0.2.0.** The columns `E` and `I` are renamed `ER_raw` and
#' `IR_raw`; the input vectors are normalised inside the function.
#'
#' @references
#' Lin, S. C., Clarkson, C., Julianto, I. M. A., Ferdianto, A., Jatmiko, &
#' Sutikna, T. (2024). A new method for quantifying flake scar organisation
#' on cores using orientation statistics. \emph{Journal of Archaeological Science},
#' \strong{167}, 105998.
#'
#' McPherron, S. P. (2018). Additional statistical and graphical methods for
#' analyzing site formation processes using artifact orientations.
#' \emph{PLOS ONE}, \strong{13}(1), e0190195.
#'
#' @examples
#' set.seed(1)
#' d <- matrix(rnorm(45), 15, 3)
#' compute_ei(d[, 1], d[, 2], d[, 3])
#'
#' @seealso [compute_spi()], [spi_ei_batch()], [random_interval()]
#' @export
compute_ei <- function(dx, dy, dz, raw = TRUE, rarefied = TRUE,
                       random_interval = TRUE, m = 3, max_comb = 2000,
                       seed = 2226) {
  .check_index_args(raw, rarefied, random_interval, m, max_comb, seed)
  u <- .unit_input(dx, dy, dz)
  n <- nrow(u$V)
  enough <- n >= 3L
  ei <- if (enough) .ei_raw(u$V) else
    c(ER = NA_real_, IR = NA_real_, lambda1 = NA_real_, lambda2 = NA_real_,
      lambda3 = NA_real_)
  rar <- if (rarefied && enough) .rarefied_fabric(u$V, m, max_comb, seed)$values else
    c(SPI = NA_real_, ER = NA_real_, IR = NA_real_)

  out <- list(n_scar = n)
  for (p in c("ER", "IR")) {
    if (raw) out[[paste0(p, "_raw")]] <- ei[[p]]
    if (rarefied) out[[paste0(p, "_rarefied")]] <- rar[[p]]
    if (random_interval) {
      b <- if (enough) .ri_bounds(p, n) else list(lower = NA_real_, upper = NA_real_)
      out[[paste0(p, "_lower")]] <- b$lower
      out[[paste0(p, "_upper")]] <- b$upper
      out[[paste0(p, "_vs_random")]] <- .ri_position(ei[[p]], b$lower, b$upper)
    }
  }
  if (raw) {
    out$lambda1 <- ei[["lambda1"]]
    out$lambda2 <- ei[["lambda2"]]
    out$lambda3 <- ei[["lambda3"]]
  }
  .say(.single_core_notes(n, "ER and IR", random_interval))
  tibble::as_tibble(out)
}

# ---- shared by compute_spi() and compute_ei() -------------------------------------

.check_index_args <- function(raw, rarefied, random_interval, m, max_comb, seed) {
  .check_flag(raw, "raw")
  .check_flag(rarefied, "rarefied")
  .check_flag(random_interval, "random_interval")
  if (!raw && !rarefied) {
    stop("At least one of `raw` and `rarefied` must be TRUE.", call. = FALSE)
  }
  .check_rarefaction_args(m, max_comb)
  .check_seed(seed)
  invisible(TRUE)
}

# Unit direction vectors of one core from three coordinate vectors; rows of
# zero length or with missing values are dropped, and the weights with them.
.unit_input <- function(dx, dy, dz, lengths = NULL) {
  if (!is.numeric(dx) || !is.numeric(dy) || !is.numeric(dz) ||
      length(dx) != length(dy) || length(dx) != length(dz)) {
    stop("`dx`, `dy` and `dz` must be numeric vectors of equal length.",
         call. = FALSE)
  }
  if (!is.null(lengths) &&
      (!is.numeric(lengths) || length(lengths) != length(dx))) {
    stop("`lengths` must be a numeric vector with one value per scar.",
         call. = FALSE)
  }
  M <- cbind(dx, dy, dz)
  storage.mode(M) <- "double"
  L <- sqrt(rowSums(M^2))
  keep <- !is.na(L) & L > .SCAR_LEN_MIN
  if (!is.null(lengths)) keep <- keep & !is.na(lengths)
  V <- M[keep, , drop = FALSE] / L[keep]
  dimnames(V) <- NULL
  list(V = V, w = if (is.null(lengths)) NULL else lengths[keep])
}

.single_core_notes <- function(n, what, random_interval) {
  c(if (n < 3L) sprintf("Fewer than 3 usable scars (%d): %s values are NA.", n, what),
    if (random_interval && n > .ri_max_n() && what != "SPI")
      sprintf("The random interval of %s is tabulated up to n = %d: NA for %d scars.",
              what, .ri_max_n(), n))
}
