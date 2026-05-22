# ==============================================================================
# stats.R
# Directional statistics for 3-D scar orientation data:
#   - compute_SPI       : Scar Pattern Index (Clarkson et al. 2006)
#   - compute_spi_angle : SPI converted to expected pairwise angle
#   - compute_EI        : Elongation / Isotropy from the orientation tensor
# ==============================================================================

#' Scar Pattern Index (Clarkson et al. 2006)
#'
#' Computes the ratio of the resultant vector magnitude to the total scar
#' length, on a scale from 0 (random orientation) to 1 (perfect alignment).
#'
#' **Default behaviour:** unweighted (Bretzke & Conard 2012) — all scars
#' contribute equally regardless of length. Pass `lengths` to get the
#' length-weighted SPI (Clarkson et al. 2006).
#'
#' @param dx,dy,dz Numeric vectors of equal length. The X, Y, Z components
#'   of the scar direction vectors. These are typically unit direction
#'   vectors returned by [align_scar_batch()] (columns `d_x`, `d_y`, `d_z`).
#' @param lengths Optional numeric vector of scar lengths. If provided,
#'   each direction vector is scaled by its length before summation,
#'   reproducing the original length-weighted SPI of Clarkson et al. (2006).
#'   If `NULL` (the default), all scars contribute equally regardless of
#'   length, following the unit-vector variant of Bretzke & Conard (2012).
#'
#' @return A single numeric value in \eqn{[0, 1]}. Values close to 1
#'   indicate strong preferred orientation; values close to 0 indicate
#'   isotropic / random patterning.
#'
#' @details
#' Two variants of SPI are supported:
#' \itemize{
#'   \item **Length-weighted** (Clarkson et al. 2006): pass raw scar
#'     displacement vectors so that longer scars contribute proportionally
#'     more to the resultant, or equivalently pass unit vectors together
#'     with `lengths`. This reproduces the original definition.
#'   \item **Unweighted** (Bretzke & Conard 2012): leave `lengths = NULL`.
#'     Every scar contributes equally regardless of length — useful when
#'     scar lengths are unreliable or when only direction matters.
#' }
#'
#' To reproduce the length-weighted definition from aligned batch output:
#' \preformatted{
#'   lens <- sqrt((aligned$e_x - aligned$s_x)^2 +
#'                (aligned$e_y - aligned$s_y)^2 +
#'                (aligned$e_z - aligned$s_z)^2)
#'   compute_SPI(aligned$d_x, aligned$d_y, aligned$d_z, lengths = lens)
#' }
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
#' @examples
#' # Perfectly aligned vectors along X
#' compute_SPI(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))   # 1
#'
#' # Random-like distribution (close to 0)
#' set.seed(1)
#' compute_SPI(rnorm(50), rnorm(50), rnorm(50))
#'
#' \dontrun{
#' # Unweighted (default) — every scar counts equally
#' compute_SPI(aligned$d_x, aligned$d_y, aligned$d_z)
#'
#' # Length-weighted — longer scars contribute more
#' lens <- get_scar_length(aligned)
#' compute_SPI(aligned$d_x, aligned$d_y, aligned$d_z, lengths = lens)
#' }
#'
#' @seealso [compute_spi_angle()], [compute_EI()], [get_scar_length()]
#' @export
compute_SPI <- function(dx, dy, dz, lengths = NULL) {
  if (!is.null(lengths)) {
    dx <- dx * lengths
    dy <- dy * lengths
    dz <- dz * lengths
  }
  resultant_magnitude <- sqrt(sum(dx)^2 + sum(dy)^2 + sum(dz)^2)
  total_length        <- sum(sqrt(dx^2 + dy^2 + dz^2))
  resultant_magnitude / total_length
}


#' Convert SPI to a scar-pattern angle (Clarkson et al. 2006)
#'
#' Converts the Scar Pattern Index into the expected pairwise angle
#' between two randomly-selected scars, following Clarkson et al.'s (2006)
#' interpretation: \eqn{\theta = \arccos(\mathrm{SPI})}.
#'
#' SPI = 1 maps to 0 degrees (parallel scars); SPI = 0 maps to 90 degrees
#' (uniformly random pairwise angles).
#'
#' @inheritParams compute_SPI
#' @param unit Either `"degrees"` (default) or `"radians"`.
#'
#' @return A single numeric value in \eqn{[0, 90]} degrees (or
#'   \eqn{[0, \pi/2]} radians).
#'
#' @details
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
#' compute_spi_angle(c(1, 0),    c(0, 1),    c(0, 0))       # 90 (orthogonal pair)
#'
#' @seealso [compute_SPI()]
#' @export
compute_spi_angle <- function(dx, dy, dz, lengths = NULL,
                              unit = c("degrees", "radians")) {
  unit  <- match.arg(unit)
  spi   <- compute_SPI(dx, dy, dz, lengths = lengths)
  # Clamp to [-1, 1] to guard against floating-point overshoot before acos.
  angle <- acos(pmin(pmax(spi, -1), 1))
  if (unit == "degrees") angle * 180 / pi else angle
}


#' Elongation (E) and Isotropy (I) from the orientation tensor
#'
#' Constructs the 3×3 orientation (fabric) tensor from unit direction
#' vectors, computes its eigenvalues
#' \eqn{\lambda_1 \ge \lambda_2 \ge \lambda_3}, and returns the shape
#' descriptors:
#' \deqn{E = 1 - \lambda_2 / \lambda_1, \quad I = \lambda_3 / \lambda_1.}
#'
#' @param ux,uy,uz Numeric vectors of equal length. The X, Y, Z components
#'   of **unit** direction vectors. Non-unit vectors are accepted without
#'   error but will produce meaningless results; normalisation is the
#'   caller's responsibility.
#'
#' @return A one-row data frame with columns:
#' \describe{
#'   \item{E}{Elongation index in \eqn{[0, 1]}: 1 = perfectly linear,
#'     0 = planar / isotropic.}
#'   \item{I}{Isotropy index in \eqn{[0, 1]}: 1 = perfectly isotropic,
#'     0 = linear / planar.}
#'   \item{lambda1, lambda2, lambda3}{The three eigenvalues, sorted in
#'     decreasing order.}
#' }
#'
#' @details
#' Negative eigenvalues caused by numerical error are clamped to zero.
#' Both `E` and `I` are returned as `NA` when \eqn{\lambda_1 \approx 0}
#' (degenerate tensor, threshold \eqn{10^{-10}}).
#'
#' The orientation tensor and the E/I descriptors follow Lin et al. (2024),
#' whose implementation adapts the fabric analysis of McPherron (2018).
#'
#' The one-row data frame return type is designed to play well with
#' [dplyr::group_modify()], allowing direct per-specimen computation.
#'
#' @references
#' Lin, S. C., Clarkson, C., Julianto, I. M. A., Ferdianto, A., & Sutikna,
#' T. (2024). A new method for quantifying flake scar organisation on cores
#' using orientation statistics. \emph{Journal of Archaeological Science},
#' \strong{167}, 105998.
#'
#' McPherron, S. P. (2018). Additional statistical and graphical methods for
#' analyzing site formation processes using artifact orientations.
#' \emph{PLOS ONE}, \strong{13}(1), e0190195.
#'
#' @examples
#' # Strongly elongated — vectors mostly aligned along X
#' compute_EI(c(1, 1, 0.9), c(0, 0, 0.1), c(0, 0, 0))
#'
#' # Isotropic — vectors evenly distributed along all three axes
#' ux <- c(1, -1,  0,  0,  0,  0)
#' uy <- c(0,  0,  1, -1,  0,  0)
#' uz <- c(0,  0,  0,  0,  1, -1)
#' compute_EI(ux, uy, uz)   # I close to 1
#'
#' @seealso [compute_SPI()], [align_morph_batch()]
#' @export
compute_EI <- function(ux, uy, uz) {
  n      <- length(ux)
  U      <- cbind(ux, uy, uz)
  T_mat  <- (t(U) %*% U) / n
  eig    <- eigen(T_mat, symmetric = TRUE)
  lambda <- sort(eig$values, decreasing = TRUE)
  lambda <- pmax(lambda, 0)
  
  data.frame(
    E       = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I       = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_),
    lambda1 = lambda[1],
    lambda2 = lambda[2],
    lambda3 = lambda[3]
  )
}