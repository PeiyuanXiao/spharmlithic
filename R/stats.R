# ==============================================================================
# stats.R
# Directional statistics for 3-D scar orientation data:
#   - compute_SPI       : Scar Pattern Index (Clarkson 2013)
#   - compute_spi_angle : SPI converted to expected pairwise angle
#   - compute_EI        : Elongation / Isotropy from the orientation tensor
# ==============================================================================

#' Scar Pattern Index (Clarkson 2013)
#'
#' Computes the ratio of the resultant vector magnitude to the total scar
#' length, on a scale from 0 (random orientation) to 1 (perfect alignment).
#'
#' @param dx,dy,dz Numeric vectors of equal length. The X, Y, Z components
#'   of the scar direction vectors. These are typically unit direction
#'   vectors returned by [align_scar_batch()] (columns `d_x`, `d_y`, `d_z`).
#' @param lengths Optional numeric vector of scar lengths. If provided,
#'   each direction vector is scaled by its length before summation,
#'   reproducing Clarkson's (2013) length-weighted SPI. If `NULL` (the
#'   default), all scars contribute equally regardless of length.
#'
#' @return A single numeric value in \eqn{[0, 1]}. Values close to 1
#'   indicate strong preferred orientation; values close to 0 indicate
#'   isotropic / random patterning.
#'
#' @details
#' Clarkson's original SPI uses raw scar displacement vectors (start-to-end),
#' so longer scars contribute proportionally more to the resultant. When
#' `lengths` is `NULL`, this function computes an unweighted variant in
#' which every scar contributes equally regardless of length — useful when
#' scar lengths are unreliable or when one is interested only in direction.
#'
#' To reproduce Clarkson's original length-weighted definition, supply scar
#' lengths via the `lengths` argument:
#' \preformatted{
#'   lens <- sqrt((aligned$e_x - aligned$s_x)^2 +
#'                (aligned$e_y - aligned$s_y)^2 +
#'                (aligned$e_z - aligned$s_z)^2)
#'   compute_SPI(aligned$d_x, aligned$d_y, aligned$d_z, lengths = lens)
#' }
#'
#' @references
#' Clarkson, C. (2013). Measuring core reduction using 3D flake scar
#' density: A test case of changing core reduction at Shum Laka
#' (NW Cameroon).
#'
#' @examples
#' # Perfectly aligned vectors along X
#' compute_SPI(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))   # 1
#'
#' # Random-like distribution (close to 0)
#' set.seed(1)
#' compute_SPI(rnorm(50), rnorm(50), rnorm(50))
#'
#' @seealso [compute_spi_angle()]
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


#' Convert SPI to a scar-pattern angle (Clarkson 2013)
#'
#' Converts the Scar Pattern Index into the expected pairwise angle
#' between two randomly-selected scars, following Clarkson's (2013)
#' interpretation: \eqn{\theta = \arccos(\mathrm{SPI})}.
#'
#' SPI = 1 maps to 0 degrees (parallel scars), SPI = 0 maps to 90 degrees
#' (uniformly random pairwise angles).
#'
#' @inheritParams compute_SPI
#' @param unit Either `"degrees"` (default) or `"radians"`.
#'
#' @return A single numeric value in \eqn{[0, 90]} degrees (or
#'   \eqn{[0, \pi/2]} radians).
#'
#' @examples
#' compute_spi_angle(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))   # 0 (parallel)
#' compute_spi_angle(c(1, 0),    c(0, 1),    c(0, 0))      # 90 (orthogonal pair)
#'
#' @references
#' Clarkson, C. (2013). Measuring core reduction using 3D flake scar
#' density: A test case of changing core reduction at Shum Laka
#' (NW Cameroon).
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
#' Constructs the 3x3 orientation (fabric) tensor from unit direction
#' vectors, computes its eigenvalues
#' \eqn{\lambda_1 \ge \lambda_2 \ge \lambda_3}, and returns the shape
#' descriptors:
#' \deqn{E = 1 - \lambda_2 / \lambda_1, \quad I = \lambda_3 / \lambda_1.}
#'
#' @param ux,uy,uz Numeric vectors of equal length. The X, Y, Z components
#'   of **unit** direction vectors.
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
#' (degenerate tensor).
#'
#' The one-row data frame return type is designed to play well with
#' [dplyr::group_modify()], allowing direct per-specimen computation.
#'
#' @examples
#' ux <- c(1, 1, 0.9)
#' uy <- c(0, 0, 0.1)
#' uz <- c(0, 0, 0.0)
#' compute_EI(ux, uy, uz)
#'
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