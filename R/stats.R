# ==============================================================================
# stats.R
# 方向统计量：Scar Pattern Index (SPI)、Elongation (E)、Isotropy (I)
# ==============================================================================

#' SPI (Clarkson method)
#'
#' Computes the SPI \eqn{R} for a set of 3-D scar direction
#' vectors following the Clarkson convention: the magnitude of the resultant
#' vector divided by the sum of all individual vector magnitudes.
#'
#' @param dx Numeric vector. X-components of the direction vectors.
#' @param dy Numeric vector. Y-components of the direction vectors.
#' @param dz Numeric vector. Z-components of the direction vectors.
#'
#' @return A single numeric value in \eqn{[0, 1]}.  Values close to 1
#'   indicate strong preferred orientation; values close to 0 indicate
#'   isotropy.
#'
#' @examples
#' # Perfectly aligned vectors along X
#' compute_SPI(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))  # returns 1
#'
#' @export
compute_SPI <- function(dx, dy, dz) {
  resultant_magnitude <- sqrt(sum(dx)^2 + sum(dy)^2 + sum(dz)^2)
  total_length        <- sum(sqrt(dx^2 + dy^2 + dz^2))
  resultant_magnitude / total_length
}


#' Elongation (E) and Isotropy (I) from the orientation tensor
#'
#' Constructs the 3×3 orientation (fabric) tensor from unit direction vectors,
#' computes its eigenvalues \eqn{\lambda_1 \ge \lambda_2 \ge \lambda_3}, and
#' returns the shape descriptors defined as:
#' \deqn{E = 1 - \lambda_2 / \lambda_1, \quad I = \lambda_3 / \lambda_1.}
#'
#' @param ux Numeric vector. X-components of unit direction vectors.
#' @param uy Numeric vector. Y-components of unit direction vectors.
#' @param uz Numeric vector. Z-components of unit direction vectors.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{E}{Elongation index \eqn{[0, 1]}: 1 = perfectly linear, 0 = planar/isotropic.}
#'   \item{I}{Isotropy index \eqn{[0, 1]}: 1 = perfectly isotropic, 0 = linear/planar.}
#'   \item{lambda1}{Largest eigenvalue.}
#'   \item{lambda2}{Middle eigenvalue.}
#'   \item{lambda3}{Smallest eigenvalue.}
#' }
#'
#' @details
#' Negative eigenvalues caused by numerical error are clamped to zero.
#' Both \code{E} and \code{I} are set to \code{NA} when
#' \eqn{\lambda_1 \approx 0} (degenerate tensor).
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

  list(
    E       = ifelse(lambda[1] > 1e-10, 1 - lambda[2] / lambda[1], NA_real_),
    I       = ifelse(lambda[1] > 1e-10,     lambda[3] / lambda[1], NA_real_),
    lambda1 = lambda[1],
    lambda2 = lambda[2],
    lambda3 = lambda[3]
  )
}
