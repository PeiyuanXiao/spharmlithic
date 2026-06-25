# ==============================================================================
# geometry_utils.R
# ==============================================================================
#' Compute a rotation matrix that maps vector \code{a} onto vector \code{b}
#'
#' Uses the Rodrigues rotation formula to construct the 3×3 rotation matrix
#' that rotates unit vector \code{a} into unit vector \code{b}.
#'
#' @param a Numeric vector of length 3. The source direction (need not be
#'   normalised; it will be normalised internally).
#' @param b Numeric vector of length 3. The target direction (need not be
#'   normalised).
#'
#' @return A 3×3 numeric rotation matrix.
#'
#' @details
#' Two degenerate cases are handled (threshold: \eqn{10^{-10}}):
#' \itemize{
#'   \item **Antiparallel** (\eqn{\cos\theta < -1 + 10^{-10}}): a 180°
#'     rotation about an arbitrary perpendicular axis is returned.
#'   \item **Parallel** (\eqn{\cos\theta > 1 - 10^{-10}}): the identity
#'     matrix is returned.
#' }
#'
#' This function is a low-level geometric utility; in typical usage it is
#' called internally by `align_scar()` and `align_morph()` rather than
#' directly.
#'
#' @examples
#' # Rotate the Z-axis onto the X-axis
#' R <- get_rot_matrix(c(0, 0, 1), c(1, 0, 0))
#' round(R %*% c(0, 0, 1), 10)   # should equal c(1, 0, 0)
#'
#' # Parallel vectors — identity matrix
#' R_id <- get_rot_matrix(c(0, 0, 1), c(0, 0, 2))
#' all.equal(R_id, diag(3))       # TRUE
#'
#' # Antiparallel vectors — 180-degree rotation
#' R_flip <- get_rot_matrix(c(0, 0, 1), c(0, 0, -1))
#' round(R_flip %*% c(0, 0, 1), 10)   # should equal c(0, 0, -1)
#'
#' @seealso [align_scar_batch()], [align_morph_batch()]
#' @export
get_rot_matrix <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  cos_theta <- sum(a * b)
  if (cos_theta < -1 + 1e-10) {
    perp <- if (abs(a[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
    v    <- perp - sum(perp * a) * a
    v    <- v / sqrt(sum(v^2))
    return(2 * outer(v, v) - diag(3))
  }
  if (cos_theta > 1 - 1e-10) return(diag(3))
  v <- c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
  v_skew <- matrix(
    c( 0,    -v[3],  v[2],
       v[3],  0,    -v[1],
       -v[2],  v[1],  0   ),
    3, 3, byrow = TRUE
  )
  diag(3) + v_skew + v_skew %*% v_skew * ((1 - cos_theta) / sum(v^2))
}