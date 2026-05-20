# ==============================================================================
# spharm_utils.R
# Pure-R helpers for working with spharm_from_*() results, plus the
# inverse-transform helper spharm_reconstruct().
# Only spharm_reconstruct() touches Python; the rest are pure R.
# ==============================================================================

#' Flatten SPHARM result list into a wide-format data frame
#'
#' Converts the per-specimen list returned by [spharm_from_directions()] or
#' [spharm_from_meshes()] into a single wide-format data frame, with one
#' row per specimen and one column per power-spectrum value and per
#' coefficient. Suitable for CSV export and downstream multivariate
#' analysis (PCA, UMAP, clustering) in R.
#'
#' @param x A list returned by [spharm_from_directions()] or
#'   [spharm_from_meshes()].
#' @param include_coeffs Logical. If `TRUE` (default), include all
#'   flattened coefficients as columns named `coeff_0001`, `coeff_0002`,
#'   ... If `FALSE`, return only the power spectrum (much smaller).
#'
#' @return A `tibble` with columns:
#' \describe{
#'   \item{ID}{Specimen identifier.}
#'   \item{power_l0, power_l1, ..., power_lN}{Raw power per spherical
#'     harmonic degree (N = lmax).}
#'   \item{coeff_0001, coeff_0002, ...}{Flattened real-valued coefficients
#'     (only when `include_coeffs = TRUE`). For [spharm_from_meshes()],
#'     which returns complex coefficients with negligible imaginary parts,
#'     only the real part is retained.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- spharm_from_directions(my_aligned_data, lmax = 20)
#'
#' # Full output (power + coefficients)
#' df_full <- spharm_to_dataframe(result)
#'
#' # Power spectrum only (smaller, easier to inspect)
#' df_power <- spharm_to_dataframe(result, include_coeffs = FALSE)
#'
#' write.csv(df_full, "spharm_results.csv", row.names = FALSE)
#' }
#'
#' @seealso [spharm_from_directions()], [spharm_from_meshes()],
#'   [spharm_reconstruct()]
#'
#' @importFrom tibble tibble as_tibble
#' @export
spharm_to_dataframe <- function(x, include_coeffs = TRUE) {
  
  if (!is.list(x) || length(x) == 0) {
    stop("`x` must be a non-empty list as returned by spharm_from_*().",
         call. = FALSE)
  }
  
  ids <- names(x)
  if (is.null(ids) || any(ids == "")) {
    stop("`x` must be a named list (one element per specimen ID).",
         call. = FALSE)
  }
  
  # Power-spectrum columns --------------------------------------------------
  power_mat <- do.call(rbind, lapply(x, function(s) s$power_spectrum))
  colnames(power_mat) <- paste0("power_l", seq_len(ncol(power_mat)) - 1L)
  
  out <- tibble::tibble(ID = ids)
  out <- cbind(out, as.data.frame(power_mat))
  
  # Coefficient columns -----------------------------------------------------
  if (include_coeffs) {
    coeff_mat <- do.call(rbind, lapply(x, function(s) {
      cf <- s$coefficients
      if (is.complex(cf)) {
        # Track A (STL) returns complex coefficients with theoretically
        # zero imaginary parts (real-valued spherical functions).
        # Explicitly take Re() to avoid spurious "discarding imaginary
        # parts" warnings from as.numeric().
        cf <- Re(cf)
      }
      as.numeric(cf)
    }))
    n_coef    <- ncol(coeff_mat)
    coef_w    <- nchar(as.character(n_coef))  # zero-pad width
    colnames(coeff_mat) <- sprintf(paste0("coeff_%0", coef_w, "d"),
                                   seq_len(n_coef))
    out <- cbind(out, as.data.frame(coeff_mat))
  }
  
  tibble::as_tibble(out)
}


#' Inverse spherical harmonic transform: coefficients to spherical density
#'
#' Given a single specimen's spherical-harmonic coefficient array (as
#' returned in `spharm_from_*()$<id>$coefficients`), reconstructs the
#' corresponding density on a Driscoll-Healy grid on the unit sphere.
#' The returned object also includes the longitude / latitude vectors and
#' an `(n_lat * n_lon, 3)` matrix of unit-sphere Cartesian coordinates,
#' for direct plotting with `plotly::plot_ly()`.
#'
#' @param coefficients Numeric or complex array of shape
#'   `(2, lmax+1, lmax+1)` — the `coefficients` element of one specimen's
#'   `spharm_from_*()` result. Both Track-A (complex) and Track-B (real)
#'   coefficients are accepted.
#' @param grid_size Integer. Latitude resolution of the reconstruction
#'   grid; longitude resolution is `2 * grid_size`. Default 64. Larger
#'   values give smoother visualisations at higher cost.
#'
#' @return A list with:
#' \describe{
#'   \item{density}{Numeric matrix of shape `(grid_size, 2 * grid_size)`
#'     - the reconstructed density on the Driscoll-Healy grid. Negative
#'     values from finite-degree truncation are clipped to zero.}
#'   \item{lon}{Numeric vector of length `2 * grid_size`. Longitude
#'     (radians, range `[0, 2*pi)`).}
#'   \item{lat}{Numeric vector of length `grid_size`. Latitude (radians,
#'     range `(-pi/2, pi/2)`; equator = 0).}
#'   \item{xyz}{Numeric matrix of shape `(grid_size * 2 * grid_size, 3)`
#'     - unit-sphere Cartesian coordinates for each grid cell, row-major
#'     order matching `as.vector(density)`. Convenient input for
#'     `plotly::plot_ly(type = "surface", surfacecolor = ...)`.}
#' }
#'
#' @references
#' Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working
#' with spherical harmonics. \emph{Geochemistry, Geophysics, Geosystems},
#' \strong{19}(8), 2574--2592.
#'
#' @examples
#' \dontrun{
#' result <- spharm_from_directions(my_aligned_data, lmax = 20)
#'
#' # Reconstruct the spherical density for one specimen
#' rec <- spharm_reconstruct(result[[1]]$coefficients, grid_size = 64)
#'
#' # Quick plotly visualisation
#' library(plotly)
#' nlat <- length(rec$lat); nlon <- length(rec$lon)
#' x_mat <- matrix(rec$xyz[, "x"], nlat, nlon)
#' y_mat <- matrix(rec$xyz[, "y"], nlat, nlon)
#' z_mat <- matrix(rec$xyz[, "z"], nlat, nlon)
#' plot_ly(x = x_mat, y = y_mat, z = z_mat,
#'         surfacecolor = rec$density,
#'         type = "surface", colorscale = "Hot",
#'         showscale = TRUE)
#'
#' # Batch reconstruction with lapply()
#' all_recs <- lapply(result, function(s) spharm_reconstruct(s$coefficients))
#' }
#'
#' @seealso [spharm_from_directions()], [spharm_from_meshes()]
#'
#' @export
spharm_reconstruct <- function(coefficients, grid_size = 64) {
  
  if (is.null(coefficients) || length(dim(coefficients)) != 3L) {
    stop("`coefficients` must be a 3-D array of shape (2, lmax+1, lmax+1).",
         call. = FALSE)
  }
  
  # If complex, take real part (Track A returns complex with ~0 imag).
  if (is.complex(coefficients)) {
    coefficients <- Re(coefficients)
  }
  
  # Call Python: inverse SH transform onto a Driscoll-Healy grid.
  density <- sh_py$kde_to_spharm$reconstruct_from_coeffs(
    coefficients,
    as.integer(grid_size)
  )
  density <- as.matrix(density)
  
  n_lat <- nrow(density)
  n_lon <- ncol(density)
  
  # DH grid sample positions (pyshtools convention: endpoint=FALSE).
  colat <- seq(0, pi,        length.out = n_lat + 1L)[-1L]
  lon   <- seq(0, 2 * pi,    length.out = n_lon + 1L)[-1L]
  
  # Build matched unit-sphere Cartesian grid, row-major to match
  # as.vector(density).
  TH <- matrix(colat, n_lat, n_lon, byrow = FALSE)
  PH <- matrix(lon,   n_lat, n_lon, byrow = TRUE)
  xyz <- cbind(
    x = as.vector(sin(TH) * cos(PH)),
    y = as.vector(sin(TH) * sin(PH)),
    z = as.vector(cos(TH))
  )
  
  list(
    density = density,
    lon     = lon,
    lat     = pi / 2 - colat,
    xyz     = xyz
  )
}