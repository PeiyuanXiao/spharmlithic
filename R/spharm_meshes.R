# ==============================================================================
# spharm_meshes.R
# Track A: STL files -> open3d simplify -> trimesh smooth -> normalize ->
#          PCA align -> spherical interpolation -> SPHARM
# ==============================================================================

#' Spherical harmonic analysis from STL meshes
#'
#' Pipeline (Track A), per STL file:
#' \enumerate{
#'   \item Read mesh via `open3d`. Files exceeding 3M faces are
#'     pre-decimated by streaming face sub-sampling to avoid memory pressure.
#'   \item Quadric-error decimation to `target_faces`.
#'   \item Laplacian smoothing (`smooth_iterations` passes).
#'   \item Volume-centroid normalization to a unit sphere.
#'   \item Area-weighted PCA alignment with energy-based sign convention.
#'   \item Cartesian-to-spherical conversion of the aligned vertices.
#'   \item Interpolation onto a regular `(grid_size, grid_size)` lat-lon grid.
#'   \item Spherical harmonic expansion via `pyshtools` (4pi normalization,
#'     zero-component normalized: all coefficients divided by `c(0, 0)`).
#' }
#'
#' Requires the **mesh** Python extension. Run
#' `install_spharmlithic_python(mesh = TRUE)` once before first use.
#'
#' @param stl_dir Character. Path to a directory containing `.stl` files.
#'   Each STL is treated as one specimen; the file basename (without
#'   extension) becomes its ID.
#' @param lmax Integer. Maximum spherical harmonic degree. Default 20.
#' @param target_faces Integer. Decimation target. Default 20000.
#' @param grid_size Integer. Latitude resolution of the interpolation
#'   grid; the longitude direction uses `2 * grid_size` points
#'   (Driscoll-Healy sampling 2). Default 256.
#' @param smooth_iterations Integer. Laplacian smoothing iterations after
#'   decimation. Default 3. Set to 0 to skip smoothing.
#' @param pre_decimate_threshold Integer. Face count above which streaming
#'   pre-decimation is triggered. Default 3000000.
#' @param pre_decimate_target Integer. Pre-decimation target face count.
#'   Default 500000.
#' @param verbose Logical. Print per-file progress. Default `TRUE`.
#'
#' @return A list with one element per successfully-processed STL, named
#'   by file basename. Each element is a list with `coefficients` and
#'   `power_spectrum` (same structure as [spharm_from_directions()]).
#'   Failed specimens are reported via warnings and omitted from the
#'   result.
#'
#' @references
#' Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working
#' with spherical harmonics. \emph{Geochemistry, Geophysics, Geosystems},
#' \strong{19}(8), 2574--2592.
#'
#' @examples
#' \dontrun{
#' # Requires mesh extension
#' install_spharmlithic_python(mesh = TRUE)
#'
#' result <- spharm_from_meshes(
#'   stl_dir      = "data/3D_models",
#'   lmax         = 20,
#'   target_faces = 20000
#' )
#'
#' df <- spharm_to_dataframe(result)
#' write.csv(df, "spharm_meshes.csv", row.names = FALSE)
#' }
#'
#' @seealso [spharm_from_directions()], [spharm_to_dataframe()], [export_spharm_html()],
#'   [install_spharmlithic_python()]
#'
#' @export
spharm_from_meshes <- function(
    stl_dir,
    lmax                    = 20,
    target_faces            = 20000,
    grid_size               = 256,
    smooth_iterations       = 3,
    pre_decimate_threshold  = 3000000,
    pre_decimate_target     = 500000,
    verbose                 = TRUE) {
  
  if (!dir.exists(stl_dir)) {
    stop("Directory does not exist: ", stl_dir, call. = FALSE)
  }
  
  # Pre-flight: warn early if mesh extension isn't installed.
  if (!reticulate::py_module_available("trimesh") ||
      !reticulate::py_module_available("open3d")) {
    stop("STL pipeline requires `trimesh` and `open3d`. Install with:\n",
         "  install_spharmlithic_python(mesh = TRUE)",
         call. = FALSE)
  }
  
  py_result <- sh_py$api$spharm_from_meshes(
    stl_dir                 = normalizePath(stl_dir, mustWork = TRUE),
    lmax                    = as.integer(lmax),
    target_faces            = as.integer(target_faces),
    grid_size               = as.integer(grid_size),
    smooth_iterations       = as.integer(smooth_iterations),
    pre_decimate_threshold  = as.integer(pre_decimate_threshold),
    pre_decimate_target     = as.integer(pre_decimate_target),
    verbose                 = as.logical(verbose)
  )
  
  # Convert. Failed specimens come back as Python None and are skipped.
  result <- list()
  for (id in names(py_result)) {
    x <- py_result[[id]]
    if (is.null(x)) {
      warning("Specimen '", id, "' failed; see backend log above.",
              call. = FALSE)
      next
    }
    result[[id]] <- list(
      coefficients   = as.array(x$coefficients),
      power_spectrum = as.numeric(x$power_spectrum)
    )
  }
  result
}