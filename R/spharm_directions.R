# ==============================================================================
# spharm_directions.R
# Track B: scar direction vectors -> vMF KDE -> Driscoll-Healy grid -> SPHARM
# ==============================================================================

#' Spherical harmonic analysis from scar direction vectors
#'
#' Pipeline (Track B):
#' \enumerate{
#'   \item Extract unit direction vectors `(d_x, d_y, d_z)` from each scar.
#'   \item Per-specimen von Mises-Fisher kernel density estimate on the sphere.
#'   \item Interpolate KDE values onto a Driscoll-Healy regular grid.
#'   \item Spherical harmonic expansion via `pyshtools` (4pi normalization).
#' }
#'
#' This function calls the bundled Python backend via `reticulate`. Make
#' sure you have run [install_spharmlithic_python()] (or
#' [use_spharmlithic_python()]) at least once.
#'
#' @param data A data frame containing direction vectors. Must contain an
#'   ID column (default `"ID"`) and three direction-component columns
#'   (default `"d_x"`, `"d_y"`, `"d_z"`, matching the output of
#'   [align_scar()] and [align_morph()]). Rows where a direction vector
#'   is zero or NA are silently dropped.
#' @param lmax Integer. Maximum spherical harmonic degree. Default 20.
#' @param bandwidth Numeric. vMF KDE bandwidth (smaller -> sharper peaks).
#'   Default 0.35. This is an empirical value; users are encouraged to
#'   conduct sensitivity analysis based on data characteristics
#'   (specimen scar count, directional concentration).
#' @param n_bearing,n_plunge Integers. Resolution of the intermediate
#'   evaluation grid for vMF KDE. Defaults 72 and 36 (5-degree spacing).
#' @param dh_size Integer. Driscoll-Healy grid latitude size; longitude is
#'   `2 * dh_size`. Default 64.
#' @param id_col,dx_col,dy_col,dz_col Column names. Defaults match
#'   [align_scar()] output.
#' @param verbose Logical. Print per-specimen progress. Default `TRUE`.
#'
#' @return A list with one element per specimen, named by ID. Each
#'   element is itself a list with:
#' \describe{
#'   \item{coefficients}{Numeric array of shape `(2, lmax+1, lmax+1)` -
#'     spherical harmonic coefficients in pyshtools 4pi normalization.
#'     Index `[1, l, m]` is the cosine coefficient (m >= 0); index
#'     `[2, l, m]` is the sine coefficient (m > 0).}
#'   \item{power_spectrum}{Numeric vector of length `lmax+1` - raw power
#'     per degree.}
#' }
#' Use [spharm_to_dataframe()] to flatten the list into a wide-format
#' data frame suitable for CSV export.
#'
#' @references
#' Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working
#' with spherical harmonics. \emph{Geochemistry, Geophysics, Geosystems},
#' \strong{19}(8), 2574--2592.
#'
#' @examples
#' \dontrun{
#' # Typical workflow
#' aligned <- align_scar_batch(my_scar_data)
#' result  <- spharm_from_directions(aligned, lmax = 20)
#'
#' # Inspect one specimen
#' result[[1]]$power_spectrum
#'
#' # Flatten for CSV / downstream analysis
#' df <- spharm_to_dataframe(result)
#' write.csv(df, "spharm_results.csv", row.names = FALSE)
#' }
#'
#' @seealso [spharm_from_meshes()], [spharm_to_dataframe()],
#'   [install_spharmlithic_python()]
#'
#' @importFrom reticulate r_to_py
#' @export
spharm_from_directions <- function(
    data,
    lmax       = 20,
    bandwidth  = 0.35,
    n_bearing  = 72,
    n_plunge   = 36,
    dh_size    = 64,
    id_col     = "ID",
    dx_col     = "d_x",
    dy_col     = "d_y",
    dz_col     = "d_z",
    verbose    = TRUE) {
  
  # ---- Validate input ---------------------------------------------------
  required <- c(id_col, dx_col, dy_col, dz_col)
  missing  <- setdiff(required, names(data))
  if (length(missing)) {
    stop("`data` is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  
  # Drop rows with zero / NA direction vectors
  d <- data[, c(id_col, dx_col, dy_col, dz_col)]
  names(d) <- c("ID", "ux", "uy", "uz")
  ok <- stats::complete.cases(d) &
    (d$ux^2 + d$uy^2 + d$uz^2) > 1e-20
  if (any(!ok)) {
    if (verbose) message("Dropping ", sum(!ok),
                         " row(s) with zero / NA direction vectors.")
    d <- d[ok, , drop = FALSE]
  }
  if (nrow(d) == 0) {
    stop("No valid direction vectors after cleaning.", call. = FALSE)
  }
  d$ID <- as.character(d$ID)
  
  # ---- Call Python backend ----------------------------------------------
  py_result <- sh_py$api$spharm_from_directions(
    df         = reticulate::r_to_py(d),
    lmax       = as.integer(lmax),
    bandwidth  = as.numeric(bandwidth),
    n_bearing  = as.integer(n_bearing),
    n_plunge   = as.integer(n_plunge),
    dh_size    = as.integer(dh_size),
    verbose    = as.logical(verbose)
  )
  
  # ---- Convert to native R structures -----------------------------------
  result <- lapply(py_result, function(x) {
    list(
      coefficients   = as.array(x$coefficients),
      power_spectrum = as.numeric(x$power_spectrum)
    )
  })
  names(result) <- names(py_result)
  result
}