# ==============================================================================
# spharm_utils.R
# Pure-R helpers for working with spharm_from_*() results.
# These do NOT require Python.
# ==============================================================================

#' Flatten SPHARM result list into a wide-format data frame
#'
#' Converts the per-specimen list returned by [spharm_from_directions()] or
#' [spharm_from_meshes()] into a single wide-format data frame, with one row
#' per specimen and one column per power-spectrum value and per coefficient.
#' Suitable for CSV export and downstream multivariate analysis (PCA, UMAP,
#' clustering) in R.
#'
#' @param x A list returned by [spharm_from_directions()] or
#'   [spharm_from_meshes()].
#' @param include_coeffs Logical. If `TRUE` (default), include all flattened
#'   coefficients as columns named `coeff_0001`, `coeff_0002`, ... If
#'   `FALSE`, return only the power spectrum (much smaller).
#'
#' @return A `tibble` with columns:
#' \describe{
#'   \item{ID}{Specimen identifier.}
#'   \item{power_l0, power_l1, ..., power_lN}{Raw power per spherical
#'     harmonic degree (where N = lmax).}
#'   \item{coeff_0001, coeff_0002, ...}{Flattened coefficients (only when
#'     `include_coeffs = TRUE`). Order follows pyshtools convention:
#'     `clm[ic, l, m]` flattened as
#'     `c(clm[1,,], clm[2,,])` row-major. Use this column order
#'     consistently if you plan to round-trip via
#'     [reshape_coeffs_to_array()].}
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
#' @seealso [spharm_from_directions()], [spharm_from_meshes()]
#'
#' @importFrom tibble tibble
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

  # Power-spectrum columns ----------------------------------------------------
  power_mat <- do.call(rbind, lapply(x, function(s) s$power_spectrum))
  colnames(power_mat) <- paste0("power_l", seq_len(ncol(power_mat)) - 1L)

  out <- tibble::tibble(ID = ids)
  out <- cbind(out, as.data.frame(power_mat))

  # Coefficient columns -------------------------------------------------------
  if (include_coeffs) {
    coeff_mat <- do.call(rbind, lapply(x, function(s) as.numeric(s$coefficients)))
    n_coef    <- ncol(coeff_mat)
    coef_w    <- nchar(as.character(n_coef))   # zero-pad width
    colnames(coeff_mat) <- sprintf(paste0("coeff_%0", coef_w, "d"),
                                   seq_len(n_coef))
    out <- cbind(out, as.data.frame(coeff_mat))
  }

  tibble::as_tibble(out)
}
