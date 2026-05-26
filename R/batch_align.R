# ==============================================================================
# batch_align.R
# Batch alignment wrappers — process a complete multi-specimen data frame
# ==============================================================================

#' Batch SVD alignment for all specimens
#'
#' A convenience wrapper around `align_scar()` that accepts a complete
#' multi-specimen data frame and handles the grouping internally.
#' Equivalent to calling `group_by(ID) %>% group_modify(~ align_scar(.x))`.
#'
#' @param data A data frame containing **all** specimens. Must contain columns
#'   `ID`, `Start_X`, `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`.
#' @param id_col Character. Name of the specimen ID column (default `"ID"`).
#'
#' @return The input data frame with nine additional columns:
#' \describe{
#'   \item{s_x, s_y, s_z}{Aligned start-point coordinates.}
#'   \item{e_x, e_y, e_z}{Aligned end-point coordinates.}
#'   \item{d_x, d_y, d_z}{Aligned unit direction vectors.}
#' }
#'
#' @details
#' Specimens are processed independently: the SVD plane and in-plane
#' rotation are estimated separately for each group defined by `id_col`.
#' Row order within each group is preserved; the output is ungrouped.
#'
#' @examples
#' \dontrun{
#' aligned <- align_scar_batch(raw_data)
#' head(aligned[, c("ID", "d_x", "d_y", "d_z")])
#' }
#'
#' @seealso [align_morph_batch()]
#'
#' @importFrom dplyr group_by group_modify ungroup
#' @export
align_scar_batch <- function(data, id_col = "ID") {
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::group_modify(~ align_scar(.x)) %>%
    dplyr::ungroup()
}


#' Batch Lin 2024 alignment for all specimens
#'
#' A convenience wrapper around `align_morph()` that accepts a complete
#' multi-specimen data frame and handles the grouping internally.
#'
#' @param data A data frame containing **all** specimens. Must contain columns
#'   `ID`, `Start_X`, `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`,
#'   `Norm_X`, `Norm_Y`, `Norm_Z`, and optionally `Length`.
#' @param id_col Character. Name of the specimen ID column (default `"ID"`).
#'
#' @return The input data frame with nine additional columns
#'   (`s_x/y/z`, `e_x/y/z`, `d_x/y/z`).
#'
#' @details
#' Specimens are processed independently: the morphological plane normal
#' and the longest-scar anchor are determined separately for each group
#' defined by `id_col`.  Row order within each group is preserved; the
#' output is ungrouped.
#'
#' @examples
#' \dontrun{
#' aligned <- align_morph_batch(raw_data)
#' head(aligned[, c("ID", "d_x", "d_y", "d_z")])
#' }
#'
#' @seealso [align_scar_batch()]
#'
#' @importFrom dplyr group_by group_modify ungroup
#' @export
align_morph_batch <- function(data, id_col = "ID") {
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::group_modify(~ align_morph(.x)) %>%
    dplyr::ungroup()
}