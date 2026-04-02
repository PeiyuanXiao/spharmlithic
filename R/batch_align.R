# ==============================================================================
# batch_align.R
# 批量对齐包装函数 — 直接传入完整数据框，一步完成
# ==============================================================================

#' Batch SVD alignment for all specimens
#'
#' A convenience wrapper around [align_scar()] that accepts a complete
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
#' @examples
#' \dontrun{
#' aligned <- align_svd_batch(my_scar)
#' }
#'
#' @seealso [align_scar()], [align_morph_batch()]
#'
#' @importFrom dplyr group_by group_modify ungroup
#' @export
align_scar_batch <- function(data, id_col = "ID") {
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::group_modify(~ align_svd(.x)) %>%
    dplyr::ungroup()
}


#' Batch Lin 2024 alignment for all specimens
#'
#' A convenience wrapper around [align_morph()] that accepts a complete
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
#' @examples
#' \dontrun{
#' aligned <- align_lin2024_batch(my_scar)
#' }
#'
#' @seealso [align_morph()], [align_scar_batch()]
#'
#' @importFrom dplyr group_by group_modify ungroup
#' @export
align_morph_batch <- function(data, id_col = "ID") {
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::group_modify(~ align_lin2024(.x)) %>%
    dplyr::ungroup()
}
