# ==============================================================================
# utils-globals.R
# Suppress R CMD check notes for dplyr NSE variables and base functions.
# ==============================================================================

# Variables used inside dplyr::mutate() are not visible to R CMD check's
# static analysis. Declaring them here suppresses the "no visible binding"
# NOTE without affecting runtime behaviour.
utils::globalVariables(c(
  # align_scar() / align_morph() — dplyr::mutate() column names
  "s_x", "s_y", "s_z",
  "e_x", "e_y", "e_z"
))

#' @importFrom rlang .data
#' @importFrom stats dist
#' @importFrom utils globalVariables
NULL