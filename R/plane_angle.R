# ==============================================================================
# plane_angle.R
# 形态平面与技术平面的偏差分析
# ==============================================================================

#' Compute angle between morphological plane and technical (SVD) plane
#'
#' For a single group (e.g. one specimen or one stratigraphic layer), this
#' function calculates the angle (in degrees) between the morphological
#' best-fit plane (provided by Geomagic, stored as `Norm_X/Y/Z`) and the
#' technical best-fit plane derived from SVD of the scar direction vectors.
#'
#' @param df_group A data frame for a **single** group. Required columns:
#' \describe{
#'   \item{Norm_X, Norm_Y, Norm_Z}{Morphological plane normal (first row used).}
#'   \item{Start_X, Start_Y, Start_Z}{Scar start-point coordinates.}
#'   \item{End_X, End_Y, End_Z}{Scar end-point coordinates.}
#' }
#'
#' @return A one-row `data.frame` with column:
#' \describe{
#'   \item{angle_deg}{Angle between the two plane normals (degrees), rounded
#'     to 2 decimal places.}
#' }
#' Returns `NULL` if fewer than three valid scar vectors are present.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' result <- my_data %>%
#'   group_by(ID) %>%
#'   group_modify(~ compute_plane_angle(.x))
#' }
#'
#' @seealso [compute_scar_plane()]
#' @export
compute_plane_angle <- function(df_group) {
  
  # --- 形态最优拟合平面（Geomagic 法线）---
  normal_geo <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal_geo <- normal_geo / sqrt(sum(normal_geo^2))
  
  # --- 技术最优拟合平面（刮痕向量 SVD）---
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  if (sum(valid) < 3) return(NULL)
  
  U          <- cbind(dx[valid] / len[valid],
                      dy[valid] / len[valid],
                      dz[valid] / len[valid])
  normal_svd <- svd(U)$v[, 3]
  normal_svd <- normal_svd / sqrt(sum(normal_svd^2))
  
  # 统一法线方向（确保两法线在同侧）
  if (sum(normal_geo * normal_svd) < 0) normal_svd <- -normal_svd
  
  # --- 夹角（度）---
  cos_angle <- min(1, max(-1, sum(normal_geo * normal_svd)))
  angle_deg <- acos(cos_angle) * 180 / pi
  
  data.frame(angle_deg = round(angle_deg, 2))
}


#' Compute mean perpendicular distance of scar endpoints from morphological plane
#'
#' For a single group (e.g. one specimen or one stratigraphic layer), this
#' function calculates the mean perpendicular distance of all scar end-points
#' from the morphological best-fit plane.
#'
#' @param df_group A data frame for a **single** group. Required columns:
#' \describe{
#'   \item{Norm_X, Norm_Y, Norm_Z}{Morphological plane normal (first row used).}
#'   \item{Pos_X, Pos_Y, Pos_Z}{A reference point on the morphological plane
#'     (first row used).}
#'   \item{Start_X, Start_Y, Start_Z}{Scar start-point coordinates.}
#'   \item{End_X, End_Y, End_Z}{Scar end-point coordinates.}
#' }
#'
#' @return A one-row `data.frame` with column:
#' \describe{
#'   \item{mean_dist}{Mean perpendicular distance of scar end-points from the
#'     morphological plane, rounded to 4 decimal places.}
#' }
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' result <- my_data %>%
#'   group_by(ID) %>%
#'   group_modify(~ compute_scar_plane(.x))
#' }
#'
#' @seealso [compute_plane_angle()]
#' @export
compute_scar_plane <- function(df_group) {
  
  # --- 形态最优拟合平面（Geomagic 法线）---
  normal_geo <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal_geo <- normal_geo / sqrt(sum(normal_geo^2))
  p0 <- as.numeric(df_group[1, c("Pos_X", "Pos_Y", "Pos_Z")])
  
  # --- 刮痕端点到形态平面的平均垂直距离 ---
  endpoints <- rbind(
    as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]),
    as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )])
  )
  dist_to_plane <- abs(
    (endpoints - matrix(p0, nrow(endpoints), 3, byrow = TRUE)) %*% normal_geo
  )
  mean_dist <- mean(dist_to_plane)
  
  data.frame(mean_dist = round(mean_dist, 4))
}