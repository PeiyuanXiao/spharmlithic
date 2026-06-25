## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = FALSE
)

# run these lines here because the block we show the reader has to be
# eval=FALSE for GitHub Actions because the workflows config file already 
# runs the Python env setup functions
library(spharmlithic)
library(readxl)

## ----install, eval = FALSE----------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("PeiyuanXiao/spharmlithic")

## ----install-python, eval = FALSE---------------------------------------------
# library(spharmlithic)
# 
# # Core install (scar pattern analysis only)
# install_spharmlithic_python()
# 
# # Full install (also enables morphological analysis from meshes)
# install_spharmlithic_python(mesh = TRUE)

## ----load-package, eval = FALSE, message = FALSE------------------------------
# library(spharmlithic)
# library(readxl)
# 
# # point reticulate to the conda environment
# use_spharmlithic_python("r-spharmlithic")

## ----example-data, eval = TRUE------------------------------------------------
scar_path <- system.file("extdata", "example_scars.xlsx",
                          package = "spharmlithic")
example   <- readxl::read_excel(scar_path)
names(example)

## ----mesh-dir, eval = TRUE----------------------------------------------------
stl_dir <- system.file("extdata", "meshes", package = "spharmlithic")
list.files(stl_dir, pattern = "\\.stl$", ignore.case = TRUE)

## ----load-scar-data, eval = TRUE----------------------------------------------
scar_path <- system.file("extdata", "example_scars.xlsx",
                          package = "spharmlithic")
raw_data  <- read_excel(scar_path)
head(raw_data)

## ----align-svd, eval = TRUE---------------------------------------------------
aligned_svd <- align_scar_batch(raw_data)
head(aligned_svd)

## ----align-lin, eval = TRUE---------------------------------------------------
aligned_lin <- align_morph_batch(raw_data)
head(aligned_lin)

## ----export-html-svd----------------------------------------------------------
# export_alignment_html_svd(raw_data,
#                           out_path = "alignment_svd.html")
# 
# export_alignment_html_lin2024(raw_data,
#                               out_path = "alignment_lin2024.html")

## ----spi, eval = TRUE---------------------------------------------------------
spi <- compute_SPI(aligned_svd$d_x, 
                   aligned_svd$d_y, 
                   aligned_svd$d_z)
spi

## ----spi-weighted, eval = TRUE------------------------------------------------
lens <- get_scar_length(aligned_svd)
spi_w <- compute_SPI(aligned_svd$d_x, 
                     aligned_svd$d_y, 
                     aligned_svd$d_z,
                     lengths = lens)
spi_w

## ----spi-angle, eval = TRUE---------------------------------------------------
spi_angle <- compute_spi_angle(aligned_svd$d_x,
                               aligned_svd$d_y,
                               aligned_svd$d_z)
spi_angle

## ----ei, eval = TRUE----------------------------------------------------------
ei <- compute_EI(aligned_svd$d_x, 
                 aligned_svd$d_y, 
                 aligned_svd$d_z)
ei

## ----per-specimen, eval = FALSE-----------------------------------------------
# library(dplyr)
# aligned_svd %>%
#   group_by(ID) %>%
#   summarise(
#     SPI       = compute_SPI(d_x, d_y, d_z),
#     SPI_angle = compute_spi_angle(d_x, d_y, d_z),
#     compute_EI(d_x, d_y, d_z)
#   )

## ----spharm-scar--------------------------------------------------------------
# sh_scar <- spharm_from_directions(aligned_svd, lmax = 20, bandwidth = 0.35)
# sh_scar

## ----df-scar------------------------------------------------------------------
# df_scar <- spharm_to_dataframe(sh_scar)
# head(df_scar)

## ----spharm-morph-------------------------------------------------------------
# stl_dir  <- system.file("extdata", "meshes", package = "spharmlithic")
# sh_morph <- spharm_from_meshes(stl_dir, lmax = 20)
# sh_morph

## ----df-morph-----------------------------------------------------------------
# df_morph <- spharm_to_dataframe(sh_morph)
# head(df_morph)

## ----lmax-morph---------------------------------------------------------------
# sh_fine <- spharm_from_meshes(stl_dir, lmax = 50)

## ----viewer-scar-only---------------------------------------------------------
# # Scar pattern only
# export_spharm_html(
#   scar     = sh_scar,
#   out_path = "viewer_scar.html"
# )

## ----viewer-both--------------------------------------------------------------
# # Both morphology and scar pattern side-by-side
# export_spharm_html(
#   morph    = sh_morph,
#   scar     = sh_scar,
#   out_path = "spharm_viewer.html"
# )

## ----viewer-meta--------------------------------------------------------------
# meta <- data.frame(
#   ID       = c("ClarksonEXP01_Levallois", "ClarksonEXP02_Discoid"),
#   Typology = c("Levallois", "Discoid")
# )
# export_spharm_html(
#   scar     = sh_scar,
#   meta     = meta,
#   out_path = "spharm_viewer.html",
#   title    = "Example specimens"
# )

## ----reconstruct--------------------------------------------------------------
# recon <- spharm_reconstruct(
#   sh_scar[["ClarksonEXP01_Levallois"]]$coefficients,
#   grid_size = 64
# )

## ----plotly-recon-------------------------------------------------------------
# library(plotly)
# nlat <- length(recon$lat)
# nlon <- length(recon$lon)
# x_mat <- matrix(recon$xyz[, "x"], nlat, nlon)
# y_mat <- matrix(recon$xyz[, "y"], nlat, nlon)
# z_mat <- matrix(recon$xyz[, "z"], nlat, nlon)
# 
# plot_ly(x = x_mat, y = y_mat, z = z_mat,
#         surfacecolor = recon$density,
#         type = "surface", colorscale = "Hot",
#         showscale = TRUE)

