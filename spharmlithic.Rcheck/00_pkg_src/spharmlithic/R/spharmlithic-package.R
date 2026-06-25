#' spharmlithic: SPHARM framework for stone artifact shape and scar pattern analysis
#'
#' @description
#' `spharmlithic` provides a complete toolkit for the geometric alignment,
#' spherical harmonic decomposition for artifact shape and flaking scar vector
#'
#' ## Main components
#'
#' ### Alignment pipelines
#' * [align_scar_batch()] — Three-step SVD pipeline
#'   (rotate -> translate -> in-plane rotate).
#' * [align_morph_batch()] — Two-step Lin 2024 pipeline
#'   (morphological normal rotation -> longest-scar translation).
#'
#' ### Classical orientation statistics
#' * [compute_SPI()] — Scar Pattern Index (Clarkson et al. 2006); supports
#'   both unit-vector and length-weighted variants.
#' * [compute_spi_angle()] — SPI converted to expected pairwise angle
#'   (Clarkson et al. 2006).
#' * [compute_EI()] — Elongation (E) and Isotropy (I) from the
#'   orientation tensor (Lin et al. 2024).
#'
#' ### Spherical harmonic analysis (Python-backed)
#' * [install_spharmlithic_python()] — One-time Python backend setup.
#' * [use_spharmlithic_python()] — Point at an existing Python env.
#' * [spharm_from_directions()] — Direction vectors -> vMF KDE -> SPHARM.
#' * [spharm_from_meshes()] — STL meshes -> spherical interpolation -> SPHARM.
#' * [spharm_to_dataframe()] — Flatten SH results for CSV export and
#'   downstream multivariate analysis (PCA, clustering).
#' * [spharm_reconstruct()] — Low-level inverse transform: coefficients ->
#'   density grid (for custom analysis or plotting).
#'
#' ### Interactive viewers and HTML export
#' * [export_spharm_html()] — Export an interactive Three.js-based 3D
#'   viewer for spherical harmonic reconstructions (morphology, scar
#'   direction, or both side-by-side).
#' * [export_alignment_html_svd()] — Four-panel SVD alignment page.
#' * [export_alignment_html_lin2024()] — Three-panel Lin 2024 alignment page.
#'
#' ### Geometry utilities
#' * [get_rot_matrix()] — Rodrigues rotation matrix between two unit vectors.
#' * [get_scar_length()] — Scar lengths from a data frame.
#'
#' @keywords internal
"_PACKAGE"