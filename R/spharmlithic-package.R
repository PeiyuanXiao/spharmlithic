#' spharmlithic: SPHARM framework for stone artifact shape and scar pattern analysis
#'
#' @description
#' `spharmlithic` provides a complete toolkit for the geometric alignment,
#' spherical harmonic decomposition for artifact shape and flaking scar vector
#'
#' ## Main components
#'
#' ### Alignment pipelines
#' * [align_scar()] / [align_scar_batch()] — Three-step SVD pipeline
#'   (rotate -> translate -> in-plane rotate).
#' * [align_morph()] / [align_morph_batch()] — Two-step Lin 2024 pipeline
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
#' * [spharm_reconstruct()] — Inverse transform: coefficients -> spherical
#'   density grid (for visualisation).
#' * [spharm_to_dataframe()] — Flatten results for CSV export.
#'
#' ### 3-D visualisation (Plotly helpers)
#' * [add_scars_3d()], [add_arrow_3d()], [add_tilted_plane_3d()],
#'   [add_plane_3d()], [make_scene()], [panel_layout()], [get_panel_json()]
#'
#' ### Interactive HTML export
#' * [export_alignment_html_svd()] — Four-panel SVD alignment page.
#' * [export_alignment_html_lin2024()] — Three-panel Lin 2024 alignment page.
#'
#' ### Visualisation builders
#' * [build_panel_scar()] — SVD four-panel figure for one specimen.
#' * [build_panel_morph()] — Lin 2024 three-panel figure for one specimen.
#'
#' ### Geometry utilities
#' * [get_rot_matrix()] — Rodrigues rotation matrix between two unit vectors.
#' * [get_scar_length()] — Scar lengths from a data frame.
#'
#' @keywords internal
"_PACKAGE"
