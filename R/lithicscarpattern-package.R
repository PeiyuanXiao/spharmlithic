#' lithicscarpattern: Lithic Flaking Scar Orientation Alignment and Analysis
#'
#' @description
#' `lithicscarpattern` provides a complete toolkit for the geometric alignment
#' and statistical analysis of 3-D lithic flaking scar orientation data.
#'
#' ## Main components
#'
#' ### Alignment pipelines
#' * [align_scar()] — Three-step SVD pipeline (rotate → translate → in-plane rotate).
#' * [align_morph()] — Two-step Lin 2024 pipeline (morphological normal rotation →
#'   longest-scar translation).
#'
#' ### Plane deviation analysis
#' * [compute_plane_angle()] — Angle between morphological and SVD planes.
#' * [compute_scar_plane()] — Mean perpendicular distance of scar endpoints from morphological plane.
#'
#' ### Scar pattern statistics
#' * [compute_SPI()] — Scar Pattern Index (Clarkson method).
#' * [compute_EI()] — Elongation (E) and Isotropy (I) from the orientation tensor.
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
#' * [build_panels_morph()] — Lin 2024 three-panel figure for one specimen.
#'
#' ### Geometry utilities
#' * [get_rot_matrix()] — Rodrigues rotation matrix between two unit vectors.
#' * [get_scar_length()] — Scar lengths from a data frame.
#'
#' @keywords internal
"_PACKAGE"
