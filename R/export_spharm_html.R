# ==============================================================================
# export_spharm_html.R
# Export interactive 3D HTML viewer for SPHARM reconstruction results
# ==============================================================================

#' Export interactive SPHARM reconstruction viewer as a self-contained HTML file
#'
#' Generates a standalone HTML file with an interactive Three.js-based 3D
#' viewer for exploring spherical harmonic reconstructions. The viewer
#' supports dual viewports (morphology and scar direction side-by-side),
#' degree-by-degree animation, multiple material presets, radial deviation
#' colormap, type-mean overlay, and OBJ/PNG export.
#'
#' Spherical harmonic synthesis is performed entirely in JavaScript on the
#' client side — the HTML file embeds only the coefficients as JSON,
#' keeping file sizes small.
#'
#' @param morph Result from [spharm_from_meshes()], or `NULL` if only scar
#'   direction data is available. A named list where each element contains
#'   `coefficients` (a `(2, lmax+1, lmax+1)` array) and `power_spectrum`.
#' @param scar Result from [spharm_from_directions()], or `NULL` if only
#'   morphology data is available. Same structure as `morph`.
#' @param meta Optional data frame with specimen metadata. Must contain an
#'   `ID` column matching the names in `morph` and/or `scar`. A `Typology`
#'   column, if present, enables the Type Mean toggle and groups the
#'   specimen dropdown. All other columns are shown in the Info panel.
#' @param out_path Character. Output file path (should end in `.html`).
#' @param lmax Integer. Maximum spherical harmonic degree to include.
#'   Default 20. Must not exceed the `lmax` used in the original analysis.
#' @param title Character. Display title in the viewer header. Default
#'   `"spharmlithic"`.
#' @param digits Integer. Number of decimal places to retain for
#'   coefficients in the JSON payload. Lower values reduce file size.
#'   Default 8.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return The output file path (invisibly), for use in pipelines.
#'
#' @details
#' At least one of `morph` or `scar` must be provided. When both are
#' supplied, specimens are matched by ID; a specimen appearing in only
#' one dataset will show "No data" in the other viewport.
#'
#' The generated HTML is fully self-contained (no local server needed)
#' and loads Three.js r128 from the cdnjs CDN. An internet connection is
#' required when first opening the file.
#'
#' @section Viewer controls:
#' \describe{
#'   \item{Specimen dropdown}{Select specimen; grouped by Typology if
#'     available.}
#'   \item{Degree slider / Play}{Reconstruct at a specific max degree,
#'     or animate from l=1 to lmax.}
#'   \item{Material selector}{Ceramic, Clay, Glass, Brushed Metal, X-Ray,
#'     or Flat.}
#'   \item{View presets}{Iso, Top, Front camera angles.}
#'   \item{Wire}{Toggle wireframe overlay.}
#'   \item{Info}{Show specimen metadata panel.}
#'   \item{Colormap}{Radial deviation blue-white-red colormap.}
#'   \item{Type Mean}{Average coefficients across specimens of the same
#'     Typology.}
#'   \item{PNG / OBJ}{Export screenshot or mesh at current degree.}
#' }
#'
#' @examples
#' \dontrun{
#' # Track B only
#' aligned <- align_scar_batch(my_scar_data)
#' scar_result <- spharm_from_directions(aligned, lmax = 20)
#' export_spharm_html(scar = scar_result, out_path = "viewer.html")
#'
#' # Both tracks with metadata
#' morph_result <- spharm_from_meshes("data/stl_files", lmax = 20)
#' meta <- data.frame(
#'   ID = names(scar_result),
#'   Typology = c("Levallois", "Discoid", "Levallois")
#' )
#' export_spharm_html(
#'   morph = morph_result,
#'   scar  = scar_result,
#'   meta  = meta,
#'   out_path = "spharm_viewer.html"
#' )
#' }
#'
#' @seealso [spharm_from_directions()], [spharm_from_meshes()],
#'   [spharm_reconstruct()]
#'
#' @export
export_spharm_html <- function(
    morph    = NULL,
    scar     = NULL,
    meta     = NULL,
    out_path,
    lmax     = 20L,
    title    = "spharmlithic",
    digits   = 8L,
    verbose  = TRUE) {
  
  # ---- Validate inputs -----------------------------------------------------
  if (is.null(morph) && is.null(scar)) {
    stop("At least one of `morph` or `scar` must be provided.", call. = FALSE)
  }
  
  lmax <- as.integer(lmax)
  if (lmax < 1L) stop("`lmax` must be >= 1.", call. = FALSE)
  
  # Collect all specimen IDs
  ids_morph <- if (!is.null(morph)) names(morph) else character(0)
  ids_scar  <- if (!is.null(scar))  names(scar)  else character(0)
  all_ids   <- sort(unique(c(ids_morph, ids_scar)))
  
  if (length(all_ids) == 0L) {
    stop("No specimen IDs found. Results must be named lists.", call. = FALSE)
  }
  
  # Validate meta
  if (!is.null(meta)) {
    if (!"ID" %in% names(meta)) {
      stop("`meta` must contain an 'ID' column.", call. = FALSE)
    }
    meta$ID <- as.character(meta$ID)
  }
  
  # ---- Build specimen records as JSON strings --------------------------------
  if (verbose) message("Building specimen data for ", length(all_ids),
                       " specimen(s)...")
  
  json_parts <- character(0)
  
  for (i in seq_along(all_ids)) {
    sid <- all_ids[i]
    
    # Extract coefficients as JSON strings
    morph_json <- "null"
    scar_json  <- "null"
    
    if (sid %in% ids_morph) {
      cj <- .cilm_to_json(morph[[sid]]$coefficients, lmax, digits)
      if (!is.null(cj)) morph_json <- cj
    }
    
    if (sid %in% ids_scar) {
      cj <- .cilm_to_json(scar[[sid]]$coefficients, lmax, digits)
      if (!is.null(cj)) scar_json <- cj
    }
    
    # Skip if neither track has data
    if (morph_json == "null" && scar_json == "null") {
      if (verbose) message("  [", i, "/", length(all_ids), "] ", sid,
                           " -- skipped (no data)")
      next
    }
    
    # Build metadata JSON
    meta_pairs <- character(0)
    if (!is.null(meta)) {
      meta_row <- meta[meta$ID == sid, , drop = FALSE]
      if (nrow(meta_row) > 0L) {
        mr <- meta_row[1L, ]
        for (col in setdiff(names(mr), "ID")) {
          val <- mr[[col]]
          if (!is.na(val)) {
            # Escape special JSON characters in value
            escaped <- gsub("\\\\", "\\\\\\\\", as.character(val))
            escaped <- gsub('"', '\\\\"', escaped)
            meta_pairs <- c(meta_pairs,
                            paste0('"', col, '":"', escaped, '"'))
          }
        }
      }
    }
    meta_json <- paste0("{", paste(meta_pairs, collapse = ","), "}")
    
    # Escape specimen ID for JSON
    escaped_id <- gsub("\\\\", "\\\\\\\\", sid)
    escaped_id <- gsub('"', '\\\\"', escaped_id)
    
    # Assemble record JSON
    record_json <- paste0(
      '{"id":"', escaped_id, '"',
      ',"morph":', morph_json,
      ',"scar":', scar_json,
      ',"meta":', meta_json,
      '}'
    )
    
    json_parts <- c(json_parts, record_json)
    if (verbose) message("  [", i, "/", length(all_ids), "] ", sid, " OK")
  }
  
  if (length(json_parts) == 0L) {
    stop("No valid specimens to export.", call. = FALSE)
  }
  
  # ---- Assemble final JSON array --------------------------------------------
  if (verbose) message("Serializing ", length(json_parts),
                       " specimen(s) to JSON...")
  data_json <- paste0("[", paste(json_parts, collapse = ","), "]")
  
  # ---- Read HTML template and substitute ------------------------------------
  template_path <- system.file("templates", "spharm_viewer.html",
                               package = "spharmlithic")
  
  # Fallback for devtools::load_all() during development
  if (template_path == "") {
    dev_path <- file.path(
      getNamespaceInfo("spharmlithic", "path"),
      "inst", "templates", "spharm_viewer.html"
    )
    if (file.exists(dev_path)) {
      template_path <- dev_path
    } else {
      stop("Cannot find HTML template. Is the package installed correctly?",
           call. = FALSE)
    }
  }
  
  html <- paste(readLines(template_path, encoding = "UTF-8"), collapse = "\n")
  
  html <- gsub("{{DATA_JSON}}",  data_json, html, fixed = TRUE)
  html <- gsub("{{LMAX}}",       as.character(lmax), html, fixed = TRUE)
  html <- gsub("{{GROUP_NAME}}", title, html, fixed = TRUE)
  
  # ---- Write output ---------------------------------------------------------
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_path, useBytes = TRUE)
  
  size_mb <- file.info(out_path)$size / (1024 * 1024)
  if (verbose) {
    message("Exported: ", out_path)
    message("  ", length(json_parts), " specimen(s), ",
            sprintf("%.1f", size_mb), " MB")
  }
  
  invisible(out_path)
}


# ---- Internal helpers -------------------------------------------------------

#' Convert a cilm coefficient array to a JSON string
#'
#' Manually builds a JSON nested array `[[[m00,m01,...],[m10,...],...],...]]`
#' from an R 3-D array of shape `(2, lmax+1, lmax+1)`.
#' This avoids `jsonlite` flattening or restructuring the nested array.
#'
#' @param coefficients A `(2, L+1, L+1)` numeric or complex array.
#' @param lmax_out Maximum degree to retain (may be <= original lmax).
#' @param digits Rounding digits.
#' @return A JSON string representing the 3-D array, or `NULL` on failure.
#' @noRd
.cilm_to_json <- function(coefficients, lmax_out, digits) {
  if (is.null(coefficients)) return(NULL)
  
  # Handle complex coefficients (Track A)
  if (is.complex(coefficients)) {
    coefficients <- Re(coefficients)
  }
  
  d <- dim(coefficients)
  if (is.null(d) || length(d) != 3L || d[1] != 2L) return(NULL)
  
  original_lmax <- d[2] - 1L
  use_lmax <- min(lmax_out, original_lmax)
  
  # Subset to requested lmax
  sub <- coefficients[, seq_len(use_lmax + 1L), seq_len(use_lmax + 1L),
                      drop = FALSE]
  sub <- round(sub, digits)
  
  # Build JSON string: cilm[i][l][m]
  # i = 0,1 (cosine, sine); l = 0..lmax; m = 0..lmax
  layers <- vapply(1:2, function(i) {
    rows <- vapply(seq_len(use_lmax + 1L), function(l) {
      vals <- as.numeric(sub[i, l, ])
      paste0("[", paste(vals, collapse = ","), "]")
    }, character(1))
    paste0("[", paste(rows, collapse = ","), "]")
  }, character(1))
  
  paste0("[", paste(layers, collapse = ","), "]")
}