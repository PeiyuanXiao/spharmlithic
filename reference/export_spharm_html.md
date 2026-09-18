# Export interactive SPHARM reconstruction viewer as a self-contained HTML file

Generates a standalone HTML file with an interactive Three.js-based 3D
viewer for exploring spherical harmonic reconstructions. The viewer
supports dual viewports (morphology and scar direction side-by-side),
degree-by-degree animation, multiple material presets, radial deviation
colormap, type-mean overlay, and OBJ/PNG export.

## Usage

``` r
export_spharm_html(
  morph = NULL,
  scar = NULL,
  meta = NULL,
  out_path,
  lmax = 20L,
  title = "spharmlithic",
  digits = 8L,
  verbose = TRUE
)
```

## Arguments

- morph:

  Result from
  [`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md),
  or `NULL` if only scar direction data is available. A named list where
  each element contains `coefficients` (a `(2, lmax+1, lmax+1)` array)
  and `power_spectrum`.

- scar:

  Result from
  [`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md),
  or `NULL` if only morphology data is available. Same structure as
  `morph`.

- meta:

  Optional data frame with specimen metadata. Must contain an `ID`
  column matching the names in `morph` and/or `scar`. A `Typology`
  column, if present, enables the Type Mean toggle and groups the
  specimen dropdown. All other columns are shown in the Info panel.

- out_path:

  Character. Output file path (should end in `.html`).

- lmax:

  Integer. Maximum spherical harmonic degree to include. Default 20.
  Must not exceed the `lmax` used in the original analysis.

- title:

  Character. Display title in the viewer header. Default
  `"spharmlithic"`.

- digits:

  Integer. Number of decimal places to retain for coefficients in the
  JSON payload. Lower values reduce file size. Default 8.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

The output file path (invisibly), for use in pipelines.

## Details

Spherical harmonic synthesis is performed entirely in JavaScript on the
client side — the HTML file embeds only the coefficients as JSON,
keeping file sizes small.

At least one of `morph` or `scar` must be provided. When both are
supplied, specimens are matched by ID; a specimen appearing in only one
dataset will show "No data" in the other viewport.

The generated HTML is fully self-contained (no local server needed) and
loads Three.js r128 from the cdnjs CDN. An internet connection is
required when first opening the file.

## Viewer controls

- Specimen dropdown:

  Select specimen; grouped by Typology if available.

- Degree slider / Play:

  Reconstruct at a specific max degree, or animate from l=1 to lmax.

- Material selector:

  Ceramic, Clay, Glass, Brushed Metal, X-Ray, or Flat.

- View presets:

  Iso, Top, Front camera angles.

- Wire:

  Toggle wireframe overlay.

- Info:

  Show specimen metadata panel.

- Colormap:

  Radial deviation blue-white-red colormap.

- Type Mean:

  Average coefficients across specimens of the same Typology.

- PNG / OBJ:

  Export screenshot or mesh at current degree.

## See also

[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md),
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md),
[`spharm_reconstruct()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_reconstruct.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Scar pattern only
aligned <- align_scar_batch(my_scar_data)
scar_result <- spharm_from_directions(aligned, lmax = 20)
export_spharm_html(scar = scar_result, out_path = "viewer.html")

# Both tracks with metadata
morph_result <- spharm_from_meshes("data/stl_files", lmax = 20)
meta <- data.frame(
  ID = names(scar_result),
  Typology = c("Levallois", "Discoid", "Levallois")
)
export_spharm_html(
  morph = morph_result,
  scar  = scar_result,
  meta  = meta,
  out_path = "spharm_viewer.html"
)
} # }
```
