# Export an interactive Lin 2024 alignment HTML page

Builds three-panel Lin 2024 alignment visualisations for every specimen
in `raw_data` and writes a standalone interactive HTML file to
`out_path`.

## Usage

``` r
export_alignment_html_lin2024(
  raw_data,
  out_path,
  id_col = "ID",
  outlines = NULL
)
```

## Arguments

- raw_data:

  A data frame with columns `Start_X/Y/Z`, `End_X/Y/Z`, `Norm_X/Y/Z`,
  `Pos_X/Y/Z`, the ID column named by `id_col`, and optionally `Length`.

- out_path:

  Character. File path for the output HTML.

- id_col:

  Name of the specimen ID column. Default `"ID"`.

- outlines:

  Optional named list of core outlines (named by ID), as returned by
  [`read_outline_iges()`](https://peiyuanxiao.github.io/spharmlithic/reference/read_outline_iges.md),
  or planes fitted with
  [`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md).
  Each outline is drawn with its centre and moved through the same
  alignment steps as the scars of its core.

## Value

Invisibly, the `out_path` string.

## Details

The exported HTML is fully self-contained (Plotly loaded from CDN) and
requires no R session to view. All specimen panels are serialised to
JSON at export time; switching specimens in the browser is instant.

For the SVD-based equivalent, see
[`export_alignment_html_svd()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_svd.md).

## See also

[`export_alignment_html_svd()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_svd.md),
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md)

## Examples

``` r
if (FALSE) { # \dontrun{
export_alignment_html_lin2024(raw_data, "output/alignment_lin2024.html")
} # }
```
