# Export an interactive SVD alignment HTML page

Builds four-panel SVD alignment visualisations for every specimen in
`raw_data` and writes a standalone interactive HTML file to `out_path`.
The page contains a drop-down menu to switch between specimens.

## Usage

``` r
export_alignment_html_svd(raw_data, out_path, id_col = "ID", outlines = NULL)
```

## Arguments

- raw_data:

  A data frame with columns `Start_X/Y/Z`, `End_X/Y/Z` and the ID column
  named by `id_col`.

- out_path:

  Character. File path for the output HTML (e.g.
  `"output/scar_alignment_svd.html"`).

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

Alignment is for viewing and for
[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md);
the scar-pattern indices of
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md)
and
[`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md)
do not change when a core is rotated or moved, and do not need it.

For the Lin 2024 morphology-based equivalent, see
[`export_alignment_html_lin2024()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_lin2024.md).

## See also

[`export_alignment_html_lin2024()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_lin2024.md),
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md)

## Examples

``` r
if (FALSE) { # \dontrun{
export_alignment_html_svd(raw_data, "output/alignment_svd.html")

# With the outlines of the peripherally exploited cores
export_alignment_html_svd(raw_data, "output/alignment_svd.html",
                          outlines = read_outline_iges(files))
} # }
```
