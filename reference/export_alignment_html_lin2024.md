# Export an interactive Lin 2024 alignment HTML page

Builds three-panel Lin 2024 alignment visualisations for every specimen
in `raw_data` and writes a standalone interactive HTML file to
`out_path`.

## Usage

``` r
export_alignment_html_lin2024(raw_data, out_path)
```

## Arguments

- raw_data:

  A data frame with columns `ID`, `Start_X/Y/Z`, `End_X/Y/Z`,
  `Norm_X/Y/Z`, `Pos_X/Y/Z`, and optionally `Length`.

- out_path:

  Character. File path for the output HTML.

## Value

Invisibly, the `out_path` string.

## Details

The exported HTML is fully self-contained (Plotly loaded from CDN) and
requires no R session to view. All specimen panels are serialised to
JSON at export time; switching specimens in the browser is instant.

For the SVD-based equivalent, see
[`export_alignment_html_svd()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_svd.md).

## See also

[`export_alignment_html_svd()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_svd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
export_alignment_html_lin2024(raw_data, "output/alignment_lin2024.html")
} # }
```
