# Export an interactive SVD alignment HTML page

Builds four-panel SVD alignment visualisations for every specimen in
`raw_data` and writes a standalone interactive HTML file to `out_path`.
The page contains a drop-down menu to switch between specimens.

## Usage

``` r
export_alignment_html_svd(raw_data, out_path)
```

## Arguments

- raw_data:

  A data frame with columns `ID`, `Start_X/Y/Z`, `End_X/Y/Z`.

- out_path:

  Character. File path for the output HTML (e.g.
  `"output/scar_alignment_svd.html"`).

## Value

Invisibly, the `out_path` string.

## Details

The exported HTML is fully self-contained (Plotly loaded from CDN) and
requires no R session to view. All specimen panels are serialised to
JSON at export time; switching specimens in the browser is instant.

For the Lin 2024 morphology-based equivalent, see
[`export_alignment_html_lin2024()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_lin2024.md).

## See also

[`export_alignment_html_lin2024()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_lin2024.md)

## Examples

``` r
if (FALSE) { # \dontrun{
export_alignment_html_svd(raw_data, "output/alignment_svd.html")
} # }
```
