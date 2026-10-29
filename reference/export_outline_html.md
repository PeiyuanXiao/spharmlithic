# Export every core and its outline to an interactive HTML page

Writes a self-contained HTML page with one 3-D view per core that has an
outline (see
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md))
and a drop-down menu to switch between them. Each view states the core's
frame check and flatness (see
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md)),
so that cores whose scars and outline do not line up can be found
quickly.

## Usage

``` r
export_outline_html(
  data,
  outlines,
  out_path,
  id_col = "ID",
  start_cols = c("Start_X", "Start_Y", "Start_Z"),
  end_cols = c("End_X", "End_Y", "End_Z"),
  min_length = 0
)
```

## Arguments

- data:

  A scar table with all cores.

- outlines:

  A named list of outlines, as for
  [`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md).

- out_path:

  Path of the HTML file to write.

- id_col:

  Name of the core ID column. Default `"ID"`.

- start_cols, end_cols:

  Names of the start- and end-point coordinate columns (X, Y, Z).
  Defaults `Start_X/Y/Z` and `End_X/Y/Z`, the columns the scars are
  recorded in. The aligned columns of
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)
  (`s_x`... and `e_x`...) give the same results: these indices do not
  change when a core is rotated or moved.

- min_length:

  Scars shorter than this (in the units of the coordinates) are left
  out. Default 0 (keep all). Use it to repeat an analysis without small
  scars, for example `min_length = 10` for 10 mm.

## Value

Invisibly, `out_path`.

## See also

[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
outlines <- read_outline_iges(list.files("outlines", full.names = TRUE))
export_outline_html(scars, outlines, "outline_check.html")
} # }
```
