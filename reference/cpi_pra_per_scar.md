# CPI and PRA of every scar of one core

Returns the Centripetality Index (CPI) and the Plane-Relative Angle
(PRA) of each scar of one core, for example to look at their
distribution. See
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md)
for the definitions.

## Usage

``` r
cpi_pra_per_scar(
  data,
  outline,
  id_col = "ID",
  start_cols = c("Start_X", "Start_Y", "Start_Z"),
  end_cols = c("End_X", "End_Y", "End_Z"),
  min_length = 0
)
```

## Arguments

- data:

  The scars of **one** core (a data frame, one row per scar). If it
  contains an `id_col` column, it must hold a single ID.

- outline:

  The outline of the core (see
  [`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md)),
  or a plane already fitted with
  [`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md).

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

A tibble with one row per usable scar: `scar` (the `Scar_ID` column if
present, otherwise the row number in `data`), `CPI` and `PRA` (degrees),
and `dist_to_outline` (the distance of the scar's start from the
outline, in the units of the coordinates).

## See also

[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md),
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md)

## Examples

``` r
t <- seq(0, 2 * pi, length.out = 13)[-13]
start <- cbind(30 * cos(t), 30 * sin(t), 0)
end <- start * 0.4
end[, 3] <- -5
one <- data.frame(Start_X = start[, 1], Start_Y = start[, 2], Start_Z = start[, 3],
                  End_X = end[, 1], End_Y = end[, 2], End_Z = end[, 3])
u <- seq(0, 2 * pi, length.out = 200)
cpi_pra_per_scar(one, cbind(30 * cos(u), 30 * sin(u), 0))
#> Warning: collapsing to unique 'x' values
#> Warning: collapsing to unique 'x' values
#> Warning: collapsing to unique 'x' values
#> # A tibble: 12 × 4
#>     scar   CPI   PRA dist_to_outline
#>    <int> <dbl> <dbl>           <dbl>
#>  1     1   180  15.5         0      
#>  2     2   180  15.5         0.00363
#>  3     3   180  15.5         0.00208
#>  4     4   180  15.5         0.00280
#>  5     5   180  15.5         0.00332
#>  6     6   180  15.5         0.00114
#>  7     7   180  15.5         0.00374
#>  8     8   180  15.5         0.00114
#>  9     9   180  15.5         0.00332
#> 10    10   180  15.5         0.00280
#> 11    11   180  15.5         0.00208
#> 12    12   180  15.5         0.00363
```
