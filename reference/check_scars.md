# Check a scar table before analysis

Reports in one message what the scar-pattern functions would make of a
scar table: how many cores and scars it holds, which scars would be left
out (missing coordinates, zero length, shorter than `min_length`), which
cores have fewer than three usable scars, and whether the `keep` columns
are constant within each core. Run it once after reading the data.

## Usage

``` r
check_scars(
  data,
  id_col = "ID",
  start_cols = c("Start_X", "Start_Y", "Start_Z"),
  end_cols = c("End_X", "End_Y", "End_Z"),
  min_length = 0,
  keep = NULL
)
```

## Arguments

- data:

  A data frame with one row per scar and the columns named by `id_col`,
  `start_cols` and `end_cols`. Scars of all cores are stacked.

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

- keep:

  Names of per-core columns (for example `"Typology"`) to carry into the
  output. Each must be constant within a core.

## Value

Invisibly, a tibble with one row per core: the ID column, the `keep`
columns, `n_rows` (rows in the table), `n_missing`, `n_zero`, `n_short`
(rows left out for each reason), `n_scar` (usable scars) and `enough`
(`n_scar` of at least 3).

## See also

[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md),
[`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md)

## Examples

``` r
d <- data.frame(ID = rep(c("a", "b"), c(4, 2)),
                Start_X = c(0, 1, 2, 3, 0, 1), Start_Y = 0, Start_Z = 0,
                End_X = c(1, 1, 3, 5, 1, 2), End_Y = c(1, 0, 1, 1, 1, 1),
                End_Z = 0)
check_scars(d)
#> 2 core(s), 6 scar row(s); usable scars per core: 2-3 (median 2.5).
#> Left out 1 scar(s) of zero length (a).
#> 1 core(s) have fewer than 3 usable scars and will be NA: b.
```
