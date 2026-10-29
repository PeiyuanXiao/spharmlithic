# SPI, ER and IR for every core

Computes the Scar Pattern Index (SPI), the elongation ratio (ER) and the
isotropy ratio (IR) of every core in a scar table, each as a raw value,
a rarefied value and a random interval at the core's scar count. This is
the batch form of
[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md)
and
[`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md).

## Usage

``` r
spi_ei_batch(
  data,
  id_col = "ID",
  start_cols = c("Start_X", "Start_Y", "Start_Z"),
  end_cols = c("End_X", "End_Y", "End_Z"),
  min_length = 0,
  keep = NULL,
  raw = TRUE,
  rarefied = TRUE,
  random_interval = TRUE,
  m = 3,
  max_comb = 2000,
  seed = 2226,
  quiet = FALSE
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

- raw:

  Logical. Return the raw SPI (`SPI_raw`)? Default `TRUE`.

- rarefied:

  Logical. Return the rarefied SPI (`SPI_rarefied`)? Default `TRUE`.

- random_interval:

  Logical. Return the random interval at the core's scar count
  (`SPI_lower`, `SPI_upper`) and the position of the raw SPI relative to
  it (`SPI_vs_random`)? Default `TRUE`.

- m:

  Number of scars per combination for the rarefied value. Default 3.

- max_comb:

  Largest number of combinations averaged. Above it, `max_comb` distinct
  combinations are drawn at random (a Monte Carlo approximation of the
  mean over all combinations). Default 2000.

- seed:

  Seed for that random draw; the caller's random stream is restored
  afterwards. `NULL` draws from the current stream. Default 2226.

- quiet:

  Logical. Suppress the summary message? Default `FALSE`.

## Value

A tibble with one row per core, in order of first appearance: the ID
column, the `keep` columns, `n_scar` (scars used), and for each of SPI,
ER and IR the columns selected by `raw`, `rarefied` and
`random_interval`: `*_raw`, `*_rarefied`, `*_lower`, `*_upper` and
`*_vs_random` (`"within"` or `"outside"`).

## Details

Rows with missing coordinates or zero length are left out, and so are
scars shorter than `min_length`. Cores with fewer than three usable
scars get `NA` values. One message summarises what was left out and
which cores needed a Monte Carlo approximation of the rarefied values
(more than `max_comb` combinations).

Use the rarefied columns to compare cores or groups of cores (for
example with
[`compare_indices()`](https://peiyuanxiao.github.io/spharmlithic/reference/compare_indices.md)):
their expectation does not change with the scar count. Use the raw
values together with the random interval to judge a single core:
`"outside"` means the raw value departs from what random directions
produce at that core's scar count (SPI and ER above, IR below the
interval: more concentrated than random; SPI below the interval: typical
of opposed flaking). See
[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md)
for details.

For each core the rarefied SPI, ER and IR are averaged over the same
combinations of scars.

## See also

[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md),
[`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md),
[`random_interval()`](https://peiyuanxiao.github.io/spharmlithic/reference/random_interval.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md),
[`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md),
[`compare_indices()`](https://peiyuanxiao.github.io/spharmlithic/reference/compare_indices.md)

## Examples

``` r
set.seed(1)
scars <- do.call(rbind, lapply(1:4, function(i) {
  n <- 8 + 4 * i
  s <- matrix(rnorm(3 * n, sd = 10), n, 3)
  d <- matrix(rnorm(3 * n), n, 3)
  d[, 1] <- d[, 1] + i / 2
  data.frame(ID = paste0("core", i), Typology = if (i <= 2) "A" else "B",
             Start_X = s[, 1], Start_Y = s[, 2], Start_Z = s[, 3],
             End_X = s[, 1] + 10 * d[, 1], End_Y = s[, 2] + 10 * d[, 2],
             End_Z = s[, 3] + 10 * d[, 3])
}))
spi_ei_batch(scars, keep = "Typology")
#> 1 core(s) have more than 2000 combinations of 3 scars; their rarefied values average 2000 random combinations (seed 2226): core4.
#> # A tibble: 4 × 18
#>   ID    Typology n_scar SPI_raw SPI_rarefied SPI_lower SPI_upper SPI_vs_random
#>   <chr> <chr>     <int>   <dbl>        <dbl>     <dbl>     <dbl> <chr>        
#> 1 core1 A            12   0.430        0.613    0.0790     0.500 within       
#> 2 core2 A            16   0.540        0.685    0.0681     0.435 outside      
#> 3 core3 B            20   0.616        0.724    0.0607     0.390 outside      
#> 4 core4 B            24   0.750        0.824    0.0553     0.357 outside      
#> # ℹ 10 more variables: ER_raw <dbl>, ER_rarefied <dbl>, ER_lower <dbl>,
#> #   ER_upper <dbl>, ER_vs_random <chr>, IR_raw <dbl>, IR_rarefied <dbl>,
#> #   IR_lower <dbl>, IR_upper <dbl>, IR_vs_random <chr>

# The same cores without scars shorter than 8 units
spi_ei_batch(scars, min_length = 8, keep = "Typology", random_interval = FALSE)
#> Left out 4 scar(s) shorter than `min_length` (core1, core2, core3, core4).
#> # A tibble: 4 × 9
#>   ID    Typology n_scar SPI_raw SPI_rarefied ER_raw ER_rarefied IR_raw
#>   <chr> <chr>     <int>   <dbl>        <dbl>  <dbl>       <dbl>  <dbl>
#> 1 core1 A            11   0.510        0.657  0.250       0.512  0.498
#> 2 core2 A            15   0.514        0.669  0.187       0.520  0.558
#> 3 core3 B            19   0.600        0.714  0.535       0.623  0.273
#> 4 core4 B            23   0.745        0.820  0.613       0.677  0.181
#> # ℹ 1 more variable: IR_rarefied <dbl>
```
