# Spatiodirectional SP-SPHARM features of every core

Describes the scar pattern of each core as a whole, with no feature
specified in advance, by 22 rotation-invariant features computed from
the removal directions of its scars, the positions of their starts, and
the correspondence between the two (Xiao et al., in prep.; extending
SP-SPHARM of Xiao et al., 2027). Use the features to compare the overall
scar patterns of groups of cores with
[`compare_groups()`](https://peiyuanxiao.github.io/spharmlithic/reference/compare_groups.md).

## Usage

``` r
spharm_spatiodirection(
  data,
  id_col = "ID",
  start_cols = c("Start_X", "Start_Y", "Start_Z"),
  end_cols = c("End_X", "End_Y", "End_Z"),
  min_length = 0,
  keep = NULL,
  features = c("dir", "start", "same", "diff"),
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

- features:

  Which blocks of features to return; any of `"dir"`, `"start"`,
  `"same"`, `"diff"` (default: all four, 22 features).

- quiet:

  Logical. Suppress the summary message? Default `FALSE`.

## Value

A tibble with one row per core: the ID column, the `keep` columns,
`n_scar` (scars used) and the features. The names of the feature columns
are stored in the attribute `"features"`. Cores with fewer than three
usable scars get `NA`.

## Details

Each scar \\i\\ is represented by two unit vectors: its removal
direction \\\hat v_i\\, and its start vector \\\hat u_i\\, the direction
from the centroid of all the core's starts to its start. For a degree
\\l\\ the features average the Legendre polynomial \\P_l\\ of the cosine
between two vectors:

- `dir_l` (l = 1..6), direction-direction pairs: the mean of \\P_l(\hat
  v_i \cdot \hat v_j)\\ over the distinct pairs of scars \\i \ne j\\.
  These describe the directional relationships among the scars; `dir_1`
  is the mean cosine between the directions of two different scars.

- `start_l` (l = 2..6), start-start pairs: the same for the start
  vectors, describing how the starts are spread around the core.

- `same_l` (l = 1..6), within-scar start-direction pairs: the mean of
  \\P_l(\hat v_i \cdot \hat u_i)\\ over the scars, describing how each
  scar's direction relates to its own start position (for example
  whether scars run inward from where they start).

- `diff_l` (l = 2..6), between-scar direction-start pairs: the mean of
  \\P_l(\hat v_i \cdot \hat u_j)\\ over the ordered pairs \\i \ne j\\,
  describing how a scar's direction relates to the starts of the other
  scars.

Each feature is multiplied by \\\sqrt{(2l+1)/2}\\, which gives every
degree the same variance under random directions, so the features can be
compared with Euclidean distances without further scaling (do **not**
standardise them again, for example in
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)).

**Scar count.** Pairing each scar with itself contributes a term that
depends only on the scar count, so those self-pairs are left out: every
feature is a mean over single scars or distinct pairs of scars (a
U-statistic), and its expectation does not change with the scar count.
For the same reason `start_1` and `diff_1`, which contain a term of
order \\1/(n-1)\\, are not returned. Compare cores with different scar
counts directly; do not add the scar count as a covariate.

The features use only angles between the vectors of a core, so they do
not change when the core is rotated or moved: raw and aligned
coordinates give the same features.

## References

Steinhardt, P. J., Nelson, D. R., & Ronchetti, M. (1983).
Bond-orientational order in liquids and glasses. *Physical Review B*,
**28**(2), 784–805.

Hoeffding, W. (1948). A class of statistics with asymptotically normal
distribution. *Annals of Mathematical Statistics*, **19**(3), 293–325.

## See also

[`compare_groups()`](https://peiyuanxiao.github.io/spharmlithic/reference/compare_groups.md),
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md)

## Examples

``` r
set.seed(1)
scars <- do.call(rbind, lapply(1:3, function(i) {
  s <- matrix(rnorm(30, sd = 10), 10, 3)
  e <- s + matrix(rnorm(30, sd = 5), 10, 3)
  data.frame(ID = paste0("core", i), Start_X = s[, 1], Start_Y = s[, 2],
             Start_Z = s[, 3], End_X = e[, 1], End_Y = e[, 2], End_Z = e[, 3])
}))
sp <- spharm_spatiodirection(scars)
sp
#> # A tibble: 3 × 24
#>   ID    n_scar   dir_1    dir_2   dir_3   dir_4   dir_5   dir_6  start_2 start_3
#>   <chr>  <int>   <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
#> 1 core1     10 -0.0615 -0.00227  0.102  -0.0536  0.0390  0.123  -0.00803  0.0776
#> 2 core2     10 -0.0358  0.0397   0.0493  0.0377 -0.124  -0.0817  0.0887   0.0593
#> 3 core3     10  0.126  -0.129   -0.0155 -0.0478  0.0755  0.0139  0.329    0.106 
#> # ℹ 14 more variables: start_4 <dbl>, start_5 <dbl>, start_6 <dbl>,
#> #   same_1 <dbl>, same_2 <dbl>, same_3 <dbl>, same_4 <dbl>, same_5 <dbl>,
#> #   same_6 <dbl>, diff_2 <dbl>, diff_3 <dbl>, diff_4 <dbl>, diff_5 <dbl>,
#> #   diff_6 <dbl>
attr(sp, "features")
#>  [1] "dir_1"   "dir_2"   "dir_3"   "dir_4"   "dir_5"   "dir_6"   "start_2"
#>  [8] "start_3" "start_4" "start_5" "start_6" "same_1"  "same_2"  "same_3" 
#> [15] "same_4"  "same_5"  "same_6"  "diff_2"  "diff_3"  "diff_4"  "diff_5" 
#> [22] "diff_6" 

# Direction features only
spharm_spatiodirection(scars, features = "dir")
#> # A tibble: 3 × 8
#>   ID    n_scar   dir_1    dir_2   dir_3   dir_4   dir_5   dir_6
#>   <chr>  <int>   <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 core1     10 -0.0615 -0.00227  0.102  -0.0536  0.0390  0.123 
#> 2 core2     10 -0.0358  0.0397   0.0493  0.0377 -0.124  -0.0817
#> 3 core3     10  0.126  -0.129   -0.0155 -0.0478  0.0755  0.0139
```
