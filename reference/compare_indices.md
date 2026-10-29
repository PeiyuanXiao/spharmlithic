# Compare groups of cores on single indices

Compares groups of cores on one index at a time, such as the rarefied
SPI, ER and IR of
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md)
or the CPI and PRA summaries of
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md).
Each index is tested for a difference in group means with Welch-type
tests, which allow the groups to have unequal variances, and every pair
of groups gets a standardised mean difference with a bootstrap
confidence interval.

## Usage

``` r
compare_indices(
  data,
  cols,
  group,
  alpha = 0.05,
  p_adjust = "holm",
  boot = 9999,
  seed = 2226
)
```

## Arguments

- data:

  A data frame with one row per core.

- cols:

  Names of the index columns to compare.

- group:

  The groups: the name of a column of `data`, or a vector with one value
  per row. The order of the factor levels (or of the sorted values) sets
  which group is `group1` in each pair.

- alpha:

  Significance level used in the notes. Default 0.05.

- p_adjust:

  Method for adjusting the pairwise p values within each index, passed
  to [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).
  Default `"holm"`.

- boot:

  Number of bootstrap resamples for the confidence interval of the
  standardised mean difference; 0 for none. Default 9999.

- seed:

  Seed set before the bootstrap of each pair; the caller's random stream
  is restored afterwards. `NULL` uses the current stream. Default 2226.

## Value

An object of class `spharmlithic_index_comparison`, a list with

- overall:

  A tibble with one row per index: `index`, `test` (`"Welch t"` for two
  groups, `"Welch ANOVA"` for more), `statistic`, `df1`, `df2` (`NA` for
  the t test), `p`, `n_groups` and `n`.

- pairwise:

  A tibble with one row per index and pair of groups: `index`, `group1`,
  `group2`, `n1`, `n2`, `mean1`, `mean2`, `smd` (standardised mean
  difference, group1 minus group2), `smd_lower`, `smd_upper` (95%
  bootstrap interval), `p` (Welch t test) and `p_adj`.

- settings:

  The arguments used.

## Details

Two groups are compared with Welch's t test, three or more with Welch's
one-way ANOVA
([`stats::oneway.test()`](https://rdrr.io/r/stats/oneway.test.html))
followed by Welch t tests between every pair, adjusted by `p_adjust`
within each index (not across indices, each index answering its own
question). Welch-type tests are used because differences in scar count
can leave groups with unequal variances.

The standardised mean difference is the difference between the group
means divided by the square root of the mean of the two group variances;
its absolute value is the effect size reported by Xiao et al. (in
prep.). The 95% interval is a percentile bootstrap stratified by group:
each group is resampled with replacement within itself.

**Which columns.** Raw SPI, ER and IR shift with the scar count; compare
groups with the rarefied values (`*_rarefied`). A message is given when
a `*_raw` column is passed. The mean and GMD of CPI and PRA do not need
rarefying.

Missing values are left out index by index, and so are groups with fewer
than two cores.

## References

Welch, B. L. (1951). On the comparison of several mean values: an
alternative approach. *Biometrika*, **38**(3/4), 330–336.

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics*, **6**(2), 65–70.

## See also

[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md),
[`compare_groups()`](https://peiyuanxiao.github.io/spharmlithic/reference/compare_groups.md)

## Examples

``` r
set.seed(1)
d <- data.frame(Typology = rep(c("A", "B", "C"), each = 10),
                SPI_rarefied = rnorm(30, rep(c(0.5, 0.6, 0.6), each = 10), 0.08),
                IR_rarefied = rnorm(30, 0.2, 0.05))
compare_indices(d, c("SPI_rarefied", "IR_rarefied"), "Typology", boot = 999)
#> Index comparison (Welch tests; pairwise p adjusted by holm within each index)
#> 
#>         index        test statistic df1   df2        p n_groups  n
#>  SPI_rarefied Welch ANOVA  6.086079   2 17.68 0.009747        3 30
#>   IR_rarefied Welch ANOVA  0.001611   2 17.18 0.998391        3 30
#> 
#> Pairwise: standardised mean difference (group1 - group2) with 95% bootstrap interval (999 resamples) 
#>         index group1 group2 n1 n2      smd smd_lower smd_upper        p   p_adj
#>  SPI_rarefied      A      B 10 10 -1.45968   -3.4402   -0.4970 0.004728 0.01418
#>  SPI_rarefied      A      C 10 10 -1.12794   -2.8638   -0.2478 0.021727 0.04345
#>  SPI_rarefied      B      C 10 10  0.37718   -0.4435    1.4689 0.410214 0.41021
#>   IR_rarefied      A      B 10 10 -0.01888   -0.9173    0.8941 0.966835 1.00000
#>   IR_rarefied      A      C 10 10 -0.02478   -0.9204    1.0702 0.956461 1.00000
#>   IR_rarefied      B      C 10 10 -0.01120   -0.8950    1.1002 0.980354 1.00000
```
