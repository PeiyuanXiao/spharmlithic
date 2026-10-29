# Compare the overall scar patterns of groups of cores

Tests whether groups of cores differ in their scar patterns as a whole,
using a set of features such as those of
[`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md).
Three permutation tests on Euclidean distances are reported side by
side, for all groups together and for every pair of groups:

## Usage

``` r
compare_groups(
  x,
  group,
  features = NULL,
  permutations = 999,
  alpha = 0.05,
  p_adjust = "holm",
  permdisp_type = c("centroid", "median"),
  seed = 2226
)
```

## Arguments

- x:

  A data frame or numeric matrix with one row per core, such as the
  output of
  [`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md).

- group:

  The groups: the name of a column of `x`, or a vector with one value
  per row of `x`.

- features:

  Names of the columns of `x` to compare. Default: the features recorded
  by
  [`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md)
  (attribute `"features"`), else the columns named like its features
  (`dir_1`, `start_2`, ...), else all numeric columns except `n_scar`.

- permutations:

  Number of permutations of each test. Default 999.

- alpha:

  Significance level used for the notes. Default 0.05.

- p_adjust:

  Method for adjusting the pairwise p values within each test, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"holm"`.

- permdisp_type:

  Distances to the group `"centroid"` (default, the original PERMDISP)
  or to the spatial `"median"`; passed to
  [`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html).

- seed:

  Seed set before each test, so each is reproducible on its own; the
  caller's random stream is restored afterwards. `NULL` uses the current
  stream. Default 2226.

## Value

An object of class `spharmlithic_group_comparison`, a list with

- overall:

  A tibble with one row per test: `test`, `statistic` (pseudo-F for
  PERMANOVA and PERMDISP, \\W_d^\*\\), `R2` (PERMANOVA) and `p`.

- pairwise:

  With three or more groups, a tibble with one row per pair: `group1`,
  `group2`, `n1`, `n2`, the statistic and p value of each test, the
  adjusted p values (`*_p_adj`) and `note`; `NULL` with two groups.

- groups:

  The number of cores per group.

- features:

  The columns compared.

- settings:

  The arguments used.

Printing the object shows both tables and the notes.

## Details

- **PERMANOVA** (Anderson 2001): do the group centroids differ?

- **PERMDISP** (Anderson 2006): do the groups differ in dispersion, the
  spread of their cores around their centroid?

- the \\W_d^\*\\ test (Hamidi et al. 2019): a Welch-type test of the
  centroids that allows the groups to differ in dispersion.

When PERMDISP is significant, PERMANOVA may mistake a difference in
dispersion for a difference in location, especially with unequal group
sizes; the location difference should then be judged by \\W_d^\*\\. The
printed result points out where this applies.

Rows with a missing feature or group are left out, and so are groups
with fewer than two cores; a message says which. The features of
[`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md)
do not change with the scar count, so the scar count should not be added
as a covariate.

The \\W_d^\*\\ statistic is computed in the package and agrees with
`WdStar::WdS()` (Hamidi et al. 2019); its p value is the proportion of
permutations of the group labels, the observed one included, with a
statistic at least as large.

## References

Anderson, M. J. (2001). A new method for non-parametric multivariate
analysis of variance. *Austral Ecology*, **26**(1), 32–46.

Anderson, M. J. (2006). Distance-based tests for homogeneity of
multivariate dispersions. *Biometrics*, **62**(1), 245–253.

Hamidi, B., Wallace, K., Vasu, C., & Alekseyenko, A. V. (2019).
\\W_d^\*\\-test: robust distance-based multivariate analysis of
variance. *Microbiome*, **7**, 51.

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics*, **6**(2), 65–70.

## See also

[`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md),
[`compare_indices()`](https://peiyuanxiao.github.io/spharmlithic/reference/compare_indices.md),
[`pairwise_permanova()`](https://peiyuanxiao.github.io/spharmlithic/reference/pairwise_permanova.md)

## Examples

``` r
set.seed(1)
x <- data.frame(group = rep(c("A", "B", "C"), each = 8),
                f1 = rnorm(24, rep(c(0, 0, 1), each = 8)),
                f2 = rnorm(24, sd = rep(c(1, 3, 1), each = 8)))
compare_groups(x, "group", permutations = 199)
#> Group comparison of 2 features (Euclidean distance, 199 permutations)
#> Groups: A (n = 8), B (n = 8), C (n = 8) 
#> 
#> Overall
#>       test statistic     R2    p
#>  PERMANOVA     1.176 0.1007 0.35
#>   PERMDISP     4.358     NA 0.03
#>        W*d     2.017     NA 0.10
#> 
#> Pairwise (p values adjusted by holm within each test)
#>  group1 group2 n1 n2 PERMANOVA_p_adj PERMDISP_p_adj Wd_p_adj note
#>       A      B  8  8            0.99          0.105     1.00     
#>       A      C  8  8            0.15          0.865     0.09     
#>       B      C  8  8            0.46          0.105     0.44     
#> 
#> Note: PERMDISP is significant overall: the groups differ in dispersion,
#>       so PERMANOVA may confuse dispersion with location; judge the
#>       location difference by W*d.
```
