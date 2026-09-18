# Pairwise PERMANOVA comparisons

Perform pairwise PERMANOVA comparisons among all levels of a grouping
variable using
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).

## Usage

``` r
pairwise_permanova(
  x,
  group,
  permutations = 999,
  p_adjust_method = "holm",
  distance = "euclidean"
)
```

## Arguments

- x:

  A numeric matrix or data frame containing observations in rows and
  variables in columns (e.g., ILR coordinates).

- group:

  A grouping variable with one value per row of `x`.

- permutations:

  Number of permutations passed to
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).
  Default is 999.

- p_adjust_method:

  Method used by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  is `"holm"`.

- distance:

  Distance metric passed to
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html). Default is
  `"euclidean"`.

## Value

A data frame with one row per pairwise comparison containing:

- group1:

  First group.

- group2:

  Second group.

- pair:

  Comparison label.

- n1:

  Sample size of group1.

- n2:

  Sample size of group2.

- R2:

  PERMANOVA R-squared.

- F:

  Pseudo-F statistic.

- p:

  Permutation p-value.

- p_adj:

  Adjusted p-value.

## Details

For each pair of groups, a distance matrix is calculated from the
supplied multivariate data and tested using PERMANOVA. P-values are
adjusted for multiple comparisons.

This function assumes that compositional data have already been
transformed to ILR coordinates. Euclidean distances computed on ILR
coordinates are equivalent to Aitchison distances in the original
compositional space.

## Examples

``` r
set.seed(123)

x <- matrix(rnorm(45), ncol = 3)
grp <- rep(c("A", "B", "C"), each = 5)

pairwise_permanova(x, grp, permutations = 99)
#>   group1 group2   pair n1 n2         R2         F    p p_adj
#> 1      A      B A vs B  5  5 0.11549097 1.0445657 0.35     1
#> 2      A      C A vs C  5  5 0.01391272 0.1128721 0.93     1
#> 3      B      C B vs C  5  5 0.05431665 0.4594912 0.69     1
```
