# Replace zeros in compositional data

Replaces zeros in each row of a compositional data matrix with a small
positive value and rescales the remaining non-zero entries so that each
row continues to sum to one.

## Usage

``` r
replace_zeros(x, delta = NULL, fraction = 0.65)
```

## Arguments

- x:

  A numeric matrix or data frame whose rows represent compositions that
  sum to one.

- delta:

  Optional numeric scalar giving the value used to replace zeros. If
  `NULL` (default), the replacement value is computed separately for
  each row as `fraction * min(non_zero_values)`.

- fraction:

  Numeric scalar in `(0, 1]` giving the fraction of the smallest
  non-zero value to use when calculating `delta`. Ignored if `delta` is
  supplied.

## Value

A numeric matrix of the same dimensions as `x`, with zeros replaced and
rows rescaled to sum to one.

## Details

This function is intended for compositional data that contain structural
zeros which would otherwise prevent log-ratio transformations such as
the centred log-ratio (CLR) transform.

By default, zeros are replaced with 65% of the smallest non-zero value
in each row.

## Examples

``` r
# Three compositions (rows sum to one); the first two contain a zero
x <- rbind(
  c(0.60, 0.30, 0.10, 0.00),
  c(0.50, 0.25, 0.25, 0.00),
  c(0.40, 0.35, 0.15, 0.10)
)

# Default: each zero becomes 65% of the smallest non-zero value in its row
replace_zeros(x)
#>         [,1]     [,2]     [,3]   [,4]
#> [1,] 0.56100 0.280500 0.093500 0.0650
#> [2,] 0.41875 0.209375 0.209375 0.1625
#> [3,] 0.40000 0.350000 0.150000 0.1000

# Use a fixed replacement value instead
replace_zeros(x, delta = 0.001)
#>        [,1]    [,2]    [,3]  [,4]
#> [1,] 0.5994 0.29970 0.09990 0.001
#> [2,] 0.4995 0.24975 0.24975 0.001
#> [3,] 0.4000 0.35000 0.15000 0.100

# Rows still sum to one
rowSums(replace_zeros(x))
#> [1] 1 1 1
```
