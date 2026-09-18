# Compute isometric log-ratio coordinates

Transform compositional data to isometric log-ratio (ILR) coordinates.

## Usage

``` r
make_ilr(x, delta = NULL, fraction = 0.65)
```

## Arguments

- x:

  A numeric matrix or data frame of compositional data.

- delta:

  Optional replacement value passed to
  [`replace_zeros()`](https://peiyuanxiao.github.io/spharmlithic/reference/replace_zeros.md).

- fraction:

  Fraction of the smallest non-zero value used by
  [`replace_zeros()`](https://peiyuanxiao.github.io/spharmlithic/reference/replace_zeros.md)
  when `delta = NULL`.

## Value

A data frame containing ILR coordinates.

## Details

Columns with zero variance are removed before transformation because
they do not contribute information and can cause numerical problems.

Zeros are replaced using
[`replace_zeros()`](https://peiyuanxiao.github.io/spharmlithic/reference/replace_zeros.md)
prior to transformation.

## Examples

``` r
x <- data.frame(
  a = c(0.5, 0.4, 0.6),
  b = c(0.3, 0.4, 0.2),
  c = c(0.2, 0.2, 0.2)
)

make_ilr(x)
#>        ilr_1
#> 1 -0.3612083
#> 2  0.0000000
#> 3 -0.7768362
```
