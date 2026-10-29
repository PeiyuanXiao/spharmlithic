# Random interval of SPI, ER or IR at a scar count

Returns the central 95% interval of an index for a core whose scar
directions are uniform on the sphere, at each requested scar count. A
raw reading outside this interval departs from what random directions
produce at that scar count; a reading inside it cannot be told apart
from random directions at that scar count.

## Usage

``` r
random_interval(parameter = c("SPI", "ER", "IR"), n)
```

## Arguments

- parameter:

  Character vector: `"SPI"`, `"ER"` and/or `"IR"` (case-insensitive).

- n:

  Scar counts (whole numbers). Counts below 3 give `NA`.

## Value

A tibble with one row per combination of `parameter` and `n`:
`parameter`, `n`, `lower` (2.5th percentile) and `upper` (97.5th
percentile).

## Details

The interval narrows as the scar count grows, so a raw value can only be
judged against the interval of its own scar count.
[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md),
[`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md)
and
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md)
look the interval up for you.

SPI is computed from its exact distribution (the resultant length of n
unit vectors uniform on the sphere), ER and IR from 5000 simulated
random cores per scar count. Values for n = 3 to 200 are tabulated; SPI
above 200 is computed on the fly; ER and IR above 200 are `NA`.

Which side of the interval indicates organisation depends on the index:
SPI and ER above the upper bound, and IR below the lower bound, indicate
scar directions more concentrated than random. SPI below the lower bound
is typical of opposed (bidirectional) flaking.

## References

Rayleigh, Lord (1919). On the problem of random vibrations, and of
random flights in one, two, or three dimensions. *Philosophical
Magazine*, **37**(220), 321–347.

## See also

[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md),
[`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md),
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md),
[`plot_random_interval()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_random_interval.md)

## Examples

``` r
random_interval("SPI", n = c(6, 20))
#> # A tibble: 2 × 4
#>   parameter     n  lower upper
#>   <chr>     <dbl>  <dbl> <dbl>
#> 1 SPI           6 0.114  0.692
#> 2 SPI          20 0.0607 0.390
random_interval(c("ER", "IR"), n = 10)
#> # A tibble: 2 × 4
#>   parameter     n  lower upper
#>   <chr>     <dbl>  <dbl> <dbl>
#> 1 ER           10 0.0800 0.653
#> 2 IR           10 0.120  0.632
```
