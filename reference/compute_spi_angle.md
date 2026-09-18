# Convert SPI to a scar-pattern angle (Clarkson et al. 2006)

Converts the Scar Pattern Index into the expected pairwise angle between
two randomly-selected scars, following Clarkson et al.'s (2006)
interpretation: \\\theta = \arccos(\mathrm{SPI})\\.

## Usage

``` r
compute_spi_angle(dx, dy, dz, lengths = NULL, unit = c("degrees", "radians"))
```

## Arguments

- dx, dy, dz:

  Numeric vectors of equal length. The X, Y, Z components of the scar
  direction vectors. These are typically unit direction vectors returned
  by
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)
  (columns `d_x`, `d_y`, `d_z`).

- lengths:

  Optional numeric vector of scar lengths. If provided, each direction
  vector is scaled by its length before summation, reproducing the
  original length-weighted SPI of Clarkson et al. (2006). If `NULL` (the
  default), all scars contribute equally regardless of length, following
  the unit-vector variant of Bretzke & Conard (2012).

- unit:

  Either `"degrees"` (default) or `"radians"`.

## Value

A single numeric value in \\\[0, 90\]\\ degrees (or \\\[0, \pi/2\]\\
radians).

## Details

SPI = 1 maps to 0 degrees (parallel scars); SPI = 0 maps to 90 degrees
(uniformly random pairwise angles).

Before applying [`acos()`](https://rdrr.io/r/base/Trig.html), the SPI
value is clamped to \\\[-1, 1\]\\ to guard against floating-point
overshoot that would otherwise produce `NaN`.

## References

Clarkson, C., Vinicius, L., & Lahr, M. M. (2006). Quantifying flake scar
patterning on cores using 3D recording techniques. *Journal of
Archaeological Science*, **33**(1), 132–142.

## See also

[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md)

## Examples

``` r
compute_spi_angle(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))   # 0 (parallel)
#> [1] 0
compute_spi_angle(c(1, 0),    c(0, 1),    c(0, 0))       # 45 (orthogonal pair)
#> [1] 45
```
