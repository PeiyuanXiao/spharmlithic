# Convert SPI to a scar-pattern angle (Clarkson et al. 2006)

Converts the raw Scar Pattern Index of one core into an angle, \\\theta
= \arccos(\mathrm{SPI})\\, following Clarkson et al. (2006). SPI = 1
maps to 0 degrees (parallel scars) and SPI = 0 to 90 degrees.

## Usage

``` r
compute_spi_angle(dx, dy, dz, lengths = NULL, unit = c("degrees", "radians"))
```

## Arguments

- dx, dy, dz:

  Numeric vectors of equal length: the X, Y and Z components of the scar
  direction vectors of **one** core, for example `End_X - Start_X` etc.,
  or the `d_x`, `d_y`, `d_z` columns returned by
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md).
  The vectors are normalised to unit length inside the function; rows of
  zero length or with missing values are left out.

- lengths:

  Optional numeric vector of scar lengths, one per vector. If given,
  each unit direction is weighted by its length, reproducing the
  length-weighted SPI of Clarkson et al. (2006). The rarefied value and
  the random interval are defined for unweighted SPI only and are `NA`
  when `lengths` is given.

- unit:

  Either `"degrees"` (default) or `"radians"`.

## Value

A single numeric value in \\\[0, 90\]\\ degrees (or \\\[0, \pi/2\]\\
radians).

## Details

Clarkson et al. (2006) read this angle as the mean angle between any two
scars. That reading is only approximate: SPI includes each scar paired
with itself. Over the distinct pairs of n unit vectors the mean cosine
is \\(n \\ \mathrm{SPI}^2 - 1) / (n - 1)\\, which is not
\\\mathrm{SPI}\\; it is the degree-1 direction feature of
[`spharm_spatiodirection()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_spatiodirection.md)
(`dir_1`, up to the factor \\\sqrt{3/2}\\). The angle is computed from
the raw SPI and, like it, shifts with the scar count.

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
