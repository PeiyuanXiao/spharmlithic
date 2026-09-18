# Scar Pattern Index (Clarkson et al. 2006)

Computes the ratio of the resultant vector magnitude to the total scar
length, on a scale from 0 (random orientation) to 1 (perfect alignment).

## Usage

``` r
compute_spi(dx, dy, dz, lengths = NULL)
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

## Value

A single numeric value in \\\[0, 1\]\\. Values close to 1 indicate
strong preferred orientation; values close to 0 indicate isotropic /
random patterning.

## Details

**Default behaviour:** unweighted (Bretzke & Conard 2012) — all scars
contribute equally regardless of length. Pass `lengths` to get the
length-weighted SPI (Clarkson et al. 2006).

Two variants of SPI are supported:

- **Length-weighted** (Clarkson et al. 2006): pass raw scar displacement
  vectors so that longer scars contribute proportionally more to the
  resultant, or equivalently pass unit vectors together with `lengths`.
  This reproduces the original definition.

- **Unweighted** (Bretzke & Conard 2012): leave `lengths = NULL`. Every
  scar contributes equally regardless of length — useful when scar
  lengths are unreliable or when only direction matters.

To reproduce the length-weighted definition from aligned batch output:


      lens <- sqrt((aligned$e_x - aligned$s_x)^2 +
                   (aligned$e_y - aligned$s_y)^2 +
                   (aligned$e_z - aligned$s_z)^2)
      compute_spi(aligned$d_x, aligned$d_y, aligned$d_z, lengths = lens)

## References

Clarkson, C., Vinicius, L., & Lahr, M. M. (2006). Quantifying flake scar
patterning on cores using 3D recording techniques. *Journal of
Archaeological Science*, **33**(1), 132–142.

Bretzke, K., & Conard, N. J. (2012). Evaluating morphological
variability in lithic assemblages using 3D models of stone artifacts.
*Journal of Archaeological Science*, **39**(12), 3741–3749.

## See also

[`compute_spi_angle()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi_angle.md),
[`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md),
[`get_scar_length()`](https://peiyuanxiao.github.io/spharmlithic/reference/get_scar_length.md)

## Examples

``` r
# Perfectly aligned vectors along X
compute_spi(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))   # 1
#> [1] 1

# Random-like distribution (close to 0)
set.seed(1)
compute_spi(rnorm(50), rnorm(50), rnorm(50))
#> [1] 0.1509766

if (FALSE) { # \dontrun{
# Unweighted (default) — every scar counts equally
compute_spi(aligned$d_x, aligned$d_y, aligned$d_z)

# Length-weighted — longer scars contribute more
lens <- get_scar_length(aligned)
compute_spi(aligned$d_x, aligned$d_y, aligned$d_z, lengths = lens)
} # }
```
