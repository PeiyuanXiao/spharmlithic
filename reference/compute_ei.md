# Elongation (E) and Isotropy (I) from the orientation tensor

Constructs the 3×3 orientation (fabric) tensor from unit direction
vectors, computes its eigenvalues \\\lambda_1 \ge \lambda_2 \ge
\lambda_3\\, and returns the shape descriptors: \$\$E = 1 - \lambda_2 /
\lambda_1, \quad I = \lambda_3 / \lambda_1.\$\$

## Usage

``` r
compute_ei(dx, dy, dz)
```

## Arguments

- dx, dy, dz:

  Numeric vectors of equal length. The X, Y, Z components of **unit**
  direction vectors, typically the `d_x`, `d_y`, `d_z` columns returned
  by
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md).
  Non-unit vectors are accepted without error but will produce
  meaningless results; normalisation is the caller's responsibility.

## Value

A one-row data frame with columns:

- E:

  Elongation index in \\\[0, 1\]\\: 1 = perfectly linear, 0 = planar /
  isotropic.

- I:

  Isotropy index in \\\[0, 1\]\\: 1 = perfectly isotropic, 0 = linear /
  planar.

- lambda1, lambda2, lambda3:

  The three eigenvalues, sorted in decreasing order.

## Details

Negative eigenvalues caused by numerical error are clamped to zero. Both
`E` and `I` are returned as `NA` when \\\lambda_1 \approx 0\\
(degenerate tensor, threshold \\10^{-10}\\).

The orientation tensor and the E/I descriptors follow Lin et al. (2024),
whose implementation adapts the fabric analysis of McPherron (2018).

The one-row data frame return type is designed to play well with
[`dplyr::group_modify()`](https://dplyr.tidyverse.org/reference/group_map.html),
allowing direct per-specimen computation.

## References

Lin, S. C., Clarkson, C., Julianto, I. M. A., Ferdianto, A., & Sutikna,
T. (2024). A new method for quantifying flake scar organisation on cores
using orientation statistics. *Journal of Archaeological Science*,
**167**, 105998.

McPherron, S. P. (2018). Additional statistical and graphical methods
for analyzing site formation processes using artifact orientations.
*PLOS ONE*, **13**(1), e0190195.

## See also

[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md),
[`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md)

## Examples

``` r
# Strongly elongated — vectors mostly aligned along X
compute_ei(c(1, 1, 0.9), c(0, 0, 0.1), c(0, 0, 0))
#>           E I lambda1     lambda2 lambda3
#> 1 0.9974723 0 0.93763 0.002370042       0

# Isotropic — vectors evenly distributed along all three axes
dx <- c(1, -1,  0,  0,  0,  0)
dy <- c(0,  0,  1, -1,  0,  0)
dz <- c(0,  0,  0,  0,  1, -1)
compute_ei(dx, dy, dz)   # I close to 1
#>   E I   lambda1   lambda2   lambda3
#> 1 0 1 0.3333333 0.3333333 0.3333333
```
