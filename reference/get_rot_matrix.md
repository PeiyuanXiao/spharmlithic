# Compute a rotation matrix that maps vector `a` onto vector `b`

Uses the Rodrigues rotation formula to construct the 3×3 rotation matrix
that rotates unit vector `a` into unit vector `b`.

## Usage

``` r
get_rot_matrix(a, b)
```

## Arguments

- a:

  Numeric vector of length 3. The source direction (need not be
  normalised; it will be normalised internally).

- b:

  Numeric vector of length 3. The target direction (need not be
  normalised).

## Value

A 3×3 numeric rotation matrix.

## Details

Two degenerate cases are handled (threshold: \\10^{-10}\\):

- **Antiparallel** (\\\cos\theta \< -1 + 10^{-10}\\): a 180° rotation
  about an arbitrary perpendicular axis is returned.

- **Parallel** (\\\cos\theta \> 1 - 10^{-10}\\): the identity matrix is
  returned.

This function is a low-level geometric utility; in typical usage it is
called internally by `align_scar()` and `align_morph()` rather than
directly.

## See also

[`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md),
[`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md)

## Examples

``` r
# Rotate the Z-axis onto the X-axis
R <- get_rot_matrix(c(0, 0, 1), c(1, 0, 0))
round(R %*% c(0, 0, 1), 10)   # should equal c(1, 0, 0)
#>      [,1]
#> [1,]    1
#> [2,]    0
#> [3,]    0

# Parallel vectors — identity matrix
R_id <- get_rot_matrix(c(0, 0, 1), c(0, 0, 2))
all.equal(R_id, diag(3))       # TRUE
#> [1] TRUE

# Antiparallel vectors — 180-degree rotation
R_flip <- get_rot_matrix(c(0, 0, 1), c(0, 0, -1))
round(R_flip %*% c(0, 0, 1), 10)   # should equal c(0, 0, -1)
#>      [,1]
#> [1,]    0
#> [2,]    0
#> [3,]   -1
```
