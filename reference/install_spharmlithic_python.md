# Install the Python backend for `spharmlithic`

Creates (or recreates) a conda environment containing the Python
packages required for spherical harmonic analysis. By default installs
the **core** environment, sufficient for the direction-vector pipeline
([`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md)).
Pass `mesh = TRUE` to additionally install `trimesh` and `open3d` for
the STL pipeline
([`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md)).

## Usage

``` r
install_spharmlithic_python(
  envname = "r-spharmlithic",
  method = c("auto", "conda", "virtualenv"),
  mesh = FALSE,
  python_version = ">=3.11,<3.13",
  new_env = identical(envname, "r-spharmlithic"),
  restart_session = TRUE
)
```

## Arguments

- envname:

  Character. Conda environment name. Default `"r-spharmlithic"`.

- method:

  Character. One of `"auto"`, `"conda"`, `"virtualenv"`. Default
  `"auto"`, which prefers conda when available (recommended, because
  `pyshtools` is easier to install via conda-forge).

- mesh:

  Logical. If `TRUE`, install the mesh-processing extension (`trimesh`,
  `open3d`) needed for
  [`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md).
  Default `FALSE` (core only).

- python_version:

  Character. Python version constraint. Default `">=3.11,<3.13"`. The
  lower bound is the oldest version the current `pyshtools`, `open3d`,
  `numpy` and `scipy` releases still publish wheels for; the upper bound
  is required because `pyshtools` has no Python 3.13 wheels. Both bounds
  now reach the installer: conda and virtualenv each take the newest
  version inside the range (3.12 at the time of writing).

- new_env:

  Logical. If `TRUE` (default when `envname` is the package default),
  remove any existing environment with the same name first.

- restart_session:

  Logical. Restart R session after install (only affects RStudio).

## Value

Invisibly, the environment name.

## Details

On macOS the conda environment uses the pthreads build of OpenBLAS, so
that it does not load a second copy of the OpenMP runtime next to the
one that comes with R (two copies abort the R session with
`OMP: Error #15`). Results are identical to the default build. If you
created the environment with an earlier version of spharmlithic, run
`install_spharmlithic_python()` again to recreate it. The mesh extension
cannot be used in a native macOS installation; see
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md).

## See also

[`use_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/use_spharmlithic_python.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Core install (direction-vector pipeline only)
install_spharmlithic_python()

# Full install (also enables STL pipeline)
install_spharmlithic_python(mesh = TRUE)
} # }
```
