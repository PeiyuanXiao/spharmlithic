# Point `spharmlithic` at an existing Python environment

Use this if you already have a conda environment with the required
Python packages (e.g. an environment created by your analysis project's
`environment.yml`, or a Docker image). Performs a lightweight check that
the required packages are importable.

## Usage

``` r
use_spharmlithic_python(envname, check_mesh = FALSE)
```

## Arguments

- envname:

  Character. Name of an existing conda environment, OR a path to a
  virtualenv root, OR a path to a Python executable.

- check_mesh:

  Logical. If `TRUE`, also verify mesh-extension packages (`trimesh`,
  `open3d`). Default `FALSE`.

## Value

Invisibly, a list with elements `envname`, `python` (path to the active
Python executable), and `available` (named logical vector of
required-package availability).

## Details

This function should be called **before** any
[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md)
or
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md)
call, and ideally before the package's first use in an R session.

The `envname` argument is resolved in the following order:

1.  **Conda environment** — if `envname` matches a name returned by
    [`reticulate::conda_list()`](https://rstudio.github.io/reticulate/reference/conda-tools.html),
    [`reticulate::use_condaenv()`](https://rstudio.github.io/reticulate/reference/use_python.html)
    is called.

2.  **Virtualenv** — if `reticulate::virtualenv_exists(envname)` is
    `TRUE`,
    [`reticulate::use_virtualenv()`](https://rstudio.github.io/reticulate/reference/use_python.html)
    is called.

3.  **Python executable path** — if `file.exists(envname)` is `TRUE`,
    [`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html)
    is called with the path directly.

If none of the above match, the function stops with an informative
error. In all cases `required = TRUE` is passed to the underlying
reticulate call, so any conflict with an already-initialised Python
session will raise an error rather than silently use a different
interpreter.

## See also

[`install_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/install_spharmlithic_python.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Reuse the conda env created by your analysis project
use_spharmlithic_python("r-spharmlithic")

# Reuse and verify mesh extension is available
use_spharmlithic_python("r-spharmlithic", check_mesh = TRUE)
} # }
```
