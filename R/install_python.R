# ==============================================================================
# install_python.R
# Python backend installation and environment-pointing helpers.
# ==============================================================================

#' Install the Python backend for `spharmlithic`
#'
#' Creates (or recreates) a conda environment containing the Python
#' packages required for spherical harmonic analysis. By default installs
#' the **core** environment, sufficient for the direction-vector pipeline
#' ([spharm_from_directions()]). Pass `mesh = TRUE` to additionally
#' install `trimesh` and `open3d` for the STL pipeline
#' ([spharm_from_meshes()]).
#'
#' @param envname Character. Conda environment name. Default
#'   `"r-spharmlithic"`.
#' @param method Character. One of `"auto"`, `"conda"`, `"virtualenv"`.
#'   Default `"auto"`, which prefers conda when available (recommended,
#'   because `pyshtools` is easier to install via conda-forge).
#' @param mesh Logical. If `TRUE`, install the mesh-processing extension
#'   (`trimesh`, `open3d`) needed for [spharm_from_meshes()]. Default
#'   `FALSE` (core only).
#' @param python_version Character. Python version constraint passed to
#'   reticulate. Default `">=3.9,<3.13"`.
#' @param new_env Logical. If `TRUE` (default when `envname` is the
#'   package default), remove any existing environment with the same name
#'   first.
#' @param restart_session Logical. Restart R session after install
#'   (only affects RStudio).
#'
#' @return Invisibly, the environment name.
#'
#' @examples
#' \dontrun{
#' # Core install (direction-vector pipeline only)
#' install_spharmlithic_python()
#'
#' # Full install (also enables STL pipeline)
#' install_spharmlithic_python(mesh = TRUE)
#' }
#'
#' @seealso [use_spharmlithic_python()]
#' @export
install_spharmlithic_python <- function(
    envname         = "r-spharmlithic",
    method          = c("auto", "conda", "virtualenv"),
    mesh            = FALSE,
    python_version  = ">=3.9,<3.13",
    new_env         = identical(envname, "r-spharmlithic"),
    restart_session = TRUE) {
  
  method <- match.arg(method)
  
  # ---- Tear down any pre-existing env with this name --------------------
  if (isTRUE(new_env)) {
    if (method %in% c("auto", "conda") &&
        length(tryCatch(reticulate::conda_python(envname),
                        error = function(e) NULL))) {
      reticulate::conda_remove(envname)
    }
    if (method %in% c("auto", "virtualenv") &&
        reticulate::virtualenv_exists(envname)) {
      reticulate::virtualenv_remove(envname, confirm = FALSE)
    }
  }
  
  # ---- Build package list -----------------------------------------------
  # Conda packages (installed via conda-forge for compiled deps)
  conda_pkgs <- c("numpy", "scipy", "pandas")
  if (mesh) conda_pkgs <- c(conda_pkgs, "trimesh")
  
  # Pip packages (pyshtools and open3d distribute as wheels)
  pip_pkgs <- c("pyshtools")
  if (mesh) pip_pkgs <- c(pip_pkgs, "open3d")
  
  # ---- Install ----------------------------------------------------------
  if (method %in% c("auto", "conda")) {
    # Create env with conda-forge channel + Python pinned version
    reticulate::conda_create(
      envname        = envname,
      packages       = conda_pkgs,
      python_version = sub("^>=", "", strsplit(python_version, ",")[[1]][1]),
      channel        = "conda-forge"
    )
    # Install pip-only packages on top
    reticulate::conda_install(
      envname  = envname,
      packages = pip_pkgs,
      pip      = TRUE
    )
  } else {
    # virtualenv path: install everything via pip
    reticulate::virtualenv_create(
      envname        = envname,
      python_version = python_version
    )
    reticulate::virtualenv_install(
      envname  = envname,
      packages = c(conda_pkgs, pip_pkgs)
    )
  }
  
  message("\nspharmlithic: Python backend installed in environment '",
          envname, "'.\n",
          "  Restart R, then call library(spharmlithic) to use it.")
  
  if (isTRUE(restart_session) &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::hasFun("restartSession")) {
    rstudioapi::restartSession()
  }
  
  invisible(envname)
}


#' Point `spharmlithic` at an existing Python environment
#'
#' Use this if you already have a conda environment with the required
#' Python packages (e.g. an environment created by your analysis project's
#' `environment.yml`, or a Docker image). Performs a lightweight check
#' that the required packages are importable.
#'
#' This function should be called **before** any [spharm_from_directions()]
#' or [spharm_from_meshes()] call, and ideally before the package's first
#' use in an R session.
#'
#' @param envname Character. Name of an existing conda environment, OR a
#'   path to a virtualenv root, OR a path to a Python executable.
#' @param check_mesh Logical. If `TRUE`, also verify mesh-extension
#'   packages (`trimesh`, `open3d`). Default `FALSE`.
#'
#' @return Invisibly, a list with elements `envname`, `python` (path to
#'   the active Python executable), and `available` (named logical
#'   vector of required-package availability).
#'
#' @details
#' The `envname` argument is resolved in the following order:
#' \enumerate{
#'   \item **Conda environment** — if `envname` matches a name returned
#'     by `reticulate::conda_list()`, `reticulate::use_condaenv()` is
#'     called.
#'   \item **Virtualenv** — if `reticulate::virtualenv_exists(envname)`
#'     is `TRUE`, `reticulate::use_virtualenv()` is called.
#'   \item **Python executable path** — if `file.exists(envname)` is
#'     `TRUE`, `reticulate::use_python()` is called with the path
#'     directly.
#' }
#' If none of the above match, the function stops with an informative
#' error. In all cases `required = TRUE` is passed to the underlying
#' reticulate call, so any conflict with an already-initialised Python
#' session will raise an error rather than silently use a different
#' interpreter.
#'
#' @examples
#' \dontrun{
#' # Reuse the conda env created by your analysis project
#' use_spharmlithic_python("r-spharmlithic")
#'
#' # Reuse and verify mesh extension is available
#' use_spharmlithic_python("r-spharmlithic", check_mesh = TRUE)
#' }
#'
#' @seealso [install_spharmlithic_python()]
#' @export
use_spharmlithic_python <- function(envname, check_mesh = FALSE) {
  
  # Try as conda env first, then virtualenv, then literal path.
  if (envname %in% reticulate::conda_list()$name) {
    reticulate::use_condaenv(envname, required = TRUE)
  } else if (reticulate::virtualenv_exists(envname)) {
    reticulate::use_virtualenv(envname, required = TRUE)
  } else if (file.exists(envname)) {
    reticulate::use_python(envname, required = TRUE)
  } else {
    stop("Could not find a conda env, virtualenv, or Python executable named '",
         envname, "'.", call. = FALSE)
  }
  
  # Force initialization so we can introspect what's loaded.
  reticulate::py_config()
  
  required_core <- c("numpy", "scipy", "pandas", "pyshtools")
  required_mesh <- c("trimesh", "open3d")
  to_check      <- if (check_mesh) c(required_core, required_mesh) else required_core
  
  available <- vapply(to_check,
                      reticulate::py_module_available,
                      logical(1))
  
  missing <- names(available)[!available]
  if (length(missing)) {
    warning("The following Python packages are NOT available in '", envname,
            "': ", paste(missing, collapse = ", "),
            "\n  Install them, or run install_spharmlithic_python(",
            if (check_mesh) "mesh = TRUE" else "",
            ").",
            call. = FALSE)
  } else {
    message("spharmlithic: using Python from '", envname,
            "'. All required packages available.")
  }
  
  invisible(list(
    envname   = envname,
    python    = reticulate::py_config()$python,
    available = available
  ))
}