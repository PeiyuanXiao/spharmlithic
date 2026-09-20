# ==============================================================================
# test-install-python.R
# Tests for install_spharmlithic_python() and use_spharmlithic_python().
# reticulate is mocked throughout: nothing is installed and Python never starts.
# Every install call passes restart_session = FALSE so that running the tests
# inside RStudio cannot restart the session.
# ==============================================================================

# ---- Helpers -----------------------------------------------------------------

# Mock the reticulate functions install_spharmlithic_python() calls. Returns
# an environment in which each mocked call's arguments are recorded by name.
local_install_mocks <- function(conda_env_exists = FALSE,
                                venv_exists = FALSE,
                                env = parent.frame()) {
  log <- new.env()
  record <- function(name) {
    force(name)
    function(...) {
      assign(name, list(...), envir = log)
      invisible(TRUE)
    }
  }
  testthat::local_mocked_bindings(
    conda_python = function(envname, ...) {
      if (conda_env_exists) "/fake/envs/python" else stop("no such env")
    },
    virtualenv_exists  = function(envname, ...) venv_exists,
    conda_remove       = record("conda_remove"),
    virtualenv_remove  = record("virtualenv_remove"),
    conda_create       = record("conda_create"),
    conda_install      = record("conda_install"),
    virtualenv_create  = record("virtualenv_create"),
    virtualenv_install = record("virtualenv_install"),
    .package = "reticulate",
    .env     = env
  )
  log
}

# Mock the reticulate functions use_spharmlithic_python() calls
local_use_mocks <- function(conda_envs = character(0),
                            venv_exists = FALSE,
                            missing = character(0),
                            env = parent.frame()) {
  log <- new.env()
  record <- function(name) {
    force(name)
    function(...) {
      assign(name, list(...), envir = log)
      invisible(TRUE)
    }
  }
  testthat::local_mocked_bindings(
    conda_list          = function(...) data.frame(name = c("base", conda_envs)),
    virtualenv_exists   = function(envname, ...) venv_exists,
    use_condaenv        = record("use_condaenv"),
    use_virtualenv      = record("use_virtualenv"),
    use_python          = record("use_python"),
    py_config           = function(...) list(python = "/fake/bin/python"),
    py_module_available = function(module) !module %in% missing,
    .package = "reticulate",
    .env     = env
  )
  log
}


# ---- install_spharmlithic_python ---------------------------------------------

test_that("conda install builds the env from conda-forge, then pip-only packages", {
  log <- local_install_mocks()
  expect_message(
    res <- install_spharmlithic_python(method = "conda", new_env = FALSE,
                                       restart_session = FALSE),
    "installed in environment 'r-spharmlithic'"
  )
  expect_identical(res, "r-spharmlithic")

  expect_identical(log$conda_create$envname, "r-spharmlithic")
  expect_identical(log$conda_create$channel, "conda-forge")
  # The range has to reach the solver as a package spec; passing it as
  # `python_version` would render as a literal `python=>=3.11,<3.13`.
  expect_identical(log$conda_create$packages[[1]], "python>=3.11,<3.13")
  expect_null(log$conda_create$python_version)
  expect_true(all(c("numpy", "scipy", "pandas") %in%
                    log$conda_create$packages))
  expect_false("trimesh" %in% log$conda_create$packages)

  expect_identical(log$conda_install$packages, "pyshtools")
  expect_true(log$conda_install$pip)
  expect_null(log$virtualenv_create)
})

test_that("only macOS gets the pthreads OpenBLAS in the conda env", {
  # The macOS CI job covers the TRUE branch, the other jobs the FALSE one
  is_macos <- Sys.info()[["sysname"]] == "Darwin"
  log <- local_install_mocks()
  suppressMessages(
    install_spharmlithic_python(method = "conda", new_env = FALSE,
                                restart_session = FALSE)
  )
  expect_identical("libopenblas=*=*pthreads*" %in% log$conda_create$packages,
                   is_macos)
})

test_that("mesh = TRUE adds trimesh (conda) and open3d (pip)", {
  log <- local_install_mocks()
  suppressMessages(
    install_spharmlithic_python(method = "conda", mesh = TRUE,
                                new_env = FALSE, restart_session = FALSE)
  )
  expect_true("trimesh" %in% log$conda_create$packages)
  expect_identical(log$conda_install$packages, c("pyshtools", "open3d"))
})

test_that("virtualenv install puts everything through pip, no conda specs", {
  log <- local_install_mocks()
  suppressMessages(
    install_spharmlithic_python(method = "virtualenv", mesh = TRUE,
                                new_env = FALSE, restart_session = FALSE)
  )
  # virtualenv_create() takes `version`. A `python_version` argument lands
  # in `...` and is discarded silently, so assert it is not used.
  expect_identical(log$virtualenv_create$version, ">=3.11,<3.13")
  expect_null(log$virtualenv_create$python_version)
  pkgs <- log$virtualenv_install$packages
  expect_identical(pkgs, c("numpy", "scipy", "pandas", "trimesh",
                           "pyshtools", "open3d"))
  expect_false(any(grepl("=", pkgs, fixed = TRUE)))
  expect_null(log$conda_create)
})

test_that("new_env = TRUE removes an existing env before creating it", {
  log <- local_install_mocks(conda_env_exists = TRUE, venv_exists = TRUE)
  suppressMessages(
    install_spharmlithic_python(method = "auto", new_env = TRUE,
                                restart_session = FALSE)
  )
  expect_identical(log$conda_remove[[1]], "r-spharmlithic")
  expect_identical(log$virtualenv_remove[[1]], "r-spharmlithic")
  expect_false(log$virtualenv_remove$confirm)
  expect_identical(log$conda_create$envname, "r-spharmlithic")
})

test_that("nothing is removed when there is no env, or for custom env names", {
  log <- local_install_mocks(conda_env_exists = FALSE)
  suppressMessages(
    install_spharmlithic_python(method = "conda", new_env = TRUE,
                                restart_session = FALSE)
  )
  expect_null(log$conda_remove)

  # new_env defaults to FALSE for any name other than "r-spharmlithic"
  log <- local_install_mocks(conda_env_exists = TRUE)
  suppressMessages(
    install_spharmlithic_python(envname = "my-env", method = "conda",
                                restart_session = FALSE)
  )
  expect_null(log$conda_remove)
  expect_identical(log$conda_create$envname, "my-env")
})

test_that("restart_session = TRUE restarts RStudio when it can", {
  skip_if_not_installed("rstudioapi")
  local_install_mocks()
  restarted <- FALSE
  local_mocked_bindings(
    hasFun         = function(...) TRUE,
    restartSession = function(...) restarted <<- TRUE,
    .package = "rstudioapi"
  )
  suppressMessages(
    install_spharmlithic_python(method = "conda", new_env = FALSE,
                                restart_session = TRUE)
  )
  expect_true(restarted)
})


# ---- use_spharmlithic_python -------------------------------------------------

test_that("use_spharmlithic_python activates a conda env and checks packages", {
  log <- local_use_mocks(conda_envs = "r-spharmlithic")
  expect_message(
    res <- use_spharmlithic_python("r-spharmlithic"),
    "All required packages available"
  )
  expect_identical(log$use_condaenv, list("r-spharmlithic", required = TRUE))
  expect_null(log$use_virtualenv)

  expect_identical(res$envname, "r-spharmlithic")
  expect_identical(res$python, "/fake/bin/python")
  expect_identical(names(res$available),
                   c("numpy", "scipy", "pandas", "pyshtools"))
  expect_true(all(res$available))
})

test_that("use_spharmlithic_python falls back to a virtualenv, then a path", {
  log <- local_use_mocks(venv_exists = TRUE)
  suppressMessages(use_spharmlithic_python("my-venv"))
  expect_identical(log$use_virtualenv, list("my-venv", required = TRUE))

  python <- tempfile()
  file.create(python)
  log <- local_use_mocks()
  suppressMessages(use_spharmlithic_python(python))
  expect_identical(log$use_python, list(python, required = TRUE))
})

test_that("use_spharmlithic_python errors when nothing matches", {
  local_use_mocks()
  expect_error(use_spharmlithic_python("no-such-env"), "Could not find")
})

test_that("use_spharmlithic_python warns about missing packages", {
  local_use_mocks(conda_envs = "r-spharmlithic",
                  missing = c("pyshtools", "open3d"))
  expect_warning(
    res <- use_spharmlithic_python("r-spharmlithic", check_mesh = TRUE),
    "pyshtools, open3d.*mesh = TRUE"
  )
  expect_identical(names(res$available),
                   c("numpy", "scipy", "pandas", "pyshtools",
                     "trimesh", "open3d"))
  expect_identical(unname(res$available[c("pyshtools", "open3d")]),
                   c(FALSE, FALSE))
})
