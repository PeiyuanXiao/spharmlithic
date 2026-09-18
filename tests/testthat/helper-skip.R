# Skip helpers --------------------------------------------------------------

skip_if_no_python_core <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(
    reticulate::py_module_available("pyshtools"),
    "Python backend (pyshtools) not available"
  )
}

skip_if_no_python_mesh <- function() {
  # Must run before anything imports open3d. On macOS the open3d wheel ships
  # its own OpenMP runtime (libomp), which clashes with the one CRAN R ships
  # (initialised by data.table), and the second copy aborts the whole R
  # process with "OMP: Error #15". The mesh pipeline is not supported
  # natively on macOS (use Docker).
  testthat::skip_on_os("mac")
  skip_if_no_python_core()
  testthat::skip_if_not(
    reticulate::py_module_available("trimesh") &&
      reticulate::py_module_available("open3d"),
    "Mesh extension (trimesh + open3d) not available"
  )
}

# Tiny synthetic dataset ----------------------------------------------------

make_test_directions <- function(n_specimens = 3, n_scars = 30, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_specimens), function(i) {
    # Random unit vectors clustered near (1, 0, 0) plus noise
    base   <- matrix(rnorm(3 * n_scars, sd = 0.3), n_scars, 3)
    base[, 1] <- base[, 1] + 1
    norms  <- sqrt(rowSums(base^2))
    u      <- base / norms
    data.frame(
      ID  = paste0("specimen_", sprintf("%03d", i)),
      d_x = u[, 1],
      d_y = u[, 2],
      d_z = u[, 3],
      stringsAsFactors = FALSE
    )
  }))
}
