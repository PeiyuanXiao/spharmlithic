# ==============================================================================
# test-spharm-meshes.R
# Tests for spharm_from_meshes() — requires Python mesh extension
# ==============================================================================

test_that("spharm_from_meshes errors on non-existent directory", {
  # This test does NOT need Python
  expect_error(
    spharm_from_meshes("/nonexistent/path/to/stls"),
    "does not exist"
  )
})

test_that("spharm_from_meshes returns expected structure", {
  skip_if_no_python_mesh()
  
  stl_dir <- system.file("extdata", "meshes", package = "spharmlithic")
  skip_if(stl_dir == "" || !dir.exists(stl_dir),
          "Example STL directory not found")
  
  stl_files <- list.files(stl_dir, pattern = "\\.stl$", ignore.case = TRUE)
  skip_if(length(stl_files) == 0, "No STL files found in example meshes")
  
  res <- spharm_from_meshes(stl_dir, lmax = 4, target_faces = 5000,
                            grid_size = 32, verbose = FALSE)
  
  expect_type(res, "list")
  expect_true(length(res) > 0)
  
  for (id in names(res)) {
    s <- res[[id]]
    expect_named(s, c("coefficients", "power_spectrum"))
    expect_equal(dim(s$coefficients), c(2L, 5L, 5L))  # (2, lmax+1, lmax+1)
    expect_length(s$power_spectrum, 5L)
    expect_true(all(s$power_spectrum >= 0))
  }
})