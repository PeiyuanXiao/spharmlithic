test_that("spharm_from_directions returns expected structure", {
  skip_if_no_python_core()

  d   <- make_test_directions(n_specimens = 2, n_scars = 25)
  res <- spharm_from_directions(d, lmax = 8, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("specimen_001", "specimen_002"))

  for (id in names(res)) {
    s <- res[[id]]
    expect_named(s, c("coefficients", "power_spectrum"))
    expect_equal(dim(s$coefficients), c(2L, 9L, 9L))   # (2, lmax+1, lmax+1)
    expect_length(s$power_spectrum, 9L)                # lmax+1
    expect_true(all(s$power_spectrum >= 0))
  }
})

test_that("spharm_from_directions errors on missing columns", {
  skip_if_no_python_core()
  d <- data.frame(ID = "x", d_x = 1, d_y = 0)  # missing d_z
  expect_error(spharm_from_directions(d), "missing required columns")
})
