# ==============================================================================
# test-geometry-utils.R
# Tests for get_rot_matrix() and get_scar_length()
# ==============================================================================

# ---- get_rot_matrix ----------------------------------------------------------

test_that("get_rot_matrix rotates source to target", {
  R <- get_rot_matrix(c(0, 0, 1), c(1, 0, 0))
  result <- as.numeric(R %*% c(0, 0, 1))
  expect_equal(result, c(1, 0, 0), tolerance = 1e-10)
})

test_that("get_rot_matrix returns identity for parallel vectors", {
  R <- get_rot_matrix(c(0, 0, 1), c(0, 0, 2))
  expect_equal(R, diag(3), tolerance = 1e-10)
})

test_that("get_rot_matrix handles antiparallel vectors", {
  R <- get_rot_matrix(c(0, 0, 1), c(0, 0, -1))
  result <- as.numeric(R %*% c(0, 0, 1))
  expect_equal(result, c(0, 0, -1), tolerance = 1e-10)
})

test_that("get_rot_matrix produces an orthogonal matrix", {
  R <- get_rot_matrix(c(1, 2, 3), c(-1, 0.5, 2))
  # R^T R should be identity
  expect_equal(t(R) %*% R, diag(3), tolerance = 1e-10)
  # det(R) should be 1 (proper rotation)
  expect_equal(det(R), 1, tolerance = 1e-10)
})

test_that("get_rot_matrix normalises inputs", {
  # Non-unit vectors should work fine
  R1 <- get_rot_matrix(c(0, 0, 5), c(10, 0, 0))
  R2 <- get_rot_matrix(c(0, 0, 1), c(1, 0, 0))
  expect_equal(R1, R2, tolerance = 1e-10)
})

test_that("get_rot_matrix handles various axis pairs", {
  pairs <- list(
    list(c(1, 0, 0), c(0, 1, 0)),
    list(c(0, 1, 0), c(0, 0, 1)),
    list(c(1, 1, 0), c(0, 1, 1)),
    list(c(1, 1, 1), c(-1, -1, -1))
  )
  for (p in pairs) {
    R <- get_rot_matrix(p[[1]], p[[2]])
    result <- as.numeric(R %*% (p[[1]] / sqrt(sum(p[[1]]^2))))
    target <- p[[2]] / sqrt(sum(p[[2]]^2))
    expect_equal(result, target, tolerance = 1e-10)
  }
})


# ---- get_scar_length ---------------------------------------------------------

test_that("get_scar_length computes from coordinates", {
  df <- data.frame(
    Start_X = c(0, 0),
    Start_Y = c(0, 0),
    Start_Z = c(0, 0),
    End_X   = c(3, 0),
    End_Y   = c(4, 5),
    End_Z   = c(0, 12)
  )
  lens <- get_scar_length(df)
  expect_equal(lens, c(5, 13))
})

test_that("get_scar_length uses Length column if present", {
  df <- data.frame(
    Length  = c(10, 20),
    Start_X = c(0, 0), Start_Y = c(0, 0), Start_Z = c(0, 0),
    End_X   = c(3, 0), End_Y   = c(4, 5), End_Z   = c(0, 12)
  )
  lens <- get_scar_length(df)
  expect_equal(lens, c(10, 20))
})

test_that("get_scar_length handles zero-length scars", {
  df <- data.frame(
    Start_X = c(1, 5),
    Start_Y = c(2, 6),
    Start_Z = c(3, 7),
    End_X   = c(1, 5),
    End_Y   = c(2, 6),
    End_Z   = c(3, 7)
  )
  lens <- get_scar_length(df)
  expect_equal(lens, c(0, 0))
})