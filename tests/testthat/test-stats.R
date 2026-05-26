# ==============================================================================
# test-stats.R
# Tests for compute_SPI(), compute_spi_angle(), compute_EI()
# ==============================================================================

# ---- compute_SPI -------------------------------------------------------------

test_that("compute_SPI returns 1 for perfectly aligned vectors", {
  # All pointing in the same direction
  spi <- compute_SPI(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))
  expect_equal(spi, 1)
})

test_that("compute_SPI returns close to 0 for opposing vectors", {
  # Two vectors perfectly opposing -> resultant = 0
  spi <- compute_SPI(c(1, -1), c(0, 0), c(0, 0))
  expect_equal(spi, 0)
})

test_that("compute_SPI returns value in [0, 1]", {
  set.seed(123)
  dx <- rnorm(100)
  dy <- rnorm(100)
  dz <- rnorm(100)
  spi <- compute_SPI(dx, dy, dz)
  expect_true(spi >= 0 && spi <= 1)
})

test_that("compute_SPI length-weighted mode works", {
  dx <- c(1, 0)
  dy <- c(0, 1)
  dz <- c(0, 0)
  lengths <- c(10, 1)
  
  spi_unweighted <- compute_SPI(dx, dy, dz)
  spi_weighted   <- compute_SPI(dx, dy, dz, lengths = lengths)
  
  # Weighted should favor the longer vector (along X), giving higher SPI
  expect_true(spi_weighted > spi_unweighted)
})

test_that("compute_SPI handles single vector", {
  spi <- compute_SPI(0.5, 0.3, 0.1)
  expect_equal(spi, 1)
})


# ---- compute_spi_angle ------------------------------------------------------

test_that("compute_spi_angle returns 0 for parallel vectors", {
  angle <- compute_spi_angle(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))
  expect_equal(angle, 0)
})

test_that("compute_spi_angle returns 45 for orthogonal pair", {
  # Two orthogonal unit vectors: resultant magnitude = sqrt(2),
  
  # total length = 2, SPI = sqrt(2)/2, acos(sqrt(2)/2) = 45 degrees
  angle <- compute_spi_angle(c(1, 0), c(0, 1), c(0, 0))
  expect_equal(angle, 45, tolerance = 1e-10)
})

test_that("compute_spi_angle respects unit argument", {
  angle_deg <- compute_spi_angle(c(1, 0), c(0, 1), c(0, 0),
                                 unit = "degrees")
  angle_rad <- compute_spi_angle(c(1, 0), c(0, 1), c(0, 0),
                                 unit = "radians")
  expect_equal(angle_deg, 45, tolerance = 1e-10)
  expect_equal(angle_rad, pi / 4, tolerance = 1e-10)
})

test_that("compute_spi_angle passes lengths through", {
  # Just check it doesn't error with lengths
  angle <- compute_spi_angle(c(1, 0), c(0, 1), c(0, 0), lengths = c(5, 1))
  expect_true(is.numeric(angle) && length(angle) == 1)
})


# ---- compute_EI --------------------------------------------------------------

test_that("compute_EI returns expected columns", {
  res <- compute_EI(c(1, 1, 0.9), c(0, 0, 0.1), c(0, 0, 0))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(c("E", "I", "lambda1", "lambda2", "lambda3") %in%
                    names(res)))
})

test_that("compute_EI: strongly elongated case", {
  # All vectors along X -> E near 1, I near 0
  res <- compute_EI(c(1, 1, 1, 1), c(0, 0, 0, 0), c(0, 0, 0, 0))
  expect_equal(res$E, 1, tolerance = 1e-10)
  expect_equal(res$I, 0, tolerance = 1e-10)
})

test_that("compute_EI: isotropic case", {
  ux <- c(1, -1,  0,  0,  0,  0)
  uy <- c(0,  0,  1, -1,  0,  0)
  uz <- c(0,  0,  0,  0,  1, -1)
  res <- compute_EI(ux, uy, uz)
  
  # All eigenvalues equal -> E = 0, I = 1
  expect_equal(res$E, 0, tolerance = 1e-10)
  expect_equal(res$I, 1, tolerance = 1e-10)
})

test_that("compute_EI: eigenvalues are sorted decreasingly", {
  set.seed(7)
  res <- compute_EI(rnorm(20), rnorm(20), rnorm(20))
  expect_true(res$lambda1 >= res$lambda2)
  expect_true(res$lambda2 >= res$lambda3)
})

test_that("compute_EI: E and I in [0, 1]", {
  set.seed(99)
  res <- compute_EI(rnorm(50), rnorm(50), rnorm(50))
  expect_true(res$E >= 0 && res$E <= 1)
  expect_true(res$I >= 0 && res$I <= 1)
})