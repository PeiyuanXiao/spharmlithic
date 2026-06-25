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


# ---- replace_zeros --------------------------------------------------------------

test_that("rows without zeros are unchanged", {
  x <- matrix(
    c(0.2, 0.3, 0.5,
      0.1, 0.4, 0.5),
    nrow = 2,
    byrow = TRUE
  )
  
  expect_equal(replace_zeros(x), x)
})

test_that("zeros are replaced and rows still sum to one", {
  x <- matrix(
    c(0.5, 0.5, 0,
      0.2, 0.3, 0.5),
    nrow = 2,
    byrow = TRUE
  )
  
  y <- replace_zeros(x)
  
  expect_true(all(y > 0))
  expect_equal(rowSums(y), c(1, 1))
})

test_that("fixed delta is applied correctly", {
  x <- matrix(c(0.5, 0.5, 0), nrow = 1)
  
  y <- replace_zeros(x, delta = 0.01)
  
  expect_equal(y[1, 3], 0.01)
  expect_equal(sum(y), 1)
})

test_that("all-zero rows generate a warning and are unchanged", {
  x <- matrix(
    c(0, 0, 0,
      0.5, 0.5, 0),
    nrow = 2,
    byrow = TRUE
  )
  
  expect_warning(
    y <- replace_zeros(x),
    "contains only zeros"
  )
  
  expect_equal(y[1, ], c(0, 0, 0))
})

test_that("non-numeric input throws an error", {
  x <- matrix(c("a", "b", "c"), nrow = 1)
  
  expect_error(
    replace_zeros(x),
    "`x` must be numeric"
  )
})

test_that("rows remain compositional after replacement", {
  x <- matrix(
    c(0.7, 0.3, 0,
      0.4, 0, 0.6,
      0, 0.2, 0.8),
    nrow = 3,
    byrow = TRUE
  )
  
  y <- replace_zeros(x)
  
  expect_equal(
    rowSums(y),
    rep(1, nrow(y)),
    tolerance = 1e-12
  )
  
  expect_true(all(y[y != 0] > 0))
})

test_that("delta that is too large throws an error", {
  x <- matrix(c(0.5, 0.5, 0), nrow = 1)
  
  expect_error(
    replace_zeros(x, delta = 1),
    "Replacement value is too large"
  )
})

test_that("fraction controls replacement value", {
  x <- matrix(c(0.5, 0.5, 0), nrow = 1)
  y <- replace_zeros(x, fraction = 0.5)
  expect_equal(y[1, 3], 0.25)
})

# ---- make_ilr ------------------------------------------------------------------------

test_that("returns a data frame", {
  x <- data.frame(
    a = c(0.5, 0.4, 0.6),
    b = c(0.3, 0.4, 0.2),
    c = c(0.2, 0.2, 0.2)
  )
  
  expect_s3_class(make_ilr(x), "data.frame")
})

test_that("removes constant columns", {
  x <- data.frame(
    a = c(0.5, 0.4, 0.6),
    b = c(0.3, 0.4, 0.2),
    c = c(0.2, 0.2, 0.2)
  )
  
  result <- make_ilr(x)
  
  expect_equal(ncol(result), 1)
})

test_that("preserves number of rows", {
  x <- data.frame(
    a = c(0.5, 0.4, 0.6),
    b = c(0.3, 0.4, 0.2),
    c = c(0.2, 0.2, 0.2)
  )
  
  expect_equal(
    nrow(make_ilr(x)),
    nrow(x)
  )
})

test_that("handles zeros via replacement", {
  x <- data.frame(
    a = c(0.5, 0.4, 0.6),
    b = c(0.5, 0.6, 0.4),
    c = c(0, 0, 0)
  )
  
  expect_no_error(make_ilr(x))
})

test_that("errors when fewer than two variable columns remain", {
  x <- data.frame(
    a = c(1, 1, 1),
    b = c(2, 2, 2),
    c = c(3, 3, 3)
  )
  
  expect_error(
    make_ilr(x),
    "At least two non-constant columns"
  )
})

test_that("errors for non-numeric input", {
  x <- data.frame(
    a = c("a", "b", "c"),
    b = c("d", "e", "f")
  )
  
  expect_error(
    make_ilr(x),
    "`x` must be numeric"
  )
})

test_that("assigns ILR column names", {
  x <- data.frame(
    a = c(0.5, 0.4, 0.6),
    b = c(0.3, 0.4, 0.2),
    c = c(0.2, 0.2, 0.2)
  )
  
  result <- make_ilr(x)
  
  expect_equal(
    colnames(result),
    "ilr_1"
  )
})

# ---- degree_diagnostics ----------------------------------------------------

test_that("returns expected columns", {
  
  x <- data.frame(
    power_l1 = c(0.5, 0.6),
    power_l2 = c(0.3, 0.2),
    power_l3 = c(0.2, 0.2)
  )
  
  result <- degree_diagnostics(
    x,
    descriptor = "test",
    max_degree = 3
  )
  
  expect_equal(
    names(result),
    c(
      "descriptor",
      "degree",
      "mean_power",
      "cv_pct",
      "cumul_pct"
    )
  )
})

test_that("returns one row per degree", {
  
  x <- data.frame(
    power_l1 = c(0.5, 0.6),
    power_l2 = c(0.3, 0.2),
    power_l3 = c(0.2, 0.2)
  )
  
  result <- degree_diagnostics(
    x,
    descriptor = "test",
    max_degree = 3
  )
  
  expect_equal(nrow(result), 3)
})

test_that("cumulative power ends at 100 percent", {
  
  x <- data.frame(
    power_l1 = c(0.5, 0.6),
    power_l2 = c(0.3, 0.2),
    power_l3 = c(0.2, 0.2)
  )
  
  result <- degree_diagnostics(
    x,
    descriptor = "test",
    max_degree = 3
  )
  
  expect_equal(
    tail(result$cumul_pct, 1),
    100,
    tolerance = 1e-10
  )
})

test_that("descriptor is propagated", {
  
  x <- data.frame(
    power_l1 = c(0.5, 0.6),
    power_l2 = c(0.3, 0.2),
    power_l3 = c(0.2, 0.2)
  )
  
  result <- degree_diagnostics(
    x,
    descriptor = "power",
    max_degree = 3
  )
  
  expect_true(all(result$descriptor == "power"))
})

test_that("zero mean power yields NA cv", {
  
  x <- data.frame(
    power_l1 = c(0, 0),
    power_l2 = c(0.5, 0.5),
    power_l3 = c(0.5, 0.5)
  )
  
  result <- degree_diagnostics(
    x,
    descriptor = "test",
    max_degree = 3
  )
  
  expect_true(is.na(result$cv_pct[1]))
})

test_that("missing power columns trigger error", {
  
  x <- data.frame(
    power_l1 = c(0.5, 0.6),
    power_l2 = c(0.3, 0.2)
  )
  
  expect_error(
    degree_diagnostics(
      x,
      descriptor = "test",
      max_degree = 3
    ),
    "Missing required columns"
  )
})

# ---- pairwise_permanova ---------------------------------------------

test_that("returns expected columns", {
  
  set.seed(1)
  
  x <- matrix(rnorm(45), ncol = 3)
  grp <- rep(c("A", "B", "C"), each = 5)
  
  result <- pairwise_permanova(
    x,
    grp,
    permutations = 9
  )
  
  expect_equal(
    names(result),
    c(
      "group1",
      "group2",
      "pair",
      "n1",
      "n2",
      "R2",
      "F",
      "p",
      "p_adj"
    )
  )
})

test_that("returns correct number of pairwise comparisons", {
  
  x <- matrix(rnorm(45), ncol = 3)
  grp <- rep(c("A", "B", "C"), each = 5)
  
  result <- pairwise_permanova(
    x,
    grp,
    permutations = 9
  )
  
  expect_equal(nrow(result), 3)
})

test_that("adjusted p values are present", {
  
  x <- matrix(rnorm(45), ncol = 3)
  grp <- rep(c("A", "B", "C"), each = 5)
  
  result <- pairwise_permanova(
    x,
    grp,
    permutations = 9
  )
  
  expect_false(any(is.na(result$p_adj)))
})

test_that("errors when dimensions do not match", {
  
  x <- matrix(rnorm(20), ncol = 2)
  grp <- rep(c("A", "B"), each = 6)
  
  expect_error(
    pairwise_permanova(x, grp),
    "Number of rows"
  )
})

test_that("errors when only one group is present", {
  
  x <- matrix(rnorm(20), ncol = 2)
  grp <- rep("A", 10)
  
  expect_error(
    pairwise_permanova(x, grp),
    "At least two groups"
  )
})

test_that("returns one comparison for two groups", {
  
  x <- matrix(rnorm(40), ncol = 2)
  grp <- rep(c("A", "B"), each = 10)
  
  result <- pairwise_permanova(
    x,
    grp,
    permutations = 9
  )
  
  expect_equal(nrow(result), 1)
})

test_that("small groups return NA statistics", {
  
  x <- matrix(rnorm(20), ncol = 2)
  grp <- c("A", "A", "B", rep("C", 7))
  
  result <- pairwise_permanova(
    x,
    grp,
    permutations = 9
  )
  
  expect_true(any(is.na(result$R2)))
})