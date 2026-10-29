# ==============================================================================
# test-stats.R
# Tests for compute_spi(), compute_spi_angle(), compute_ei(), and the
# compositional and group helpers of the power-spectrum workflow
# ==============================================================================

# ---- compute_spi -------------------------------------------------------------

test_that("compute_spi returns a one-row tibble with raw, rarefied and interval", {
  set.seed(1)
  res <- compute_spi(rnorm(12), rnorm(12), rnorm(12))
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1)
  expect_named(res, c("n_scar", "SPI_raw", "SPI_rarefied", "SPI_lower",
                      "SPI_upper", "SPI_vs_random"))
  expect_equal(res$n_scar, 12L)
  expect_true(res$SPI_vs_random %in% c("within", "outside"))
})

test_that("compute_spi returns 1 for perfectly aligned vectors", {
  res <- compute_spi(c(1, 1, 1), c(0, 0, 0), c(0, 0, 0))
  expect_equal(res$SPI_raw, 1)
  expect_equal(res$SPI_rarefied, 1)
  expect_equal(res$SPI_vs_random, "outside")
})

test_that("compute_spi returns 0 for opposing vectors", {
  res <- compute_spi(c(1, -1, 1, -1), c(0, 0, 0, 0), c(0, 0, 0, 0))
  expect_equal(res$SPI_raw, 0)
})

test_that("compute_spi normalises its input: scaling a vector changes nothing", {
  set.seed(2)
  d <- matrix(rnorm(30), 10, 3)
  a <- compute_spi(d[, 1], d[, 2], d[, 3])
  b <- compute_spi(d[, 1] * 1:10, d[, 2] * 1:10, d[, 3] * 1:10)
  expect_equal(a, b)
})

test_that("compute_spi returns value in [0, 1]", {
  set.seed(123)
  res <- compute_spi(rnorm(100), rnorm(100), rnorm(100))
  expect_true(res$SPI_raw >= 0 && res$SPI_raw <= 1)
  expect_true(res$SPI_rarefied >= 0 && res$SPI_rarefied <= 1)
})

test_that("compute_spi length-weighted mode: raw only, with a message", {
  dx <- c(1, 0, 0)
  dy <- c(0, 1, 0)
  dz <- c(0, 0, 1)
  un <- compute_spi(dx, dy, dz)
  expect_message(w <- compute_spi(dx, dy, dz, lengths = c(10, 1, 1)),
                 "unweighted SPI only")
  expect_true(w$SPI_raw > un$SPI_raw)
  expect_true(is.na(w$SPI_rarefied))
  expect_true(is.na(w$SPI_lower))
  expect_true(is.na(w$SPI_vs_random))
})

test_that("compute_spi gives NA below three scars", {
  expect_message(res <- compute_spi(c(1, 0), c(0, 1), c(0, 0)), "Fewer than 3")
  expect_equal(res$n_scar, 2L)
  expect_true(is.na(res$SPI_raw))
  expect_true(is.na(res$SPI_rarefied))
})

test_that("compute_spi drops zero-length and missing vectors", {
  res <- compute_spi(c(1, 0, 1, 0, NA), c(0, 0, 1, 1, 1), c(0, 0, 0, 1, 1))
  expect_equal(res$n_scar, 3L)
})

test_that("compute_spi output columns follow the flags", {
  d <- c(1, 0.5, 0.2, 0.9)
  expect_named(compute_spi(d, rev(d), d^2, rarefied = FALSE, random_interval = FALSE),
               c("n_scar", "SPI_raw"))
  expect_named(compute_spi(d, rev(d), d^2, raw = FALSE),
               c("n_scar", "SPI_rarefied", "SPI_lower", "SPI_upper", "SPI_vs_random"))
  expect_error(compute_spi(d, d, d, raw = FALSE, rarefied = FALSE), "At least one")
})

test_that("rarefied SPI equals raw SPI when m equals the scar count", {
  set.seed(4)
  d <- matrix(rnorm(15), 5, 3)
  res <- compute_spi(d[, 1], d[, 2], d[, 3], m = 5)
  expect_equal(res$SPI_rarefied, res$SPI_raw)
})

test_that("rarefied SPI is the mean SPI over all three-scar combinations", {
  set.seed(5)
  d <- matrix(rnorm(21), 7, 3)
  u <- d / sqrt(rowSums(d^2))
  cmb <- utils::combn(7, 3)
  brute <- mean(apply(cmb, 2, function(i) sqrt(sum(colMeans(u[i, ])^2))))
  expect_equal(compute_spi(d[, 1], d[, 2], d[, 3])$SPI_rarefied, brute)
})

test_that("compute_spi restores the caller's random stream", {
  set.seed(6)
  d <- matrix(rnorm(90), 30, 3)    # C(30, 3) > 2000: a random draw is needed
  set.seed(10)
  before <- .Random.seed
  compute_spi(d[, 1], d[, 2], d[, 3])
  expect_identical(.Random.seed, before)
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
  angle <- compute_spi_angle(c(1, 0), c(0, 1), c(0, 0), lengths = c(5, 1))
  expect_true(is.numeric(angle) && length(angle) == 1)
  expect_true(angle < 45)
})


# ---- compute_ei --------------------------------------------------------------

test_that("compute_ei returns expected columns", {
  res <- compute_ei(c(1, 1, 0.9), c(0, 0, 0.1), c(0, 0, 0))
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1)
  expect_named(res, c("n_scar",
                      paste0("ER_", c("raw", "rarefied", "lower", "upper", "vs_random")),
                      paste0("IR_", c("raw", "rarefied", "lower", "upper", "vs_random")),
                      "lambda1", "lambda2", "lambda3"))
})

test_that("compute_ei: strongly elongated case", {
  # All vectors along X -> ER near 1, IR near 0
  res <- compute_ei(c(1, 1, 1, 1), c(0, 0, 0, 0), c(0, 0, 0, 0))
  expect_equal(res$ER_raw, 1, tolerance = 1e-10)
  expect_equal(res$IR_raw, 0, tolerance = 1e-10)
})

test_that("compute_ei: isotropic case", {
  dx <- c(1, -1,  0,  0,  0,  0)
  dy <- c(0,  0,  1, -1,  0,  0)
  dz <- c(0,  0,  0,  0,  1, -1)
  res <- compute_ei(dx, dy, dz)

  # All eigenvalues equal -> ER = 0, IR = 1
  expect_equal(res$ER_raw, 0, tolerance = 1e-10)
  expect_equal(res$IR_raw, 1, tolerance = 1e-10)
})

test_that("compute_ei: eigenvalues are sorted decreasingly", {
  set.seed(7)
  res <- compute_ei(rnorm(20), rnorm(20), rnorm(20))
  expect_true(res$lambda1 >= res$lambda2)
  expect_true(res$lambda2 >= res$lambda3)
})

test_that("compute_ei: ER and IR in [0, 1]", {
  set.seed(99)
  res <- compute_ei(rnorm(50), rnorm(50), rnorm(50))
  for (v in c(res$ER_raw, res$IR_raw, res$ER_rarefied, res$IR_rarefied)) {
    expect_true(v >= 0 && v <= 1)
  }
})

test_that("rarefied ER and IR are the means over all three-scar combinations", {
  set.seed(8)
  d <- matrix(rnorm(24), 8, 3)
  u <- d / sqrt(rowSums(d^2))
  cmb <- utils::combn(8, 3)
  er_ir <- apply(cmb, 2, function(i) {
    l <- sort(eigen(crossprod(u[i, ]) / 3, symmetric = TRUE)$values, decreasing = TRUE)
    c(1 - l[2] / l[1], l[3] / l[1])
  })
  res <- compute_ei(d[, 1], d[, 2], d[, 3])
  expect_equal(res$ER_rarefied, mean(er_ir[1, ]), tolerance = 1e-12)
  expect_equal(res$IR_rarefied, mean(er_ir[2, ]), tolerance = 1e-12)
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
