# ==============================================================================
# test-spatiodirection.R
# Tests for spharm_spatiodirection()
# ==============================================================================

test_that("spharm_spatiodirection returns 22 features per core", {
  d <- make_cores()
  sp <- spharm_spatiodirection(d, keep = "Typology", quiet = TRUE)
  feats <- c(paste0("dir_", 1:6), paste0("start_", 2:6), paste0("same_", 1:6),
             paste0("diff_", 2:6))
  expect_equal(attr(sp, "features"), feats)
  expect_named(sp, c("ID", "Typology", "n_scar", feats))
  expect_equal(nrow(sp), 3L)
  expect_false(anyNA(sp[, feats]))
})

test_that("`features` selects blocks of features", {
  sp <- spharm_spatiodirection(make_cores(), features = c("dir", "same"), quiet = TRUE)
  expect_equal(attr(sp, "features"), c(paste0("dir_", 1:6), paste0("same_", 1:6)))
  expect_error(spharm_spatiodirection(make_cores(), features = "xyz"), "one or more of")
})

test_that("dir_1 is the mean pairwise cosine, sqrt(3/2) (n SPI^2 - 1) / (n - 1)", {
  d <- make_cores(n = 15)
  sp <- spharm_spatiodirection(d, quiet = TRUE)
  spi <- spi_ei_batch(d, rarefied = FALSE, random_interval = FALSE, quiet = TRUE)$SPI_raw
  n <- 15
  expect_equal(sp$dir_1, sqrt(3 / 2) * (n * spi^2 - 1) / (n - 1), tolerance = 1e-12)
})

test_that("degree-l features are pair means of Legendre polynomials", {
  d <- make_cores(n = 9)
  sp <- spharm_spatiodirection(d, quiet = TRUE)
  v <- as.matrix(d[, c("End_X", "End_Y", "End_Z")] - d[, c("Start_X", "Start_Y", "Start_Z")])
  v <- v / sqrt(rowSums(v^2))
  x <- tcrossprod(v)[upper.tri(diag(9))]
  P2 <- (3 * x^2 - 1) / 2
  expect_equal(sp$dir_2, sqrt(5 / 2) * mean(P2), tolerance = 1e-12)
  s <- as.matrix(d[, c("Start_X", "Start_Y", "Start_Z")])
  u <- sweep(s, 2, colMeans(s))
  u <- u / sqrt(rowSums(u^2))
  expect_equal(sp$same_1, sqrt(3 / 2) * mean(rowSums(v * u)), tolerance = 1e-12)
})

test_that("features do not change when a core is rotated and moved", {
  d <- make_cores()
  d2 <- rigid_move(d, random_rotation(), c(10, 20, -5))
  a <- spharm_spatiodirection(d, quiet = TRUE)
  b <- spharm_spatiodirection(d2, quiet = TRUE)
  f <- attr(a, "features")
  expect_equal(as.matrix(b[, f]), as.matrix(a[, f]), tolerance = 1e-10)
})

test_that("cores with fewer than three scars are NA", {
  d <- make_cores(n = c(10, 2))
  expect_message(sp <- spharm_spatiodirection(d), "fewer than 3")
  expect_true(all(is.na(unlist(sp[2, attr(sp, "features")]))))
})

test_that("for random directions the features centre on zero", {
  set.seed(9)
  d <- do.call(rbind, lapply(1:150, function(i) {
    s <- matrix(rnorm(24), 8, 3)
    e <- s + matrix(rnorm(24), 8, 3)
    data.frame(ID = i, Start_X = s[, 1], Start_Y = s[, 2], Start_Z = s[, 3],
               End_X = e[, 1], End_Y = e[, 2], End_Z = e[, 3])
  }))
  sp <- spharm_spatiodirection(d, features = c("dir", "same"), quiet = TRUE)
  m <- colMeans(sp[, attr(sp, "features")])
  expect_true(all(abs(m) < 0.06))
})
