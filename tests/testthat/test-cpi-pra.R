# ==============================================================================
# test-cpi-pra.R
# Tests for fit_intersection_plane(), cpi_pra_batch(), cpi_pra_per_scar() and
# read_outline_iges()
# ==============================================================================

# ---- fit_intersection_plane ------------------------------------------------------

test_that("the plane of a flat circle is z = 0 with the centre at the origin", {
  pl <- fit_intersection_plane(disc_outline())
  expect_s3_class(pl, "intersection_plane")
  expect_equal(abs(pl$normal), c(0, 0, 1), tolerance = 1e-12)
  expect_equal(pl$centre, c(0, 0, 0), tolerance = 1e-10)
  expect_equal(pl$flatness, 0, tolerance = 1e-12)
  expect_equal(nrow(pl$loop), 201L)          # closed: first point repeated
})

test_that("the centre is the area centroid, not the mean of the points", {
  # a square sampled densely on one side only
  side <- seq(0, 10, length.out = 200)
  sq <- rbind(cbind(side, 0, 0), cbind(10, 10, 0), cbind(0, 10, 0))
  pl <- fit_intersection_plane(sq)
  expect_equal(pl$centre[1:2], c(5, 5), tolerance = 1e-10)
  expect_true(mean(sq[, 2]) < 1)
})

test_that("fit_intersection_plane warns about points out of order", {
  set.seed(1)
  ol <- disc_outline()[sample(200), ]
  expect_warning(fit_intersection_plane(ol), "crosses itself")
})

test_that("fit_intersection_plane warns about an outline with a missing stretch", {
  arc <- disc_outline()[1:150, ]           # three quarters of the circle
  expect_warning(fit_intersection_plane(arc), "not closed")
  expect_silent(fit_intersection_plane(disc_outline()[seq(1, 200, 20), ]))
})

test_that("fit_intersection_plane rejects bad outlines", {
  expect_error(fit_intersection_plane(matrix(1:6, 2, 3)), "four distinct")
  expect_error(fit_intersection_plane(matrix(1:8, 4, 2)), "three")
  expect_error(fit_intersection_plane(data.frame(a = 1:5, b = "x", c = 1:5)), "three numeric")
})

# ---- cpi_pra_batch -----------------------------------------------------------------

test_that("centripetal scars read CPI 180 and their dip as PRA", {
  d <- make_disc(cpi = 180, pra = 20)
  res <- cpi_pra_batch(d, list(disc = disc_outline()), quiet = TRUE)
  expect_equal(res$CPI_mean, 180, tolerance = 1e-8)
  expect_equal(res$PRA_mean, 20, tolerance = 1e-8)
  expect_equal(res$CPI_GMD, 0, tolerance = 1e-8)
  expect_equal(res$PRA_GMD, 0, tolerance = 1e-8)
  expect_true(res$frame_ok)
})

test_that("CPI reads 90 along the edge and 0 away from the centre; PRA 0 when flat", {
  edge <- cpi_pra_batch(make_disc(cpi = 90, pra = 0), list(disc = disc_outline()), quiet = TRUE)
  out <- cpi_pra_batch(make_disc(cpi = 0, pra = 0), list(disc = disc_outline()), quiet = TRUE)
  expect_equal(edge$CPI_mean, 90, tolerance = 1e-8)
  expect_equal(out$CPI_mean, 0, tolerance = 1e-6)
  expect_equal(edge$PRA_mean, 0, tolerance = 1e-8)
})

test_that("GMD equals the mean absolute difference over all pairs", {
  set.seed(2)
  d <- make_disc(k = 9, cpi = 150, pra = 30)
  d$End_Z <- d$End_Z + rnorm(9)
  per <- cpi_pra_per_scar(d, disc_outline())
  res <- cpi_pra_batch(d, list(disc = disc_outline()), quiet = TRUE)
  expect_equal(res$PRA_GMD, mean(stats::dist(per$PRA)))
  expect_equal(res$CPI_GMD, mean(stats::dist(per$CPI)))
  expect_equal(res$PRA_mean, mean(per$PRA))
})

test_that("CPI and PRA do not change when the core and its outline are moved", {
  d <- make_disc(k = 10, cpi = 140, pra = 35)
  R <- random_rotation()
  shift <- c(5, -3, 12)
  d2 <- rigid_move(d, R, shift)
  ol2 <- sweep(disc_outline() %*% t(R), 2, shift, "+")
  a <- cpi_pra_batch(d, list(disc = disc_outline()), quiet = TRUE)
  b <- cpi_pra_batch(d2, list(disc = ol2), quiet = TRUE)
  expect_equal(b$CPI_mean, a$CPI_mean, tolerance = 1e-8)
  expect_equal(b$PRA_GMD, a$PRA_GMD, tolerance = 1e-8)
})

test_that("cpi_pra_batch flags scars and outline in different frames", {
  d <- rigid_move(make_disc(), diag(3), c(100, 0, 0))
  expect_warning(res <- cpi_pra_batch(d, list(disc = disc_outline()), quiet = TRUE),
                 "frame check")
  expect_false(res$frame_ok)
})

test_that("cpi_pra_batch reports cores without an outline and unmatched outlines", {
  d <- rbind(make_disc(id = "a"), make_disc(id = "b"))
  ol <- list(a = disc_outline(), zz = disc_outline())
  expect_message(res <- cpi_pra_batch(d, ol), "1 core\\(s\\) have no outline")
  expect_message(cpi_pra_batch(d, ol), "match no core ID: zz")
  expect_equal(res$ID, "a")
  expect_error(cpi_pra_batch(d, list(q = disc_outline())), "No outline matches")
  expect_error(cpi_pra_batch(d, disc_outline()), "named by core ID")
})

test_that("cpi_pra_batch accepts a fitted plane in place of an outline", {
  d <- make_disc()
  a <- cpi_pra_batch(d, list(disc = disc_outline()), quiet = TRUE)
  b <- cpi_pra_batch(d, list(disc = fit_intersection_plane(disc_outline())), quiet = TRUE)
  expect_equal(a, b)
})

test_that("cpi_pra_per_scar returns one row per scar and refuses several cores", {
  d <- make_disc(k = 8)
  d$Scar_ID <- paste0("s", 1:8)
  per <- cpi_pra_per_scar(d, disc_outline())
  expect_named(per, c("scar", "CPI", "PRA", "dist_to_outline"))
  expect_equal(per$scar, d$Scar_ID)
  expect_true(all(per$dist_to_outline < 1e-6))
  expect_error(cpi_pra_per_scar(rbind(d, transform(d, ID = "x")), disc_outline()),
               "more than one core")
})

# ---- read_outline_iges ---------------------------------------------------------------

test_that("read_outline_iges chains B-spline segments into one closed loop", {
  sq <- rbind(c(0, 0, 0), c(10, 0, 0), c(10, 10, 0), c(0, 10, 0), c(0, 5, 0))
  f <- tempfile(fileext = ".igs")
  write_iges_polyline(f, sq)
  ol <- read_outline_iges(f, n_per_seg = 20)
  expect_named(ol, sub("\\.igs$", "", basename(f)))
  loop <- ol[[1]]
  expect_equal(colnames(loop), c("X", "Y", "Z"))
  expect_equal(unname(loop[1, ]), unname(loop[nrow(loop), ]))
  expect_equal(attr(loop, "n_segments"), 2L)
  pl <- fit_intersection_plane(loop)
  expect_equal(pl$centre[1:2], c(5, 5), tolerance = 1e-8)
})

test_that("read_outline_iges uses `ids` and leaves out unreadable files", {
  sq <- rbind(c(0, 0, 0), c(10, 0, 0), c(10, 10, 0), c(0, 10, 0))
  f1 <- tempfile(fileext = ".igs")
  f2 <- tempfile(fileext = ".igs")
  write_iges_polyline(f1, sq)
  write_iges_polyline(f2, sq, extra_type = 100L)
  expect_warning(ol <- read_outline_iges(c(f1, f2), ids = c("A", "B")),
                 "circular arc")
  expect_named(ol, "A")
  expect_error(read_outline_iges(f1, ids = c("A", "B")), "one ID per file")
})
