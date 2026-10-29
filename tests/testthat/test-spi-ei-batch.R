# ==============================================================================
# test-spi-ei-batch.R
# Tests for spi_ei_batch(), random_interval() and check_scars()
# ==============================================================================

# ---- random_interval ---------------------------------------------------------

test_that("random_interval returns one row per parameter and n", {
  ri <- random_interval(c("SPI", "ir"), n = c(5, 20))
  expect_s3_class(ri, "tbl_df")
  expect_named(ri, c("parameter", "n", "lower", "upper"))
  expect_equal(ri$parameter, c("SPI", "SPI", "IR", "IR"))
  expect_true(all(ri$lower < ri$upper))
})

test_that("the SPI interval narrows as the scar count grows", {
  ri <- random_interval("SPI", n = c(3, 6, 20, 100))
  expect_true(all(diff(ri$upper - ri$lower) < 0))
})

test_that("random_interval: NA below three scars and for ER/IR above the table", {
  expect_true(is.na(random_interval("SPI", 2)$lower))
  expect_message(ri <- random_interval("ER", 300), "tabulated up to")
  expect_true(is.na(ri$lower))
})

test_that("the exact SPI distribution continues the table smoothly", {
  ri <- random_interval("SPI", c(200, 201))
  expect_equal(ri$upper[2], ri$upper[1], tolerance = 0.01)
  expect_true(ri$upper[2] < ri$upper[1])
})

test_that("random_interval rejects unknown parameters and non-counts", {
  expect_error(random_interval("XYZ", 5), "Unknown parameter")
  expect_error(random_interval("SPI", 5.5), "whole numbers")
})

# ---- spi_ei_batch ---------------------------------------------------------------

test_that("spi_ei_batch returns one row per core with all columns", {
  d <- make_cores()
  res <- spi_ei_batch(d, keep = "Typology", quiet = TRUE)
  expect_equal(res$ID, c("core1", "core2", "core3"))
  expect_equal(res$Typology, c("A", "B", "A"))
  expect_equal(res$n_scar, c(12L, 8L, 20L))
  for (p in c("SPI", "ER", "IR")) {
    expect_true(all(paste0(p, c("_raw", "_rarefied", "_lower", "_upper", "_vs_random")) %in%
                      names(res)))
  }
})

test_that("spi_ei_batch agrees with compute_spi and compute_ei core by core", {
  d <- make_cores()
  res <- spi_ei_batch(d, quiet = TRUE)
  for (i in seq_len(nrow(res))) {
    g <- d[d$ID == res$ID[i], ]
    dx <- g$End_X - g$Start_X; dy <- g$End_Y - g$Start_Y; dz <- g$End_Z - g$Start_Z
    s <- compute_spi(dx, dy, dz)
    e <- compute_ei(dx, dy, dz)
    expect_equal(res$SPI_raw[i], s$SPI_raw)
    expect_equal(res$SPI_rarefied[i], s$SPI_rarefied)
    expect_equal(res$SPI_vs_random[i], s$SPI_vs_random)
    expect_equal(res$ER_rarefied[i], e$ER_rarefied)
    expect_equal(res$IR_upper[i], e$IR_upper)
  }
})

test_that("spi_ei_batch gives the same values for raw and aligned coordinates", {
  d <- make_cores()
  a <- align_scar_batch(d)
  raw <- spi_ei_batch(d, quiet = TRUE)
  ali <- spi_ei_batch(a, start_cols = c("s_x", "s_y", "s_z"),
                      end_cols = c("e_x", "e_y", "e_z"), quiet = TRUE)
  expect_equal(ali$SPI_raw, raw$SPI_raw, tolerance = 1e-12)
  expect_equal(ali$ER_raw, raw$ER_raw, tolerance = 1e-12)
  expect_equal(ali$IR_rarefied, raw$IR_rarefied, tolerance = 1e-12)
})

test_that("spi_ei_batch reports dropped scars and cores below three scars", {
  d <- make_cores(n = c(6, 3))
  d$End_X[d$ID == "core2"][1] <- d$Start_X[d$ID == "core2"][1]
  d$End_Y[d$ID == "core2"][1] <- d$Start_Y[d$ID == "core2"][1]
  d$End_Z[d$ID == "core2"][1] <- d$Start_Z[d$ID == "core2"][1]
  d$Start_X[2] <- NA
  expect_message(res <- spi_ei_batch(d), "zero length")
  expect_equal(res$n_scar, c(5L, 2L))
  expect_true(all(is.na(unlist(res[2, c("SPI_raw", "ER_rarefied", "IR_lower")]))))
  expect_message(spi_ei_batch(d), "fewer than 3 usable scars are NA: core2")
})

test_that("spi_ei_batch: min_length leaves out short scars", {
  d <- make_cores(n = c(10))
  L <- sqrt((d$End_X - d$Start_X)^2 + (d$End_Y - d$Start_Y)^2 + (d$End_Z - d$Start_Z)^2)
  res <- spi_ei_batch(d, min_length = stats::median(L), quiet = TRUE)
  expect_equal(res$n_scar, sum(L >= stats::median(L)))
})

test_that("spi_ei_batch flags Monte Carlo cores and restores the random stream", {
  d <- make_cores(n = c(30))
  set.seed(11)
  before <- .Random.seed
  expect_message(spi_ei_batch(d), "random combinations")
  expect_identical(.Random.seed, before)
})

test_that("spi_ei_batch checks its input", {
  d <- make_cores()
  expect_error(spi_ei_batch(d[, -3]), "lacks column")
  expect_error(spi_ei_batch(d, keep = "Start_X"), "varies within core")
  expect_error(spi_ei_batch(d, min_length = -1), "non-negative")
})

# ---- check_scars ------------------------------------------------------------------

test_that("check_scars summarises the table in one message", {
  d <- make_cores(n = c(5, 2))
  d$Start_X[1] <- NA
  expect_message(res <- check_scars(d, keep = "Typology"), "2 core\\(s\\)")
  expect_equal(res$n_rows, c(5L, 2L))
  expect_equal(res$n_missing, c(1L, 0L))
  expect_equal(res$n_scar, c(4L, 2L))
  expect_equal(res$enough, c(TRUE, FALSE))
})
