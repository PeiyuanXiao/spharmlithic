test_that("spharm_to_dataframe handles minimal input", {
  fake_result <- list(
    A = list(
      coefficients   = array(rnorm(2 * 3 * 3), dim = c(2, 3, 3)),
      power_spectrum = c(1, 0.5, 0.25)
    ),
    B = list(
      coefficients   = array(rnorm(2 * 3 * 3), dim = c(2, 3, 3)),
      power_spectrum = c(1, 0.4, 0.3)
    )
  )

  df_full <- spharm_to_dataframe(fake_result)
  expect_s3_class(df_full, "tbl_df")
  expect_equal(df_full$ID, c("A", "B"))
  expect_true(all(c("power_l0", "power_l1", "power_l2") %in% names(df_full)))
  expect_true(any(grepl("^coeff_", names(df_full))))
  expect_equal(df_full$power_l0, c(1, 1))

  df_pwr <- spharm_to_dataframe(fake_result, include_coeffs = FALSE)
  expect_false(any(grepl("^coeff_", names(df_pwr))))
})

test_that("spharm_to_dataframe rejects bad input", {
  expect_error(spharm_to_dataframe(list()))
  expect_error(spharm_to_dataframe(list(list(power_spectrum = 1))))  # unnamed
})

test_that("spharm_to_dataframe rejects complex coefficients", {
  old_result <- list(A = list(
    coefficients   = array(complex(real = 1, imaginary = 0.1), dim = c(2, 3, 3)),
    power_spectrum = c(1, 0.5, 0.25)
  ))
  expect_error(spharm_to_dataframe(old_result), "Complex")
})

test_that("spharm_reconstruct rejects complex coefficients", {
  cf <- array(complex(real = 1, imaginary = 0.1), dim = c(2, 2, 2))
  expect_error(spharm_reconstruct(cf), "Complex")
})


# ---- Known-answer tests for the spherical harmonic transforms ---------------
# Coefficients use pyshtools 4pi normalization: C00 = 1 is the constant 1,
# C10 = 1 is sqrt(3) * cos(colat) = sqrt(3) * z, C11 = 1 is sqrt(3) * x.

test_that("spharm_reconstruct places density on the correct lat/lon grid", {
  skip_if_no_python_core()

  cf <- array(0, dim = c(2, 2, 2))
  cf[1, 1, 1] <- 2   # C00
  cf[1, 2, 1] <- 1   # C10
  rec <- spharm_reconstruct(cf, grid_size = 16)
  expect_equal(rec$lat[1], pi / 2)
  expect_equal(rec$lon[1], 0)
  expect_equal(rec$density[, 1], 2 + sqrt(3) * sin(rec$lat), tolerance = 1e-8)

  cf <- array(0, dim = c(2, 2, 2))
  cf[1, 1, 1] <- 2   # C00
  cf[1, 2, 2] <- 1   # C11
  rec <- spharm_reconstruct(cf, grid_size = 16)
  expect_equal(as.vector(rec$density), 2 + sqrt(3) * rec$xyz[, "x"],
               tolerance = 1e-8)
})

test_that("mesh expansion returns real coefficients that reconstruct the surface", {
  skip_if_no_python_core()

  # Radius on the (n, n) grid built by the mesh pipeline. The y and x*y
  # terms are pure sine terms, which the real part of complex
  # coefficients would lose.
  n     <- 32
  colat <- (seq_len(n) - 1) * pi / n
  lon   <- (seq_len(n) - 1) * 2 * pi / n
  x <- outer(sin(colat), cos(lon))
  y <- outer(sin(colat), sin(lon))
  r <- 1 + 0.15 * x + 0.2 * y + 0.2 * x * y

  cf <- spharmlithic:::sh_py$spherical_harmonics$compute_spherical_harmonics(r)
  expect_false(is.complex(cf))
  expect_equal(dim(cf), c(2L, n / 2, n / 2))
  expect_equal(cf[2, 2, 2], 0.2 / sqrt(3), tolerance = 1e-8)  # S11: the y term

  rec <- spharm_reconstruct(cf, grid_size = 32)
  xyz <- rec$xyz
  expect_equal(
    as.vector(rec$density),
    1 + 0.15 * xyz[, "x"] + 0.2 * xyz[, "y"] + 0.2 * xyz[, "x"] * xyz[, "y"],
    tolerance = 1e-8
  )
})
