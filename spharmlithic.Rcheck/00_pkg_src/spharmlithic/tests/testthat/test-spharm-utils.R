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
