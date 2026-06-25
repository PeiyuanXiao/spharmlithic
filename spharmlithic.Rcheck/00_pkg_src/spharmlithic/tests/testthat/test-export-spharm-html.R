# ==============================================================================
# test-export-spharm-html.R
# Tests for export_spharm_html() — JSON logic only, no Python required
# ==============================================================================

# ---- Internal helper: .cilm_to_json -----------------------------------------

test_that(".cilm_to_json converts a 3D array to valid JSON", {
  coeff <- array(1:18, dim = c(2, 3, 3))
  json <- spharmlithic:::.cilm_to_json(coeff, lmax_out = 2, digits = 4)
  
  expect_type(json, "character")
  expect_true(nchar(json) > 0)
  
  # Should parse as valid JSON
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  expect_length(parsed, 2)         # two layers (cosine, sine)
  expect_length(parsed[[1]], 3)    # lmax+1 = 3 rows
  expect_length(parsed[[1]][[1]], 3) # lmax+1 = 3 columns
})

test_that(".cilm_to_json truncates to requested lmax", {
  coeff <- array(rnorm(2 * 6 * 6), dim = c(2, 6, 6))  # lmax=5
  json <- spharmlithic:::.cilm_to_json(coeff, lmax_out = 3, digits = 4)
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  
  # Should have lmax_out+1 = 4 rows per layer
  
  expect_length(parsed[[1]], 4)
  expect_length(parsed[[2]], 4)
})

test_that(".cilm_to_json handles complex coefficients", {
  coeff <- array(complex(real = rnorm(18), imaginary = 1e-15),
                 dim = c(2, 3, 3))
  json <- spharmlithic:::.cilm_to_json(coeff, lmax_out = 2, digits = 4)
  expect_type(json, "character")
  expect_true(nchar(json) > 0)
})

test_that(".cilm_to_json returns NULL for invalid input", {
  expect_null(spharmlithic:::.cilm_to_json(NULL, 2, 4))
  expect_null(spharmlithic:::.cilm_to_json(matrix(1:4, 2, 2), 2, 4))
  # Wrong first dim
  expect_null(spharmlithic:::.cilm_to_json(array(1:27, dim = c(3, 3, 3)), 2, 4))
})


# ---- export_spharm_html input validation ------------------------------------

test_that("export_spharm_html errors when both morph and scar are NULL", {
  expect_error(
    export_spharm_html(morph = NULL, scar = NULL, out_path = tempfile()),
    "At least one"
  )
})

test_that("export_spharm_html errors on lmax < 1", {
  fake <- list(A = list(coefficients = array(0, dim = c(2, 2, 2)),
                        power_spectrum = c(0, 0)))
  expect_error(
    export_spharm_html(scar = fake, out_path = tempfile(), lmax = 0),
    "lmax"
  )
})

test_that("export_spharm_html errors when results are unnamed", {
  fake <- list(list(coefficients = array(0, dim = c(2, 2, 2)),
                    power_spectrum = c(0, 0)))
  expect_error(
    export_spharm_html(scar = fake, out_path = tempfile()),
    "No specimen IDs"
  )
})

test_that("export_spharm_html errors when meta lacks ID column", {
  fake <- list(A = list(coefficients = array(0, dim = c(2, 3, 3)),
                        power_spectrum = c(0, 0, 0)))
  expect_error(
    export_spharm_html(scar = fake, meta = data.frame(name = "A"),
                       out_path = tempfile()),
    "ID"
  )
})