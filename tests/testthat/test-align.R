# ==============================================================================
# test-align.R
# Tests for align_scar_batch() and align_morph_batch()
# ==============================================================================

# make_test_scar_data() lives in helper-data.R

# ---- align_scar_batch --------------------------------------------------------

test_that("align_scar_batch returns expected columns and dimensions", {
  raw <- make_test_scar_data()
  res <- align_scar_batch(raw)
  
  # Should add 9 alignment columns
  added_cols <- c("s_x", "s_y", "s_z", "e_x", "e_y", "e_z",
                  "d_x", "d_y", "d_z")
  expect_true(all(added_cols %in% names(res)))
  
  # Row count preserved
  expect_equal(nrow(res), nrow(raw))
  
  # Both IDs present
  expect_equal(sort(unique(res$ID)), c("S001", "S002"))
})

test_that("align_scar_batch produces unit direction vectors", {
  raw <- make_test_scar_data()
  res <- align_scar_batch(raw)
  
  norms <- sqrt(res$d_x^2 + res$d_y^2 + res$d_z^2)
  # All non-zero-length scars should have unit d vectors
  valid <- norms > 1e-10
  expect_true(all(abs(norms[valid] - 1) < 1e-8))
})

test_that("align_scar_batch respects custom id_col", {
  raw <- make_test_scar_data()
  names(raw)[names(raw) == "ID"] <- "specimen"
  res <- align_scar_batch(raw, id_col = "specimen")
  
  expect_true("s_x" %in% names(res))
  expect_equal(nrow(res), nrow(raw))
})


# ---- align_morph_batch -------------------------------------------------------

test_that("align_morph_batch returns expected columns and dimensions", {
  raw <- make_test_scar_data()
  res <- align_morph_batch(raw)
  
  added_cols <- c("s_x", "s_y", "s_z", "e_x", "e_y", "e_z",
                  "d_x", "d_y", "d_z")
  expect_true(all(added_cols %in% names(res)))
  expect_equal(nrow(res), nrow(raw))
})

test_that("align_morph_batch produces unit direction vectors", {
  raw <- make_test_scar_data()
  res <- align_morph_batch(raw)
  
  norms <- sqrt(res$d_x^2 + res$d_y^2 + res$d_z^2)
  valid <- norms > 1e-10
  expect_true(all(abs(norms[valid] - 1) < 1e-8))
})

test_that("align_morph_batch anchors longest scar start at origin XY", {
  raw <- make_test_scar_data()
  res <- align_morph_batch(raw)
  
  for (sid in unique(res$ID)) {
    sub <- res[res$ID == sid, ]
    lens <- sqrt((sub$e_x - sub$s_x)^2 +
                   (sub$e_y - sub$s_y)^2 +
                   (sub$e_z - sub$s_z)^2)
    longest_idx <- which.max(lens)
    # The longest scar's start-point should be at x=0, y=0
    expect_equal(sub$s_x[longest_idx], 0, tolerance = 1e-10)
    expect_equal(sub$s_y[longest_idx], 0, tolerance = 1e-10)
  }
})