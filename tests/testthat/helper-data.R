# Shared test data -----------------------------------------------------------

# Minimal raw scar data frame: two specimens with known start/end points
make_test_scar_data <- function(n_scars = 5, seed = 42) {
  set.seed(seed)
  # Two specimens with known start/end points
  specimens <- list(
    list(
      id = "S001",
      pos = c(10, 20, 30), norm = c(0, 0, 1),
      starts = matrix(rnorm(n_scars * 3), n_scars, 3),
      ends   = matrix(rnorm(n_scars * 3) + 2, n_scars, 3)
    ),
    list(
      id = "S002",
      pos = c(-5, 10, 15), norm = c(0.5, 0.5, 0.707),
      starts = matrix(rnorm(n_scars * 3, sd = 2), n_scars, 3),
      ends   = matrix(rnorm(n_scars * 3, sd = 2) + 3, n_scars, 3)
    )
  )

  do.call(rbind, lapply(specimens, function(sp) {
    data.frame(
      ID      = sp$id,
      Pos_X   = sp$pos[1],
      Pos_Y   = sp$pos[2],
      Pos_Z   = sp$pos[3],
      Norm_X  = sp$norm[1],
      Norm_Y  = sp$norm[2],
      Norm_Z  = sp$norm[3],
      Scar_ID = paste0("s", seq_len(n_scars)),
      Start_X = sp$starts[, 1],
      Start_Y = sp$starts[, 2],
      Start_Z = sp$starts[, 3],
      End_X   = sp$ends[, 1],
      End_Y   = sp$ends[, 2],
      End_Z   = sp$ends[, 3],
      Typology = if (sp$id == "S001") "Levallois" else "Discoid",
      stringsAsFactors = FALSE
    )
  }))
}
