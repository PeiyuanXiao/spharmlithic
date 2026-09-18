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

test_that(".cilm_to_json rejects complex coefficients", {
  coeff <- array(complex(real = rnorm(18), imaginary = 0.1),
                 dim = c(2, 3, 3))
  expect_error(
    spharmlithic:::.cilm_to_json(coeff, lmax_out = 2, digits = 4),
    "Complex"
  )
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


# ---- export_spharm_html output -----------------------------------------------

# Stand-in for a spharm_from_*() result: a named list of specimens, each with
# a (2, lmax + 1, lmax + 1) coefficient array
fake_spharm <- function(ids, lmax = 6, seed = 1) {
  set.seed(seed)
  n <- lmax + 1
  stats::setNames(lapply(ids, function(id) {
    cf <- array(round(rnorm(2 * n * n), 3), dim = c(2, n, n))
    cf[1, 1, 1] <- 1
    list(coefficients = cf, power_spectrum = rep(0, n))
  }), ids)
}

# The specimen records the viewer embeds as `const DATA=[...];`
read_viewer_data <- function(path) {
  html <- readLines(path, encoding = "UTF-8", warn = FALSE)
  line <- trimws(html[startsWith(trimws(html), "const DATA=")])
  jsonlite::fromJSON(sub(";$", "", sub("^const DATA=", "", line)),
                     simplifyVector = FALSE)
}

test_that("export_spharm_html embeds every specimen from both tracks", {
  morph <- fake_spharm(c("A", "B"))
  scar  <- fake_spharm(c("B", "C"), seed = 2)
  out   <- tempfile(fileext = ".html")

  res <- expect_invisible(
    export_spharm_html(morph, scar, out_path = out, verbose = FALSE)
  )
  expect_identical(res, out)
  expect_false(any(grepl("{{", readLines(out, warn = FALSE), fixed = TRUE)))

  recs <- read_viewer_data(out)
  expect_identical(vapply(recs, `[[`, character(1), "id"), c("A", "B", "C"))
  expect_null(recs[[1]]$scar)     # A is morphology only
  expect_null(recs[[3]]$morph)    # C is scar only

  # Coefficients arrive as cilm[layer][l][m]
  cf <- morph$B$coefficients
  expect_equal(recs[[2]]$morph[[1]][[3]][[2]], cf[1, 3, 2])
  expect_equal(recs[[2]]$morph[[2]][[4]][[3]], cf[2, 4, 3])
  expect_equal(recs[[2]]$scar[[1]][[1]][[1]], 1)
})

test_that("export_spharm_html truncates to lmax and fills in the title", {
  out <- tempfile(fileext = ".html")
  export_spharm_html(scar = fake_spharm("A", lmax = 10), out_path = out,
                     lmax = 4, title = "Clarkson cores", verbose = FALSE)

  cilm <- read_viewer_data(out)[[1]]$scar
  expect_length(cilm[[1]], 5)
  expect_length(cilm[[1]][[5]], 5)

  html <- readLines(out, encoding = "UTF-8", warn = FALSE)
  expect_true(any(grepl("const LMAX=4;", html, fixed = TRUE)))
  expect_true(any(grepl("Clarkson cores", html, fixed = TRUE)))
})

test_that("export_spharm_html writes metadata and escapes quotes", {
  scar <- fake_spharm(c("A", "B", 'Core "7"', "D"))
  meta <- data.frame(
    ID       = c("A", "B", 'Core "7"'),
    Site     = c('Say "hi"', "back\\slash", NA),
    Typology = c("Levallois", "Discoid", "Discoid")
  )
  out <- tempfile(fileext = ".html")
  export_spharm_html(scar = scar, meta = meta, out_path = out,
                     verbose = FALSE)

  recs <- read_viewer_data(out)
  by_id <- stats::setNames(recs, vapply(recs, `[[`, character(1), "id"))
  expect_identical(by_id$A$meta$Site, 'Say "hi"')
  expect_identical(by_id$B$meta$Site, "back\\slash")
  expect_null(by_id[['Core "7"']]$meta$Site)          # NA is left out
  expect_identical(by_id[['Core "7"']]$meta$Typology, "Discoid")
  expect_length(by_id$D$meta, 0)                      # not in meta
})

test_that("export_spharm_html skips specimens without coefficients", {
  morph <- fake_spharm(c("A", "B"))
  morph$B$coefficients <- NULL
  out <- tempfile(fileext = ".html")

  msgs <- capture_messages(export_spharm_html(morph, out_path = out))
  expect_true(any(grepl("B -- skipped (no data)", msgs, fixed = TRUE)))
  expect_identical(vapply(read_viewer_data(out), `[[`, character(1), "id"),
                   "A")

  morph$A$coefficients <- NULL
  expect_error(
    export_spharm_html(morph, out_path = tempfile(), verbose = FALSE),
    "No valid specimens"
  )
})

test_that("export_spharm_html creates the output folder and reports progress", {
  out  <- file.path(tempfile(), "nested", "viewer.html")
  msgs <- capture_messages(
    export_spharm_html(scar = fake_spharm("A"), out_path = out)
  )
  expect_true(file.exists(out))
  expect_true(any(grepl("Exported: ", msgs, fixed = TRUE)))
})


# ---- Viewer JavaScript: Legendre functions ----------------------------------

test_that("viewer computePlm() uses pyshtools 4pi normalization", {
  skip_if_not_installed("V8")

  html  <- readLines(system.file("templates", "spharm_viewer.html",
                                 package = "spharmlithic"),
                     encoding = "UTF-8")
  start <- which(startsWith(html, "function computePlm"))
  end   <- start - 1 + which(startsWith(html[start:length(html)], "}"))[1]
  ctx   <- V8::v8()
  ctx$eval(paste(html[start:end], collapse = "\n"))

  th  <- c(0.3, 1.1, 2.4)
  plm <- function(l, m) {
    ctx$get(sprintf("Array.from(computePlm(2, %s)[%d][%d])",
                    jsonlite::toJSON(th, digits = NA), l, m))
  }

  # Closed forms of the 4pi-normalized functions (no Condon-Shortley phase)
  expect_equal(plm(1, 0), sqrt(3) * cos(th))
  expect_equal(plm(1, 1), sqrt(3) * sin(th))
  expect_equal(plm(2, 0), sqrt(5) / 2 * (3 * cos(th)^2 - 1))
  expect_equal(plm(2, 1), sqrt(15) * sin(th) * cos(th))
  expect_equal(plm(2, 2), sqrt(15) / 2 * sin(th)^2)
})