# ==============================================================================
# test-export-alignment-html.R
# Tests for export_alignment_html_svd() / export_alignment_html_lin2024() and
# the Plotly helpers they use — no Python required
# ==============================================================================

# ---- Helpers -----------------------------------------------------------------

# Parse the `allPanels["<id>"] = {...};` line the page embeds for one specimen
# into a list of panels (p0, p1, ...), each holding Plotly `data` and `layout`
read_panels <- function(html, id) {
  prefix <- sprintf('allPanels["%s"] = ', id)
  line   <- trimws(html[startsWith(trimws(html), prefix)])
  json   <- sub(";$", "", substring(line, nchar(prefix) + 1))
  jsonlite::fromJSON(json, simplifyVector = FALSE)
}

count_types <- function(panel, type) {
  sum(vapply(panel$data, function(tr) tr$type, character(1)) == type)
}

# Line traces drawn in one colour (scars, arrows)
lines_in <- function(panel, color) {
  Filter(function(tr) identical(tr$line$color, color), panel$data)
}

direction_of <- function(tr) {
  d <- c(diff(unlist(tr$x)), diff(unlist(tr$y)), diff(unlist(tr$z)))
  d / sqrt(sum(d^2))
}

length_of <- function(tr) {
  sqrt(diff(unlist(tr$x))^2 + diff(unlist(tr$y))^2 + diff(unlist(tr$z))^2)
}


# ---- export_alignment_html_svd -----------------------------------------------

test_that("export_alignment_html_svd writes one page covering every specimen", {
  raw <- make_test_scar_data()
  out <- tempfile(fileext = ".html")

  res <- expect_invisible(export_alignment_html_svd(raw, out))
  expect_identical(res, out)
  expect_true(file.exists(out))

  html <- readLines(out, warn = FALSE)
  expect_true(any(grepl("plotly-2.27.0.min.js", html, fixed = TRUE)))
  expect_true(any(grepl("Core Alignment Pipeline (SVD normal)", html,
                        fixed = TRUE)))
  for (id in c("S001", "S002")) {
    option <- sprintf('<option value="%s">%s</option>', id, id)
    expect_true(any(grepl(option, html, fixed = TRUE)))
    expect_named(read_panels(html, id), c("p0", "p1", "p2", "p3"))
  }
})

test_that("SVD panels draw every scar, the normal, the plane and the main axis", {
  raw <- make_test_scar_data(n_scars = 5)
  out <- tempfile(fileext = ".html")
  export_alignment_html_svd(raw, out)
  p <- read_panels(readLines(out, warn = FALSE), "S001")

  # Each scar is a line plus a cone; the normal arrow is one more of each
  for (k in c("p0", "p1", "p2")) {
    expect_identical(count_types(p[[k]], "scatter3d"), 6L)
    expect_identical(count_types(p[[k]], "cone"), 6L)
    expect_identical(count_types(p[[k]], "mesh3d"), 1L)
    expect_length(lines_in(p[[k]], "steelblue"), 5)
  }
  # The last panel adds the PCA main axis (orange arrow)
  expect_identical(count_types(p$p3, "scatter3d"), 7L)
  expect_identical(count_types(p$p3, "cone"), 7L)
  expect_length(lines_in(p$p3, "orange"), 1)
})

test_that("SVD alignment brings the normal onto Z and the main axis onto X", {
  raw <- make_test_scar_data()
  out <- tempfile(fileext = ".html")
  export_alignment_html_svd(raw, out)
  p <- read_panels(readLines(out, warn = FALSE), "S001")

  # Step 0 shows the SVD normal of the unit scar directions (z >= 0)
  df <- raw[raw$ID == "S001", ]
  d  <- as.matrix(df[, c("End_X", "End_Y", "End_Z")]) -
    as.matrix(df[, c("Start_X", "Start_Y", "Start_Z")])
  normal <- svd(d / sqrt(rowSums(d^2)))$v[, 3]
  if (normal[3] < 0) normal <- -normal
  expect_equal(direction_of(lines_in(p$p0, "red")[[1]]), normal,
               tolerance = 1e-3)

  # After step 1 the normal points along Z; after step 3 the main axis is X
  for (k in c("p1", "p2", "p3")) {
    expect_equal(direction_of(lines_in(p[[k]], "red")[[1]]), c(0, 0, 1),
                 tolerance = 1e-3)
  }
  expect_equal(direction_of(lines_in(p$p3, "orange")[[1]]), c(1, 0, 0),
               tolerance = 1e-3)

  # Once centred, the reference plane lies in z = 0
  plane <- Filter(function(tr) tr$type == "mesh3d", p$p3$data)[[1]]
  expect_equal(unlist(plane$z), rep(0, 4))
})

test_that("SVD uses Norm_X/Y/Z when a specimen has fewer than three scars", {
  raw <- make_test_scar_data()
  two <- raw[raw$ID == "S002", ][1:2, ]
  out <- tempfile(fileext = ".html")
  export_alignment_html_svd(two, out)
  p <- read_panels(readLines(out, warn = FALSE), "S002")

  norm <- c(0.5, 0.5, 0.707)
  expect_equal(direction_of(lines_in(p$p0, "red")[[1]]),
               norm / sqrt(sum(norm^2)), tolerance = 1e-3)
})


# ---- export_alignment_html_lin2024 -------------------------------------------

test_that("export_alignment_html_lin2024 writes one page covering every specimen", {
  raw <- make_test_scar_data()
  out <- tempfile(fileext = ".html")

  res <- expect_invisible(export_alignment_html_lin2024(raw, out))
  expect_identical(res, out)

  html <- readLines(out, warn = FALSE)
  expect_true(any(grepl("Core Alignment Pipeline (Lin 2024)", html,
                        fixed = TRUE)))
  for (id in c("S001", "S002")) {
    option <- sprintf('<option value="%s">%s</option>', id, id)
    expect_true(any(grepl(option, html, fixed = TRUE)))
    expect_named(read_panels(html, id), c("p0", "p1", "p2"))
  }
})

test_that("Lin 2024 panels highlight the longest scar and anchor its start", {
  raw <- make_test_scar_data(n_scars = 5)
  out <- tempfile(fileext = ".html")
  export_alignment_html_lin2024(raw, out)
  p <- read_panels(readLines(out, warn = FALSE), "S002")

  longest <- max(get_scar_length(raw[raw$ID == "S002", ]))
  for (k in c("p0", "p1", "p2")) {
    pink <- lines_in(p[[k]], "pink")
    expect_length(pink, 1)
    expect_identical(pink[[1]]$line$width, 6L)
    expect_length(lines_in(p[[k]], "steelblue"), 4)
    # Rotation and translation keep the highlighted scar's length
    expect_equal(length_of(pink[[1]]), longest, tolerance = 1e-3)
  }

  # Step 1 aligns the morphological normal with Z
  expect_equal(direction_of(lines_in(p$p1, "red")[[1]]), c(0, 0, 1),
               tolerance = 1e-3)

  # Step 2 moves the longest scar's start point to x = y = 0
  start <- vapply(lines_in(p$p2, "pink")[[1]][c("x", "y")],
                  function(v) as.numeric(v[[1]]), numeric(1))
  expect_equal(unname(start), c(0, 0))
})


# ---- Panel layout ------------------------------------------------------------

test_that("every exported panel keeps its step title and true proportions", {
  raw <- make_test_scar_data()
  pages <- list(
    svd     = list(fun = export_alignment_html_svd,     n = 4),
    lin2024 = list(fun = export_alignment_html_lin2024, n = 3)
  )
  for (page in pages) {
    out <- tempfile(fileext = ".html")
    page$fun(raw, out)
    p <- read_panels(readLines(out, warn = FALSE), "S001")
    for (i in seq_len(page$n)) {
      layout <- p[[i]]$layout
      expect_match(layout$title$text, sprintf("^<b>Step %d", i - 1))
      expect_identical(layout$scene$aspectmode, "data")
      expect_identical(layout$scene$zaxis$title, "Z")
      expect_identical(layout$paper_bgcolor, "#f5f7fa")
    }
  }
})


# ---- 3-D helpers ---------------------------------------------------------------

test_that("add_tilted_plane_3d draws a square perpendicular to the normal", {
  center <- c(1, 2, 3)
  # The second normal is close to the X axis, which takes the other branch
  for (normal in list(c(0, 0, 1), c(1, 0, 0.1), c(0.3, -0.4, 0.8))) {
    fig <- add_tilted_plane_3d(plotly::plot_ly(), center, normal,
                               half_size = 2)
    tr  <- plotly::plotly_build(fig)$x$data[[1]]
    corners <- cbind(tr$x, tr$y, tr$z)

    n <- normal / sqrt(sum(normal^2))
    expect_equal(as.numeric(sweep(corners, 2, center) %*% n), rep(0, 4))
    sides <- sqrt(rowSums((corners - corners[c(2, 3, 4, 1), ])^2))
    expect_equal(sides, rep(4, 4))
  }
})

test_that("add_arrow_3d scales the arrow to the requested length", {
  fig <- add_arrow_3d(plotly::plot_ly(), origin = c(1, 1, 1),
                      direction = c(0, 0, 10), scale = 2)
  shaft <- plotly::plotly_build(fig)$x$data[[1]]
  expect_equal(as.numeric(shaft$x), c(1, 1))
  expect_equal(as.numeric(shaft$y), c(1, 1))
  expect_equal(as.numeric(shaft$z), c(1, 3))
})
