# ==============================================================================
# test-plots.R
# Tests for plot_core(), export_outline_html(), plot_random_interval() and the
# outlines of the alignment viewers
# ==============================================================================

test_that("plot_core draws the scars, and the outline when given", {
  d <- make_disc(k = 6)
  p0 <- plotly::plotly_build(plot_core(d))
  p1 <- plotly::plotly_build(plot_core(d, outline = disc_outline()))
  expect_s3_class(plot_core(d), "plotly")
  # two traces per scar; with an outline also the outline, centre and plane
  expect_equal(length(p0$x$data), 12L)
  expect_equal(length(p1$x$data), 15L)
  expect_match(p1$x$layout$title$text, "start_dist")
})

test_that("plot_core picks one core by ID", {
  d <- rbind(make_disc(id = "a"), make_disc(k = 5, id = "b"))
  expect_error(plot_core(d), "choose one with `id`")
  p <- plotly::plotly_build(plot_core(d, id = "b"))
  expect_equal(length(p$x$data), 10L)
  expect_error(plot_core(d, id = "zz"), "No core with ID")
})

test_that("export_outline_html writes one page with every core that has an outline", {
  d <- rbind(make_disc(id = "a"), make_disc(id = "b"), make_disc(id = "c"))
  f <- tempfile(fileext = ".html")
  export_outline_html(d, list(a = disc_outline(), c = disc_outline()), f)
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, 'allPanels\\["a"\\]')
  expect_match(html, 'allPanels\\["c"\\]')
  expect_false(grepl('allPanels\\["b"\\]', html))
})

test_that("plot_random_interval returns a ggplot", {
  skip_if_not_installed("ggplot2")
  fab <- spi_ei_batch(make_cores(), rarefied = FALSE, keep = "Typology", quiet = TRUE)
  p <- plot_random_interval(fab, "SPI", group = "Typology")
  expect_s3_class(p, "ggplot")
  expect_error(plot_random_interval(fab[, "ID"], "SPI"), "n_scar")
})

test_that("the alignment viewers draw a moved outline", {
  d <- make_disc(k = 8, id = "S1")
  d$Norm_X <- 0; d$Norm_Y <- 0; d$Norm_Z <- 1
  d$Pos_X <- 0; d$Pos_Y <- 0; d$Pos_Z <- 0
  f <- tempfile(fileext = ".html")
  export_alignment_html_svd(d, f, outlines = list(S1 = disc_outline()))
  expect_match(paste(readLines(f, warn = FALSE), collapse = ""), "Black line")
  ps <- build_panel_scar("S1", d, outline = disc_outline())
  n_plain <- length(plotly::plotly_build(build_panel_scar("S1", d)$p3)$x$data)
  expect_equal(length(plotly::plotly_build(ps$p3)$x$data), n_plain + 2L)
  f2 <- tempfile(fileext = ".html")
  export_alignment_html_lin2024(d, f2, outlines = list(S1 = disc_outline()))
  expect_true(file.exists(f2))
})
