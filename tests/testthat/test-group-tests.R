# ==============================================================================
# test-group-tests.R
# Tests for compare_groups() and compare_indices()
# ==============================================================================

# Three groups of 8 cores; group B is more dispersed, group C is shifted.
make_groups <- function() {
  set.seed(1)
  X <- cbind(rnorm(24, rep(c(0, 0, 1), each = 8)), rnorm(24, sd = rep(c(1, 3, 1), each = 8)),
             rnorm(24))
  data.frame(group = rep(c("A", "B", "C"), each = 8), f1 = X[, 1], f2 = X[, 2], f3 = X[, 3])
}

# ---- compare_groups ---------------------------------------------------------------

test_that("compare_groups returns overall and pairwise tables", {
  x <- make_groups()
  res <- compare_groups(x, "group", permutations = 99)
  expect_s3_class(res, "spharmlithic_group_comparison")
  expect_equal(res$overall$test, c("PERMANOVA", "PERMDISP", "W*d"))
  expect_equal(nrow(res$pairwise), 3L)
  expect_true(all(c("PERMANOVA_p_adj", "PERMDISP_p_adj", "Wd_p_adj", "note") %in%
                    names(res$pairwise)))
  expect_equal(res$features, c("f1", "f2", "f3"))
  expect_equal(unname(res$groups), c(8L, 8L, 8L))
  expect_output(print(res), "Pairwise")
})

test_that("W*d statistic and p value agree with WdStar (reference values)", {
  # reference: WdStar 2.4.0, WdS(dist(X), f) and, after set.seed(42),
  # WdS.test(dist(X), f, nrep = 99)
  x <- make_groups()
  res <- compare_groups(x, "group", permutations = 99, seed = 42)
  expect_equal(res$overall$statistic[3], 1.1019225557363121, tolerance = 1e-14)
  expect_equal(res$overall$p[3], 0.42)
  ab <- res$pairwise[res$pairwise$group1 == "A" & res$pairwise$group2 == "B", ]
  expect_equal(ab$Wd, 0.020058076699376783, tolerance = 1e-14)
})

test_that("two-group W*d is |mean1 - mean2|^2 / (s1 / n1 + s2 / n2)", {
  x <- make_groups()
  x <- x[x$group != "C", ]
  X <- as.matrix(x[, -1])
  a <- X[x$group == "A", ]
  b <- X[x$group == "B", ]
  ref <- sum((colMeans(a) - colMeans(b))^2) /
    (sum(apply(a, 2, stats::var)) / 8 + sum(apply(b, 2, stats::var)) / 8)
  res <- compare_groups(x, "group", permutations = 9)
  expect_equal(res$overall$statistic[3], ref, tolerance = 1e-12)
  expect_null(res$pairwise)
})

test_that("compare_groups finds the features of spharm_spatiodirection", {
  d <- make_cores(n = c(8, 9, 10, 11, 12, 13), lean = c(0, 0, 0, 1, 1, 1))
  sp <- spharm_spatiodirection(d, quiet = TRUE)
  sp$grp <- rep(c("x", "y"), each = 3)
  res <- compare_groups(sp, "grp", permutations = 9)
  expect_equal(res$features, attr(sp, "features"))
  sp2 <- as.data.frame(sp)               # attribute lost: found by name
  attr(sp2, "features") <- NULL
  expect_equal(compare_groups(sp2, "grp", permutations = 9)$features, attr(sp, "features"))
})

test_that("compare_groups drops incomplete rows and single-core groups", {
  x <- make_groups()
  x$f1[1] <- NA
  x$group[24] <- "D"
  expect_message(res <- compare_groups(x, "group", permutations = 9),
                 "fewer than 2 cores: D")
  expect_equal(sum(res$groups), 22L)
  expect_error(compare_groups(x[1:3, ], "group", permutations = 9), "At least two groups")
})

test_that("compare_groups restores the caller's random stream", {
  x <- make_groups()
  set.seed(3)
  before <- .Random.seed
  compare_groups(x, "group", permutations = 9)
  expect_identical(.Random.seed, before)
})

# ---- compare_indices ----------------------------------------------------------------

test_that("compare_indices matches Welch tests and the SMD formula", {
  set.seed(4)
  d <- data.frame(g = rep(c("a", "b", "c"), c(8, 9, 10)), v = rnorm(27, rep(1:3, c(8, 9, 10))))
  res <- compare_indices(d, "v", "g", boot = 199)
  expect_s3_class(res, "spharmlithic_index_comparison")
  expect_equal(res$overall$test, "Welch ANOVA")
  expect_equal(res$overall$p, stats::oneway.test(v ~ g, d, var.equal = FALSE)$p.value)
  a <- d$v[d$g == "a"]; b <- d$v[d$g == "b"]
  ab <- res$pairwise[res$pairwise$group1 == "a" & res$pairwise$group2 == "b", ]
  expect_equal(ab$p, stats::t.test(a, b)$p.value)
  expect_equal(ab$smd, (mean(a) - mean(b)) / sqrt((stats::var(a) + stats::var(b)) / 2))
  expect_true(ab$smd_lower < ab$smd && ab$smd < ab$smd_upper)
  expect_equal(res$pairwise$p_adj, stats::p.adjust(res$pairwise$p, "holm"))
  expect_output(print(res), "standardised mean difference")
})

test_that("compare_indices uses Welch t for two groups and respects level order", {
  set.seed(5)
  d <- data.frame(g = factor(rep(c("b", "a"), each = 6), levels = c("b", "a")),
                  v = rnorm(12))
  res <- compare_indices(d, "v", "g", boot = 0)
  expect_equal(res$overall$test, "Welch t")
  expect_equal(res$pairwise$group1, "b")
  expect_true(is.na(res$pairwise$smd_lower))
})

test_that("compare_indices warns about raw columns and checks its input", {
  d <- data.frame(g = rep(c("a", "b"), each = 5), SPI_raw = 1:10 / 10)
  expect_message(compare_indices(d, "SPI_raw", "g", boot = 0), "rarefied")
  expect_error(compare_indices(d, "nope", "g"), "lacks column")
  expect_error(compare_indices(d, "SPI_raw", "nope"), "names no column")
})
