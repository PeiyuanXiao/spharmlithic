# ==============================================================================
# group_tests.R
# Comparing groups of cores:
#   - compare_groups()  : whole scar patterns (many features at once):
#                         PERMANOVA, PERMDISP and the W*d test
#   - compare_indices() : one index at a time: Welch tests and standardised
#                         mean differences
# ==============================================================================

#' Compare the overall scar patterns of groups of cores
#'
#' Tests whether groups of cores differ in their scar patterns as a whole,
#' using a set of features such as those of [spharm_spatiodirection()]. Three
#' permutation tests on Euclidean distances are reported side by side, for
#' all groups together and for every pair of groups:
#'
#' * **PERMANOVA** (Anderson 2001): do the group centroids differ?
#' * **PERMDISP** (Anderson 2006): do the groups differ in dispersion, the
#'   spread of their cores around their centroid?
#' * the \eqn{W_d^*} test (Hamidi et al. 2019): a Welch-type test of the
#'   centroids that allows the groups to differ in dispersion.
#'
#' When PERMDISP is significant, PERMANOVA may mistake a difference in
#' dispersion for a difference in location, especially with unequal group
#' sizes; the location difference should then be judged by \eqn{W_d^*}. The
#' printed result points out where this applies.
#'
#' @param x A data frame or numeric matrix with one row per core, such as
#'   the output of [spharm_spatiodirection()].
#' @param group The groups: the name of a column of `x`, or a vector with one
#'   value per row of `x`.
#' @param features Names of the columns of `x` to compare. Default: the
#'   features recorded by [spharm_spatiodirection()] (attribute
#'   `"features"`), else the columns named like its features (`dir_1`,
#'   `start_2`, ...), else all numeric columns except `n_scar`.
#' @param permutations Number of permutations of each test. Default 999.
#' @param alpha Significance level used for the notes. Default 0.05.
#' @param p_adjust Method for adjusting the pairwise p values within each
#'   test, passed to [stats::p.adjust()]. Default `"holm"`.
#' @param permdisp_type Distances to the group `"centroid"` (default, the
#'   original PERMDISP) or to the spatial `"median"`; passed to
#'   [vegan::betadisper()].
#' @param seed Seed set before each test, so each is reproducible on its own;
#'   the caller's random stream is restored afterwards. `NULL` uses the
#'   current stream. Default 2226.
#'
#' @return An object of class `spharmlithic_group_comparison`, a list with
#'   \describe{
#'     \item{overall}{A tibble with one row per test: `test`, `statistic`
#'       (pseudo-F for PERMANOVA and PERMDISP, \eqn{W_d^*}), `R2` (PERMANOVA) and
#'       `p`.}
#'     \item{pairwise}{With three or more groups, a tibble with one row per
#'       pair: `group1`, `group2`, `n1`, `n2`, the statistic and p value of
#'       each test, the adjusted p values (`*_p_adj`) and `note`; `NULL` with
#'       two groups.}
#'     \item{groups}{The number of cores per group.}
#'     \item{features}{The columns compared.}
#'     \item{settings}{The arguments used.}
#'   }
#'   Printing the object shows both tables and the notes.
#'
#' @details
#' Rows with a missing feature or group are left out, and so are groups with
#' fewer than two cores; a message says which. The features of
#' [spharm_spatiodirection()] do not change with the scar count, so the
#' scar count should not be added as a covariate.
#'
#' The \eqn{W_d^*} statistic is computed in the package and agrees with
#' `WdStar::WdS()` (Hamidi et al. 2019); its p value is the proportion of
#' permutations of the group labels, the observed one included, with a
#' statistic at least as large.
#'
#' @references
#' Anderson, M. J. (2001). A new method for non-parametric multivariate
#' analysis of variance. \emph{Austral Ecology}, \strong{26}(1), 32--46.
#'
#' Anderson, M. J. (2006). Distance-based tests for homogeneity of
#' multivariate dispersions. \emph{Biometrics}, \strong{62}(1), 245--253.
#'
#' Hamidi, B., Wallace, K., Vasu, C., & Alekseyenko, A. V. (2019).
#' \eqn{W_d^*}-test: robust distance-based multivariate analysis of variance.
#' \emph{Microbiome}, \strong{7}, 51.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics}, \strong{6}(2), 65--70.
#'
#' @examples
#' set.seed(1)
#' x <- data.frame(group = rep(c("A", "B", "C"), each = 8),
#'                 f1 = rnorm(24, rep(c(0, 0, 1), each = 8)),
#'                 f2 = rnorm(24, sd = rep(c(1, 3, 1), each = 8)))
#' compare_groups(x, "group", permutations = 199)
#'
#' @seealso [spharm_spatiodirection()], [compare_indices()],
#'   [pairwise_permanova()]
#' @export
compare_groups <- function(x, group, features = NULL, permutations = 999,
                           alpha = 0.05, p_adjust = "holm",
                           permdisp_type = c("centroid", "median"), seed = 2226) {
  permdisp_type <- match.arg(permdisp_type)
  .check_seed(seed)
  .check_test_args(permutations, alpha, p_adjust)
  g <- .group_values(x, group)
  features <- .feature_cols(x, features)
  X <- as.matrix(as.data.frame(x)[, features, drop = FALSE])
  storage.mode(X) <- "double"

  ok <- stats::complete.cases(X) & !is.na(g)
  notes <- character(0)
  if (any(!ok)) {
    notes <- c(notes, sprintf("Left out %d core(s) with a missing feature or group.",
                              sum(!ok)))
  }
  X <- X[ok, , drop = FALSE]
  g <- g[ok]
  tab <- table(g)
  small <- names(tab)[tab < 2L]
  if (length(small)) {
    notes <- c(notes, sprintf("Left out %d group(s) with fewer than 2 cores: %s.",
                              length(small), paste(small, collapse = ", ")))
    keep <- !(as.character(g) %in% small)
    X <- X[keep, , drop = FALSE]
    g <- g[keep]
  }
  g <- droplevels(g)
  if (nlevels(g) < 2L) {
    .say(notes)
    stop("At least two groups with two or more cores each are needed.", call. = FALSE)
  }
  .say(notes)

  overall <- .three_tests(X, g, permutations, permdisp_type, seed)
  pairwise <- NULL
  if (nlevels(g) > 2L) {
    prs <- utils::combn(levels(g), 2)
    pairwise <- do.call(rbind, lapply(seq_len(ncol(prs)), function(k) {
      sel <- g %in% prs[, k]
      r <- .three_tests(X[sel, , drop = FALSE], droplevels(g[sel]), permutations,
                        permdisp_type, seed)
      data.frame(group1 = prs[1, k], group2 = prs[2, k],
                 n1 = sum(g == prs[1, k]), n2 = sum(g == prs[2, k]),
                 PERMANOVA_F = r$statistic[1], PERMANOVA_R2 = r$R2[1],
                 PERMANOVA_p = r$p[1], PERMDISP_F = r$statistic[2],
                 PERMDISP_p = r$p[2], Wd = r$statistic[3], Wd_p = r$p[3],
                 stringsAsFactors = FALSE)
    }))
    for (t in c("PERMANOVA", "PERMDISP", "Wd")) {
      pairwise[[paste0(t, "_p_adj")]] <- stats::p.adjust(pairwise[[paste0(t, "_p")]],
                                                         method = p_adjust)
    }
    pairwise$note <- .disp_note(pairwise$PERMANOVA_p_adj, pairwise$PERMDISP_p_adj, alpha)
    pairwise <- tibble::as_tibble(pairwise)
  }

  structure(list(overall = tibble::as_tibble(overall), pairwise = pairwise,
                 groups = c(tab[levels(g)]), features = features,
                 settings = list(permutations = permutations, alpha = alpha,
                                 p_adjust = p_adjust, permdisp_type = permdisp_type,
                                 seed = seed, distance = "euclidean")),
            class = "spharmlithic_group_comparison")
}

#' @export
print.spharmlithic_group_comparison <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Group comparison of", length(x$features), "features",
      sprintf("(Euclidean distance, %d permutations)\n", s$permutations))
  cat("Groups:", paste(sprintf("%s (n = %d)", names(x$groups), x$groups),
                       collapse = ", "), "\n\n")
  cat("Overall\n")
  print(as.data.frame(x$overall), digits = digits, row.names = FALSE)
  notes <- character(0)
  o <- x$overall
  if (o$p[2] <= s$alpha) {
    notes <- c(notes, paste0(
      "PERMDISP is significant overall: the groups differ in dispersion, so ",
      "PERMANOVA may confuse dispersion with location; judge the location ",
      "difference by W*d."))
  }
  if (!is.null(x$pairwise)) {
    cat(sprintf("\nPairwise (p values adjusted by %s within each test)\n", s$p_adjust))
    pw <- as.data.frame(x$pairwise)
    show <- pw[, c("group1", "group2", "n1", "n2", "PERMANOVA_p_adj",
                   "PERMDISP_p_adj", "Wd_p_adj", "note")]
    print(show, digits = digits, row.names = FALSE)
    k <- sum(nzchar(pw$note))
    if (k) {
      notes <- c(notes, sprintf(paste0(
        "PERMDISP is significant for %d of %d pairs: for these pairs PERMANOVA ",
        "may confuse dispersion with location; judge them by W*d."), k, nrow(pw)))
    }
  }
  if (length(notes)) {
    cat("\n", paste(unlist(lapply(notes, strwrap, initial = "Note: ", prefix = "      ")),
                    collapse = "\n"), "\n", sep = "")
  }
  invisible(x)
}

#' Compare groups of cores on single indices
#'
#' Compares groups of cores on one index at a time, such as the rarefied SPI,
#' ER and IR of [spi_ei_batch()] or the CPI and PRA summaries of
#' [cpi_pra_batch()]. Each index is tested for a difference in group means
#' with Welch-type tests, which allow the groups to have unequal variances,
#' and every pair of groups gets a standardised mean difference with a
#' bootstrap confidence interval.
#'
#' @param data A data frame with one row per core.
#' @param cols Names of the index columns to compare.
#' @param group The groups: the name of a column of `data`, or a vector with
#'   one value per row. The order of the factor levels (or of the sorted
#'   values) sets which group is `group1` in each pair.
#' @param alpha Significance level used in the notes. Default 0.05.
#' @param p_adjust Method for adjusting the pairwise p values within each
#'   index, passed to [stats::p.adjust()]. Default `"holm"`.
#' @param boot Number of bootstrap resamples for the confidence interval of
#'   the standardised mean difference; 0 for none. Default 9999.
#' @param seed Seed set before the bootstrap of each pair; the caller's random
#'   stream is restored afterwards. `NULL` uses the current stream.
#'   Default 2226.
#'
#' @return An object of class `spharmlithic_index_comparison`, a list with
#'   \describe{
#'     \item{overall}{A tibble with one row per index: `index`, `test`
#'       (`"Welch t"` for two groups, `"Welch ANOVA"` for more), `statistic`,
#'       `df1`, `df2` (`NA` for the t test), `p`, `n_groups` and `n`.}
#'     \item{pairwise}{A tibble with one row per index and pair of groups:
#'       `index`, `group1`, `group2`, `n1`, `n2`, `mean1`, `mean2`, `smd`
#'       (standardised mean difference, group1 minus group2), `smd_lower`,
#'       `smd_upper` (95% bootstrap interval), `p` (Welch t test) and
#'       `p_adj`.}
#'     \item{settings}{The arguments used.}
#'   }
#'
#' @details
#' Two groups are compared with Welch's t test, three or more with Welch's
#' one-way ANOVA ([stats::oneway.test()]) followed by Welch t tests between
#' every pair, adjusted by `p_adjust` within each index (not across indices,
#' each index answering its own question). Welch-type tests are used because
#' differences in scar count can leave groups with unequal variances.
#'
#' The standardised mean difference is the difference between the group
#' means divided by the square root of the mean of the two group variances;
#' its absolute value is the effect size reported by Xiao et al. (in prep.).
#' The 95% interval is a percentile bootstrap stratified by group: each
#' group is resampled with replacement within itself.
#'
#' **Which columns.** Raw SPI, ER and IR shift with the scar count; compare
#' groups with the rarefied values (`*_rarefied`). A message is given when a
#' `*_raw` column is passed. The mean and GMD of CPI and PRA do not need
#' rarefying.
#'
#' Missing values are left out index by index, and so are groups with fewer
#' than two cores.
#'
#' @references
#' Welch, B. L. (1951). On the comparison of several mean values: an
#' alternative approach. \emph{Biometrika}, \strong{38}(3/4), 330--336.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics}, \strong{6}(2), 65--70.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(Typology = rep(c("A", "B", "C"), each = 10),
#'                 SPI_rarefied = rnorm(30, rep(c(0.5, 0.6, 0.6), each = 10), 0.08),
#'                 IR_rarefied = rnorm(30, 0.2, 0.05))
#' compare_indices(d, c("SPI_rarefied", "IR_rarefied"), "Typology", boot = 999)
#'
#' @seealso [spi_ei_batch()], [cpi_pra_batch()], [compare_groups()]
#' @export
compare_indices <- function(data, cols, group, alpha = 0.05, p_adjust = "holm",
                            boot = 9999, seed = 2226) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.character(cols) || !length(cols)) {
    stop("`cols` must name one or more index columns.", call. = FALSE)
  }
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("`data` lacks column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  for (cc in cols) {
    if (!is.numeric(data[[cc]])) stop("Column `", cc, "` must be numeric.", call. = FALSE)
  }
  .check_test_args(1, alpha, p_adjust)
  if (length(boot) != 1L || !is.numeric(boot) || is.na(boot) || boot < 0 ||
      boot != round(boot)) {
    stop("`boot` must be a whole number >= 0.", call. = FALSE)
  }
  .check_seed(seed)
  g_all <- .group_values(data, group)

  notes <- character(0)
  raw_cols <- cols[grepl("_raw$", cols)]
  if (length(raw_cols)) {
    notes <- c(notes, sprintf(paste("%s: raw values shift with the scar count;",
                                    "compare groups with the rarefied values",
                                    "(*_rarefied) instead."),
                              paste(raw_cols, collapse = ", ")))
  }
  overall <- list()
  pairwise <- list()
  for (cc in cols) {
    v <- data[[cc]]
    ok <- !is.na(v) & !is.na(g_all)
    v <- v[ok]
    g <- g_all[ok]
    tab <- table(g)
    small <- names(tab)[tab < 2L]
    if (length(small)) {
      notes <- c(notes, sprintf("%s: left out group(s) with fewer than 2 cores: %s.",
                                cc, paste(small, collapse = ", ")))
      v <- v[!(as.character(g) %in% small)]
      g <- g[!(as.character(g) %in% small)]
    }
    g <- droplevels(g)
    k <- nlevels(g)
    if (k < 2L) {
      notes <- c(notes, sprintf("%s: fewer than two groups with two or more values; not tested.", cc))
      next
    }
    if (k == 2L) {
      tt <- stats::t.test(v[g == levels(g)[1]], v[g == levels(g)[2]], var.equal = FALSE)
      overall[[cc]] <- data.frame(index = cc, test = "Welch t",
                                  statistic = unname(tt$statistic),
                                  df1 = unname(tt$parameter), df2 = NA_real_,
                                  p = tt$p.value, n_groups = k, n = length(v),
                                  stringsAsFactors = FALSE)
    } else {
      ow <- stats::oneway.test(v ~ g, var.equal = FALSE)
      overall[[cc]] <- data.frame(index = cc, test = "Welch ANOVA",
                                  statistic = unname(ow$statistic),
                                  df1 = unname(ow$parameter[1]),
                                  df2 = unname(ow$parameter[2]),
                                  p = ow$p.value, n_groups = k, n = length(v),
                                  stringsAsFactors = FALSE)
    }
    prs <- utils::combn(levels(g), 2)
    pw <- do.call(rbind, lapply(seq_len(ncol(prs)), function(j) {
      a <- v[g == prs[1, j]]
      b <- v[g == prs[2, j]]
      tt <- stats::t.test(a, b, var.equal = FALSE)
      ci <- .smd_ci(a, b, boot, seed)
      data.frame(index = cc, group1 = prs[1, j], group2 = prs[2, j],
                 n1 = length(a), n2 = length(b), mean1 = mean(a), mean2 = mean(b),
                 smd = .smd(a, b), smd_lower = ci[1], smd_upper = ci[2],
                 p = tt$p.value, stringsAsFactors = FALSE)
    }))
    pw$p_adj <- stats::p.adjust(pw$p, method = p_adjust)
    pairwise[[cc]] <- pw
  }
  .say(notes)
  if (!length(overall)) stop("No index could be tested.", call. = FALSE)
  structure(list(overall = tibble::as_tibble(do.call(rbind, unname(overall))),
                 pairwise = tibble::as_tibble(do.call(rbind, unname(pairwise))),
                 settings = list(alpha = alpha, p_adjust = p_adjust, boot = boot,
                                 seed = seed)),
            class = "spharmlithic_index_comparison")
}

#' @export
print.spharmlithic_index_comparison <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Index comparison (Welch tests; pairwise p adjusted by", s$p_adjust,
      "within each index)\n\n")
  print(as.data.frame(x$overall), digits = digits, row.names = FALSE)
  cat("\nPairwise: standardised mean difference (group1 - group2)",
      if (s$boot > 0) sprintf("with 95%% bootstrap interval (%d resamples)", s$boot), "\n")
  pw <- as.data.frame(x$pairwise)
  print(pw[, c("index", "group1", "group2", "n1", "n2", "smd", "smd_lower",
               "smd_upper", "p", "p_adj")], digits = digits, row.names = FALSE)
  invisible(x)
}

# ---- helpers ---------------------------------------------------------------------

.check_test_args <- function(permutations, alpha, p_adjust) {
  if (length(permutations) != 1L || !is.numeric(permutations) || is.na(permutations) ||
      permutations < 1 || permutations != round(permutations)) {
    stop("`permutations` must be a whole number >= 1.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !p_adjust %in% stats::p.adjust.methods) {
    stop("`p_adjust` must be one of: ", paste(stats::p.adjust.methods, collapse = ", "),
         ".", call. = FALSE)
  }
  invisible(TRUE)
}

# The group of every row of x, as a factor.
.group_values <- function(x, group) {
  nr <- NROW(x)
  if (is.character(group) && length(group) == 1L && nr != 1L) {
    if (!is.data.frame(x) || !group %in% names(x)) {
      stop("`group` names no column of the data: ", group, ".", call. = FALSE)
    }
    group <- x[[group]]
  }
  if (length(group) != nr) {
    stop("`group` must name a column of the data or give one value per row.",
         call. = FALSE)
  }
  if (is.factor(group)) group else factor(group)
}

# The feature columns of x.
.feature_cols <- function(x, features) {
  nm <- colnames(x)
  if (is.null(nm)) stop("`x` must have column names.", call. = FALSE)
  if (is.null(features)) {
    a <- attr(x, "features")
    features <- if (!is.null(a) && all(a %in% nm)) a else {
      sd <- nm[grepl("^(dir|start|same|diff)_[0-9]+$", nm)]
      if (length(sd)) sd else {
        num <- if (is.data.frame(x)) nm[vapply(x, is.numeric, TRUE)] else nm
        setdiff(num, c("n_scar", "n_cpi"))
      }
    }
  }
  miss <- setdiff(features, nm)
  if (length(miss)) {
    stop("`x` lacks feature column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  if (!length(features)) stop("No feature columns to compare.", call. = FALSE)
  features
}

# PERMANOVA, PERMDISP and W*d of the rows of X between the levels of g.
.three_tests <- function(X, g, permutations, permdisp_type, seed) {
  d <- stats::dist(X)
  grp <- g
  pa <- .with_seed(seed, vegan::adonis2(d ~ grp, data = data.frame(grp = grp),
                                        permutations = permutations))
  bd <- vegan::betadisper(d, grp, type = permdisp_type)
  pd <- .with_seed(seed, vegan::permutest(bd, permutations = permutations))
  wd <- .with_seed(seed, .wd_test(d, grp, permutations))
  data.frame(test = c("PERMANOVA", "PERMDISP", "W*d"),
             statistic = c(pa$F[1], pd$tab$F[1], wd$statistic),
             R2 = c(pa$R2[1], NA_real_, NA_real_),
             p = c(pa[["Pr(>F)"]][1], pd$tab[["Pr(>F)"]][1], wd$p),
             stringsAsFactors = FALSE)
}

.disp_note <- function(p_perm, p_disp, alpha) {
  ifelse(is.na(p_disp) | p_disp > alpha, "",
         ifelse(!is.na(p_perm) & p_perm <= alpha,
                "PERMANOVA and PERMDISP significant: read W*d",
                "dispersion differs: read W*d"))
}

# ---- the W*d test (Hamidi et al. 2019) ---------------------------------------------
# For groups g = 1..k with n_g members, from the squared distances D2:
#   SS_gh  sums of squared distances within and between groups (each pair once)
#   s_g    = SS_gg / (n_g (n_g - 1))         the dispersion of group g
#   W      = sum_g n_g / s_g
#   Ws     = sum_{g < h} (n_g + n_h) / (s_g s_h) *
#              [ (SS_gg + SS_hh + 2 SS_gh) / (n_g + n_h) - SS_gg / n_g - SS_hh / n_h ]
#   h      = sum_g (1 - n_g / s_g / W)^2 / (n_g - 1)
#   W*d    = Ws / W / (k - 1) / (1 + 2 (k - 2) / (k^2 - 1) * h)
# With two groups W*d = |mean_1 - mean_2|^2 / (s_1 / n_1 + s_2 / n_2).

.wd_stat <- function(D2, f) {
  K <- vapply(levels(f), function(lev) f == lev, logical(length(f)))
  storage.mode(K) <- "double"
  SS <- t(K) %*% D2 %*% K / 2
  ns <- colSums(K)
  s2 <- diag(SS) / ns / (ns - 1)
  W <- sum(ns / s2)
  k <- length(ns)
  prs <- utils::combn(k, 2)
  Ws <- sum(vapply(seq_len(ncol(prs)), function(j) {
    idx <- prs[, j]
    sum(ns[idx]) / prod(s2[idx]) *
      (sum(SS[idx, idx]) / sum(ns[idx]) - sum(diag(SS[idx, idx]) / ns[idx]))
  }, 0))
  h <- sum((1 - ns / s2 / W)^2 / (ns - 1))
  Ws / W / (k - 1) / (1 + (2 * (k - 2) / (k^2 - 1)) * h)
}

# Permutation test: the labels are permuted `nrep` times; the p value counts
# the observed statistic among the permutations.
.wd_test <- function(d, f, nrep) {
  D2 <- as.matrix(d)^2
  N <- length(f)
  st <- c(.wd_stat(D2, f), replicate(nrep, .wd_stat(D2, f[sample(N)])))
  list(statistic = st[1], p = sum(st >= st[1]) / (nrep + 1))
}

# ---- standardised mean difference ---------------------------------------------------

.smd <- function(a, b) (mean(a) - mean(b)) / sqrt((stats::var(a) + stats::var(b)) / 2)

# 95% percentile interval of .smd() from a bootstrap stratified by group.
.smd_ci <- function(a, b, boot, seed) {
  if (boot < 1) return(c(NA_real_, NA_real_))
  bs <- .with_seed(seed, vapply(seq_len(boot), function(i) {
    .smd(a[sample.int(length(a), replace = TRUE)],
         b[sample.int(length(b), replace = TRUE)])
  }, 0))
  bs <- bs[is.finite(bs)]
  if (!length(bs)) return(c(NA_real_, NA_real_))
  stats::quantile(bs, c(0.025, 0.975), names = FALSE)
}
