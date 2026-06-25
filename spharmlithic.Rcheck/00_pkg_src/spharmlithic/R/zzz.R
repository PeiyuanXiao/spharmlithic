# ==============================================================================
# zzz.R
# Package hooks:
#   .onLoad   — silent setup: lazily imports the bundled Python module.
#   .onAttach — user-facing startup message if the backend is unavailable.
# ==============================================================================

#' @noRd
sh_py <- NULL

.onLoad <- function(libname, pkgname) {
  py_path <- system.file("python", package = pkgname)
  
  # Attempt a delay-loaded import of the bundled Python module.
  # Wrapped in tryCatch so that .onLoad() always completes successfully:
  # if Python is already initialised (e.g. after a devtools::load_all()
  # reload in the same session), reticulate skips delay_load and attempts
  # an immediate import; that may fail if the path is not yet resolvable.
  # In that case sh_py stays NULL; the error surfaces only when a
  # Python-backed function is actually called.
  sh_py <<- tryCatch(
    reticulate::import_from_path(
      "spharmlithic_py",
      path = py_path,
      delay_load = list(
        environment = "r-spharmlithic",
        on_error    = function(e) NULL
      )
    ),
    error = function(e) NULL
  )
  
  invisible()
}

# Guard called by every Python-backed function. Turns the cryptic
# "attempt to apply non-function" (raised when sh_py is NULL because the
# backend was never set up) into an actionable message.
#' @noRd
.ensure_backend <- function() {
  if (is.null(sh_py)) {
    stop(
      "spharmlithic: Python backend not initialised.\n",
      "  Run install_spharmlithic_python() to set it up, or\n",
      "  use_spharmlithic_python('r-spharmlithic') to point at an existing env.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.onAttach <- function(libname, pkgname) {
  if (is.null(sh_py)) {
    packageStartupMessage(
      "spharmlithic: Python backend not initialised.\n",
      "  Run install_spharmlithic_python() to set it up,\n",
      "  or use_spharmlithic_python('my_env') to point at an existing env.\n",
      "  If you are reloading within an active session, call ",
      "use_spharmlithic_python() to restore the backend."
    )
  }
}

# Re-export the magrittr pipe so users don't need library(magrittr).
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' Replace zeros in compositional data
#'
#' Replaces zeros in each row of a compositional data matrix with a small
#' positive value and rescales the remaining non-zero entries so that each row
#' continues to sum to one.
#'
#' This function is intended for compositional data that contain structural
#' zeros which would otherwise prevent log-ratio transformations such as the
#' centred log-ratio (CLR) transform.
#'
#' By default, zeros are replaced with 65% of the smallest non-zero value in
#' each row.
#'
#' @param x A numeric matrix or data frame whose rows represent compositions
#'   that sum to one.
#' @param delta Optional numeric scalar giving the value used to replace zeros.
#'   If `NULL` (default), the replacement value is computed separately for each
#'   row as `fraction * min(non_zero_values)`.
#' @param fraction Numeric scalar in `(0, 1]` giving the fraction of the
#'   smallest non-zero value to use when calculating `delta`. Ignored if
#'   `delta` is supplied.
#'
#' @return A numeric matrix of the same dimensions as `x`, with zeros replaced
#'   and rows rescaled to sum to one.
#'
#' @export
replace_zeros <- function(x, delta = NULL, fraction = 0.65) {
  
  x <- as.matrix(x)
  
  if (!is.numeric(x)) {
    stop("`x` must be numeric.", call. = FALSE)
  }
  
  if (!is.numeric(fraction) ||
      length(fraction) != 1 ||
      fraction <= 0 ||
      fraction > 1) {
    stop("`fraction` must be a single number in (0, 1].", call. = FALSE)
  }
  
  for (i in seq_len(nrow(x))) {
    
    row_i <- x[i, ]
    zero_idx <- row_i == 0
    
    if (!any(zero_idx)) {
      next
    }
    
    if (all(zero_idx)) {
      warning(
        sprintf("Row %d contains only zeros and was left unchanged.", i),
        call. = FALSE
      )
      next
    }
    
    d <- if (is.null(delta)) {
      min(row_i[!zero_idx]) * fraction
    } else {
      delta
    }
    
    n_zero <- sum(zero_idx)
    
    if (n_zero * d >= 1) {
      stop(
        sprintf(
          "Replacement value is too large in row %d: n_zero * delta >= 1.",
          i
        ),
        call. = FALSE
      )
    }
    
    row_i[zero_idx] <- d
    row_i[!zero_idx] <- row_i[!zero_idx] * (1 - n_zero * d)
    
    x[i, ] <- row_i
    
  }
  
  x
}

#' Compute isometric log-ratio coordinates
#'
#' Transform compositional data to isometric log-ratio (ILR) coordinates.
#'
#' Columns with zero variance are removed before transformation because they
#' do not contribute information and can cause numerical problems.
#'
#' Zeros are replaced using [replace_zeros()] prior to transformation.
#'
#' @param x A numeric matrix or data frame of compositional data.
#' @param delta Optional replacement value passed to [replace_zeros()].
#' @param fraction Fraction of the smallest non-zero value used by
#'   [replace_zeros()] when `delta = NULL`.
#'
#' @return A data frame containing ILR coordinates.
#'
#' @examples
#' x <- data.frame(
#'   a = c(0.5, 0.4, 0.6),
#'   b = c(0.3, 0.4, 0.2),
#'   c = c(0.2, 0.2, 0.2)
#' )
#'
#' make_ilr(x)
#'
#' @export
make_ilr <- function(x, delta = NULL, fraction = 0.65) {
  
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be numeric.", call. = FALSE)
  }
  
  keep <- apply(
    x,
    2,
    function(v) stats::sd(v, na.rm = TRUE) > 0
  )
  
  if (sum(keep) < 2) {
    stop(
      "At least two non-constant columns are required for ILR transformation.",
      call. = FALSE
    )
  }
  
  x <- x[, keep, drop = FALSE]
  
  ilr_mat <- compositions::ilr(
    replace_zeros(
      x,
      delta = delta,
      fraction = fraction
    )
  )
  
  ilr_df <- as.data.frame(ilr_mat)
  colnames(ilr_df) <- paste0("ilr_", seq_len(ncol(ilr_df)))
  rownames(ilr_df) <- rownames(x)
  ilr_df
}

#' Summarise power spectrum diagnostics by degree
#'
#' Calculate per-degree coefficients of variation and cumulative mean power
#' from a power spectrum.
#'
#' The coefficient of variation (CV) quantifies across-specimen variability
#' for each harmonic degree. The cumulative power indicates the proportion of
#' total mean power explained by successive degrees.
#'
#' @param power_df A data frame containing columns named
#'   `"power_l1"`, `"power_l2"`, ..., `"power_lN"`.
#' @param descriptor Character string identifying the descriptor type
#'   (e.g., `"shape"`, `"amplitude"`, `"power"`).
#' @param max_degree Maximum harmonic degree to analyse.
#'
#' @return A data frame with one row per degree containing:
#' \describe{
#'   \item{descriptor}{Descriptor name supplied by the user.}
#'   \item{degree}{Harmonic degree.}
#'   \item{mean_power}{Mean power at that degree.}
#'   \item{cv_pct}{Coefficient of variation (%).}
#'   \item{cumul_pct}{Cumulative percentage of total mean power.}
#' }
#'
#' @examples
#' x <- data.frame(
#'   power_l1 = c(0.4, 0.5, 0.6),
#'   power_l2 = c(0.3, 0.3, 0.2),
#'   power_l3 = c(0.3, 0.2, 0.2)
#' )
#'
#' degree_diagnostics(x, "example", max_degree = 3)
#'
#' @export
degree_diagnostics <- function(power_df, descriptor, max_degree = 20) {
  cols <- paste0("power_l", seq_len(max_degree))
  
  missing_cols <- setdiff(cols, names(power_df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }
  
  mat <- as.matrix(power_df[, cols, drop = FALSE])
  
  if (!is.numeric(mat)) {
    stop("Power spectrum columns must be numeric.", call. = FALSE)
  }
  
  mu <- colMeans(mat, na.rm = TRUE)
  
  cv <- vapply(seq_along(mu), function(i) {
    if (mu[i] == 0) {
      NA_real_
    } else {
      stats::sd(mat[, i], na.rm = TRUE) / mu[i] * 100
    }
  }, numeric(1))
  
  data.frame(
    descriptor = rep(descriptor, max_degree),
    degree = seq_len(max_degree),
    mean_power = mu,
    cv_pct = cv,
    cumul_pct = cumsum(mu) / sum(mu) * 100,
    row.names = NULL
  )
}

#' Pairwise PERMANOVA comparisons
#'
#' Perform pairwise PERMANOVA comparisons among all levels of a grouping
#' variable using [vegan::adonis2()].
#'
#' For each pair of groups, a distance matrix is calculated from the supplied
#' multivariate data and tested using PERMANOVA. P-values are adjusted for
#' multiple comparisons.
#' 
#' @details
#' This function assumes that compositional data have already been transformed
#' to ILR coordinates. Euclidean distances computed on ILR coordinates are
#' equivalent to Aitchison distances in the original compositional space.
#'
#' @param x A numeric matrix or data frame containing observations in rows and
#'   variables in columns (e.g., ILR coordinates).
#' @param group A grouping variable with one value per row of `x`.
#' @param permutations Number of permutations passed to
#'   [vegan::adonis2()]. Default is 999.
#' @param p_adjust_method Method used by [stats::p.adjust()].
#'   Default is `"holm"`.
#' @param distance Distance metric passed to [stats::dist()].
#'   Default is `"euclidean"`.
#'
#' @return A data frame with one row per pairwise comparison containing:
#' \describe{
#'   \item{group1}{First group.}
#'   \item{group2}{Second group.}
#'   \item{pair}{Comparison label.}
#'   \item{n1}{Sample size of group1.}
#'   \item{n2}{Sample size of group2.}
#'   \item{R2}{PERMANOVA R-squared.}
#'   \item{F}{Pseudo-F statistic.}
#'   \item{p}{Permutation p-value.}
#'   \item{p_adj}{Adjusted p-value.}
#' }
#'
#' @examples
#' set.seed(123)
#' 
#' x <- matrix(rnorm(45), ncol = 3)
#' grp <- rep(c("A", "B", "C"), each = 5)
#'
#' pairwise_permanova(x, grp, permutations = 99)
#'
#' @export
pairwise_permanova <- function(
    x,
    group,
    permutations = 999,
    p_adjust_method = "holm",
    distance = "euclidean"
) {
  
  x <- as.matrix(x)
  
  if (!is.numeric(x)) {
    stop("`x` must be numeric.", call. = FALSE)
  }
  
  if (nrow(x) != length(group)) {
    stop(
      "Number of rows in `x` must equal length of `group`.",
      call. = FALSE
    )
  }
  
  group <- droplevels(factor(group))
  
  if (nlevels(group) < 2) {
    stop(
      "At least two groups are required.",
      call. = FALSE
    )
  }
  
  combs <- utils::combn(levels(group), 2, simplify = TRUE)
  
  results <- lapply(seq_len(ncol(combs)), function(i) {
    
    g1 <- combs[1, i]
    g2 <- combs[2, i]
    
    keep <- group %in% c(g1, g2)
    
    sub_x <- x[keep, , drop = FALSE]
    sub_group <- droplevels(group[keep])
    
    n1 <- sum(sub_group == g1)
    n2 <- sum(sub_group == g2)
    
    if (n1 < 2 || n2 < 2) {
      return(data.frame(
        group1 = g1,
        group2 = g2,
        pair = paste(g1, "vs", g2),
        n1 = n1,
        n2 = n2,
        R2 = NA_real_,
        F = NA_real_,
        p = NA_real_
      ))
    }
    
    d <- stats::dist(sub_x, method = distance)
    
    fit <- vegan::adonis2(
      d ~ sub_group,
      data = data.frame(sub_group),
      permutations = permutations
    )
    
    data.frame(
      group1 = g1,
      group2 = g2,
      pair = paste(g1, "vs", g2),
      n1 = n1,
      n2 = n2,
      R2 = fit$R2[1],
      F = fit$F[1],
      p = fit[["Pr(>F)"]][1]
    )
  })
  
  out <- do.call(rbind, results)
  
  out$p_adj <- stats::p.adjust(
    out$p,
    method = p_adjust_method
  )
  
  rownames(out) <- NULL
  
  out
}
