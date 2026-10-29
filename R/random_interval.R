# ==============================================================================
# random_interval.R
# The random interval of SPI, ER and IR at a scar count: the central 95% of
# the values a core of n scars takes when its directions are uniform on the
# sphere (no directional organisation).
#
#   SPI      the exact finite-n distribution of the resultant length of n unit
#            vectors (Rayleigh's random flight in three dimensions; Rayleigh
#            1919, Treloar 1946), for any n
#   ER, IR   5000 simulated random cores per n, n = 3..200
#
# For n = 3..200 all three are read from the table .random_interval_table
# (R/sysdata.rda, built by data-raw/random_interval.R from the tables of
# Xiao et al., in prep.); SPI above 200 is computed from its exact
# distribution; ER and IR above 200 are not available (NA).
# ==============================================================================

.RI_PARAMS <- c("SPI", "ER", "IR")
.RI_PROBS  <- c(lower = 0.025, upper = 0.975)

#' Random interval of SPI, ER or IR at a scar count
#'
#' Returns the central 95% interval of an index for a core whose scar
#' directions are uniform on the sphere, at each requested scar count. A raw
#' reading outside this interval departs from what random directions produce
#' at that scar count; a reading inside it cannot be told apart from random
#' directions at that scar count.
#'
#' The interval narrows as the scar count grows, so a raw value can only be
#' judged against the interval of its own scar count. [compute_spi()],
#' [compute_ei()] and [spi_ei_batch()] look the interval up for you.
#'
#' @param parameter Character vector: `"SPI"`, `"ER"` and/or `"IR"`
#'   (case-insensitive).
#' @param n Scar counts (whole numbers). Counts below 3 give `NA`.
#'
#' @return A tibble with one row per combination of `parameter` and `n`:
#'   `parameter`, `n`, `lower` (2.5th percentile) and `upper` (97.5th
#'   percentile).
#'
#' @details
#' SPI is computed from its exact distribution (the resultant length of n
#' unit vectors uniform on the sphere), ER and IR from 5000 simulated random
#' cores per scar count. Values for n = 3 to 200 are tabulated; SPI above 200
#' is computed on the fly; ER and IR above 200 are `NA`.
#'
#' Which side of the interval indicates organisation depends on the index:
#' SPI and ER above the upper bound, and IR below the lower bound, indicate
#' scar directions more concentrated than random. SPI below the lower bound
#' is typical of opposed (bidirectional) flaking.
#'
#' @references
#' Rayleigh, Lord (1919). On the problem of random vibrations, and of random
#' flights in one, two, or three dimensions. *Philosophical Magazine*,
#' **37**(220), 321--347.
#'
#' @examples
#' random_interval("SPI", n = c(6, 20))
#' random_interval(c("ER", "IR"), n = 10)
#'
#' @seealso [compute_spi()], [compute_ei()], [spi_ei_batch()],
#'   [plot_random_interval()]
#' @export
random_interval <- function(parameter = c("SPI", "ER", "IR"), n) {
  parameter <- .match_ri_param(parameter)
  if (missing(n) || !is.numeric(n) || !length(n)) {
    stop("`n` must be one or more scar counts.", call. = FALSE)
  }
  ok <- is.na(n) | n == round(n)
  if (!all(ok)) stop("`n` must be whole numbers (scar counts).", call. = FALSE)
  g <- expand.grid(n = n, parameter = parameter, stringsAsFactors = FALSE)
  b <- .ri_bounds(g$parameter, g$n)
  if (any(g$parameter %in% c("ER", "IR") & !is.na(g$n) & g$n > .ri_max_n())) {
    message("The random interval of ER and IR is tabulated up to n = ",
            .ri_max_n(), "; larger scar counts give NA.")
  }
  tibble::tibble(parameter = g$parameter, n = g$n, lower = b$lower,
                 upper = b$upper)
}

.match_ri_param <- function(parameter) {
  if (!is.character(parameter) || !length(parameter)) {
    stop("`parameter` must be \"SPI\", \"ER\" or \"IR\".", call. = FALSE)
  }
  i <- match(toupper(parameter), .RI_PARAMS)
  if (anyNA(i)) {
    stop("Unknown parameter(s): ", paste(parameter[is.na(i)], collapse = ", "),
         ". Use \"SPI\", \"ER\" or \"IR\".", call. = FALSE)
  }
  .RI_PARAMS[i]
}

.ri_max_n <- function() max(.random_interval_table$n)

# Lower and upper bound for each (parameter, n), vectorised.
.ri_bounds <- function(parameter, n) {
  parameter <- rep_len(parameter, length(n))
  tab <- .random_interval_table
  lower <- rep(NA_real_, length(n))
  upper <- rep(NA_real_, length(n))
  ok <- !is.na(n) & n >= 3
  hit <- match(paste(parameter, n), paste(tab$parameter, tab$n))
  intab <- ok & !is.na(hit)
  lower[intab] <- tab$lower[hit[intab]]
  upper[intab] <- tab$upper[hit[intab]]
  exact <- ok & is.na(hit) & parameter == "SPI"
  for (i in which(exact)) {
    q <- .spi_quantiles_cached(n[i])
    lower[i] <- q[["lower"]]
    upper[i] <- q[["upper"]]
  }
  list(lower = lower, upper = upper)
}

# Within or outside the interval; NA where either is missing.
.ri_position <- function(value, lower, upper) {
  out <- ifelse(value >= lower & value <= upper, "within", "outside")
  out[is.na(value) | is.na(lower) | is.na(upper)] <- NA_character_
  out
}

# ---- the exact distribution of SPI ------------------------------------------------
# SPI = |sum v_i| / n for n unit vectors uniform on S^2. With the characteristic
# function of one step phi(k) = sin(k)/k,
#
#     P(r <= R) = (2/pi) * int_0^Inf phi(k)^n * ( sin(kR)/k - R cos(kR) ) dk .
#
# The integrand oscillates, so it is integrated one pi-wide piece at a time up
# to a cut-off K set from the tail bound (1 + R) K^(1-n) / (n - 1). Near k = 0
# the bracket cancels to R^3 k^2 / 3 and the series is used instead.

.spi_integrand <- function(k, n, R) {
  out <- numeric(length(k))
  small <- k < 1e-4
  kk <- k[!small]
  out[!small] <- (sin(kk) / kk)^n * (sin(kk * R) / kk - R * cos(kk * R))
  out[small] <- R^3 * k[small]^2 / 3
  out
}

# P(SPI <= s) for n unit vectors
.spi_cdf <- function(s, n) {
  if (s <= 0) return(0)
  if (s >= 1) return(1)
  R <- n * s
  K <- max(20, ((1 + R) / ((n - 1) * 1e-8))^(1 / (n - 1)))
  edges <- unique(c(seq(0, K, by = pi), K))
  tot <- sum(vapply(seq_len(length(edges) - 1L), function(j)
    stats::integrate(.spi_integrand, edges[j], edges[j + 1L], n = n, R = R,
                     rel.tol = 1e-10, abs.tol = 1e-13,
                     subdivisions = 200L)$value, numeric(1)))
  min(1, max(0, 2 / pi * tot))
}

# the p-quantile of SPI at n
.spi_quantile <- function(p, n, lower = 1e-9, upper = 1 - 1e-9) {
  stats::uniroot(function(s) .spi_cdf(s, n) - p, c(lower, upper),
                 tol = 1e-10)$root
}

# several quantiles at one n, bracketing each root with the previous one
.spi_quantiles <- function(probs, n) {
  o <- order(probs)
  out <- numeric(length(probs))
  lo <- 1e-9
  for (i in o) {
    out[i] <- .spi_quantile(probs[i], n, lower = lo)
    lo <- max(1e-9, out[i] - 1e-6)
  }
  out
}

.spi_cache <- new.env(parent = emptyenv())

.spi_quantiles_cached <- function(n) {
  key <- as.character(n)
  if (is.null(.spi_cache[[key]])) {
    q <- .spi_quantiles(unname(.RI_PROBS), n)
    .spi_cache[[key]] <- stats::setNames(q, names(.RI_PROBS))
  }
  .spi_cache[[key]]
}
