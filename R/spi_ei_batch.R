# ==============================================================================
# spi_ei_batch.R
# SPI, ER and IR of every core in a scar table: raw, rarefied, and the random
# interval at each core's scar count.
# ==============================================================================

#' SPI, ER and IR for every core
#'
#' Computes the Scar Pattern Index (SPI), the elongation ratio (ER) and the
#' isotropy ratio (IR) of every core in a scar table, each as a raw value, a
#' rarefied value and a random interval at the core's scar count. This is the
#' batch form of [compute_spi()] and [compute_ei()].
#'
#' @param data A data frame with one row per scar and the columns named by
#'   `id_col`, `start_cols` and `end_cols`. Scars of all cores are stacked.
#' @param id_col Name of the core ID column. Default `"ID"`.
#' @param start_cols,end_cols Names of the start- and end-point coordinate
#'   columns (X, Y, Z). Defaults `Start_X/Y/Z` and `End_X/Y/Z`, the columns
#'   the scars are recorded in. The aligned columns of [align_scar_batch()]
#'   (`s_x`... and `e_x`...) give the same results: these indices do not
#'   change when a core is rotated or moved.
#' @param min_length Scars shorter than this (in the units of the
#'   coordinates) are left out. Default 0 (keep all). Use it to repeat an
#'   analysis without small scars, for example `min_length = 10` for 10 mm.
#' @param keep Names of per-core columns (for example `"Typology"`) to carry
#'   into the output. Each must be constant within a core.
#' @param quiet Logical. Suppress the summary message? Default `FALSE`.
#' @inheritParams compute_spi
#'
#' @return A tibble with one row per core, in order of first appearance: the
#'   ID column, the `keep` columns, `n_scar` (scars used), and for each of
#'   SPI, ER and IR the columns selected by `raw`, `rarefied` and
#'   `random_interval`: `*_raw`, `*_rarefied`, `*_lower`, `*_upper` and
#'   `*_vs_random` (`"within"` or `"outside"`).
#'
#' @details
#' Rows with missing coordinates or zero length are left out, and so are
#' scars shorter than `min_length`. Cores with fewer than three usable scars
#' get `NA` values. One message summarises what was left out and which cores
#' needed a Monte Carlo approximation of the rarefied values (more than
#' `max_comb` combinations).
#'
#' Use the rarefied columns to compare cores or groups of cores (for
#' example with [compare_indices()]): their expectation does not change with
#' the scar count. Use the raw values together with the random interval to
#' judge a single core: `"outside"` means the raw value departs from what
#' random directions produce at that core's scar count (SPI and ER above,
#' IR below the interval: more concentrated than random; SPI below the
#' interval: typical of opposed flaking). See [compute_spi()] for details.
#'
#' For each core the rarefied SPI, ER and IR are averaged over the same
#' combinations of scars.
#'
#' @examples
#' set.seed(1)
#' scars <- do.call(rbind, lapply(1:4, function(i) {
#'   n <- 8 + 4 * i
#'   s <- matrix(rnorm(3 * n, sd = 10), n, 3)
#'   d <- matrix(rnorm(3 * n), n, 3)
#'   d[, 1] <- d[, 1] + i / 2
#'   data.frame(ID = paste0("core", i), Typology = if (i <= 2) "A" else "B",
#'              Start_X = s[, 1], Start_Y = s[, 2], Start_Z = s[, 3],
#'              End_X = s[, 1] + 10 * d[, 1], End_Y = s[, 2] + 10 * d[, 2],
#'              End_Z = s[, 3] + 10 * d[, 3])
#' }))
#' spi_ei_batch(scars, keep = "Typology")
#'
#' # The same cores without scars shorter than 8 units
#' spi_ei_batch(scars, min_length = 8, keep = "Typology", random_interval = FALSE)
#'
#' @seealso [compute_spi()], [compute_ei()], [random_interval()],
#'   [cpi_pra_batch()], [spharm_spatiodirection()], [compare_indices()]
#' @export
spi_ei_batch <- function(data, id_col = "ID", start_cols = c("Start_X", "Start_Y", "Start_Z"),
                         end_cols = c("End_X", "End_Y", "End_Z"), min_length = 0,
                         keep = NULL, raw = TRUE, rarefied = TRUE,
                         random_interval = TRUE, m = 3, max_comb = 2000,
                         seed = 2226, quiet = FALSE) {
  .check_scar_table(data, id_col, start_cols, end_cols, keep, min_length = min_length)
  .check_index_args(raw, rarefied, random_interval, m, max_comb, seed)
  .check_flag(quiet, "quiet")
  sc <- .scar_cores(data, id_col, start_cols, end_cols, min_length)
  cores <- sc$cores
  n <- vapply(cores, function(x) x$n, 0L)

  rows <- lapply(cores, function(x) {
    if (x$n < 3L) {
      return(list(raw = c(SPI = NA_real_, ER = NA_real_, IR = NA_real_),
                  rar = c(SPI = NA_real_, ER = NA_real_, IR = NA_real_),
                  mc = FALSE))
    }
    ei <- .ei_raw(x$v_hat)
    rr <- if (rarefied) .rarefied_fabric(x$v_hat, m, max_comb, seed) else
      list(values = c(SPI = NA_real_, ER = NA_real_, IR = NA_real_), enumerated = TRUE)
    list(raw = c(SPI = .spi_raw(x$v_hat), ER = ei[["ER"]], IR = ei[["IR"]]),
         rar = rr$values, mc = identical(rr$enumerated, FALSE))
  })

  out <- .core_frame(data, id_col, keep, sc$ids)
  out$n_scar <- unname(n)
  for (p in .RI_PARAMS) {
    value <- vapply(rows, function(r) r$raw[[p]], 0)
    if (raw) out[[paste0(p, "_raw")]] <- unname(value)
    if (rarefied) out[[paste0(p, "_rarefied")]] <- unname(vapply(rows, function(r) r$rar[[p]], 0))
    if (random_interval) {
      b <- .ri_bounds(p, ifelse(n >= 3L, n, NA_integer_))
      out[[paste0(p, "_lower")]] <- b$lower
      out[[paste0(p, "_upper")]] <- b$upper
      out[[paste0(p, "_vs_random")]] <- .ri_position(unname(value), b$lower, b$upper)
    }
  }

  few <- names(n)[n < 3L]
  mc <- names(rows)[vapply(rows, function(r) r$mc, TRUE)]
  big <- names(n)[n > .ri_max_n()]
  .say(c(.dropped_notes(cores),
         if (length(few)) sprintf("%d core(s) with fewer than 3 usable scars are NA: %s.",
                                  length(few), .fmt_ids(few)),
         if (length(mc)) sprintf(paste("%d core(s) have more than %d combinations of %d scars;",
                                       "their rarefied values average %d random combinations",
                                       "(seed %s): %s."),
                                 length(mc), as.integer(max_comb), as.integer(m),
                                 as.integer(max_comb), if (is.null(seed)) "NULL" else seed,
                                 .fmt_ids(mc)),
         if (random_interval && length(big))
           sprintf(paste("The random interval of ER and IR is tabulated up to n = %d;",
                         "it is NA for %d core(s) with more scars: %s."),
                   .ri_max_n(), length(big), .fmt_ids(big))),
       quiet)
  tibble::as_tibble(out)
}
