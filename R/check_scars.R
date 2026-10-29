# ==============================================================================
# check_scars.R
# One look at a scar table before analysing it.
# ==============================================================================

#' Check a scar table before analysis
#'
#' Reports in one message what the scar-pattern functions would make of a
#' scar table: how many cores and scars it holds, which scars would be left
#' out (missing coordinates, zero length, shorter than `min_length`), which
#' cores have fewer than three usable scars, and whether the `keep` columns
#' are constant within each core. Run it once after reading the data.
#'
#' @inheritParams spi_ei_batch
#'
#' @return Invisibly, a tibble with one row per core: the ID column, the
#'   `keep` columns, `n_rows` (rows in the table), `n_missing`, `n_zero`,
#'   `n_short` (rows left out for each reason), `n_scar` (usable scars) and
#'   `enough` (`n_scar` of at least 3).
#'
#' @examples
#' d <- data.frame(ID = rep(c("a", "b"), c(4, 2)),
#'                 Start_X = c(0, 1, 2, 3, 0, 1), Start_Y = 0, Start_Z = 0,
#'                 End_X = c(1, 1, 3, 5, 1, 2), End_Y = c(1, 0, 1, 1, 1, 1),
#'                 End_Z = 0)
#' check_scars(d)
#'
#' @seealso [spi_ei_batch()], [cpi_pra_batch()], [spharm_spatiodirection()]
#' @export
check_scars <- function(data, id_col = "ID", start_cols = c("Start_X", "Start_Y", "Start_Z"),
                        end_cols = c("End_X", "End_Y", "End_Z"), min_length = 0,
                        keep = NULL) {
  .check_scar_table(data, id_col, start_cols, end_cols, keep, min_length = min_length)
  sc <- .scar_cores(data, id_col, start_cols, end_cols, min_length)
  cores <- sc$cores
  cnt <- function(f) unname(vapply(cores, function(x) x[[f]], 0L))
  out <- .core_frame(data, id_col, keep, sc$ids)
  out$n_rows <- unname(lengths(split(seq_len(nrow(data)),
                                     factor(as.character(data[[id_col]]),
                                            levels = names(cores)))))
  out$n_missing <- cnt("n_missing")
  out$n_zero <- cnt("n_zero")
  out$n_short <- cnt("n_short")
  out$n_scar <- cnt("n")
  out$enough <- out$n_scar >= 3L

  n <- out$n_scar
  few <- names(cores)[!out$enough]
  .say(c(sprintf("%d core(s), %d scar row(s); usable scars per core: %d-%d (median %g).",
                 length(cores), nrow(data), min(n), max(n), stats::median(n)),
         .dropped_notes(cores),
         if (length(few)) sprintf("%d core(s) have fewer than 3 usable scars and will be NA: %s.",
                                  length(few), .fmt_ids(few))
         else "Every core has at least 3 usable scars.",
         if (length(keep)) sprintf("Constant within every core: %s.",
                                   paste(setdiff(keep, id_col), collapse = ", "))))
  invisible(tibble::as_tibble(out))
}
