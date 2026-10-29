# data-raw/random_interval.R -- builds R/sysdata.rda: the table of the random
# interval of SPI, ER and IR that random_interval(), compute_spi(),
# compute_ei() and spi_ei_batch() read.
#
#   Rscript data-raw/random_interval.R        # from the package root
#
# Source: results/baseline/baseline_quantiles_fine.csv of the article
# repository (Xiao, Li and Marwick, in prep.; CRF_ROOT, default
# H:/Core_reduction_framework), written there by analysis/scripts/baseline.R:
#   SPI      the exact distribution of the resultant length of n unit vectors
#            uniform on the sphere (Rayleigh's random flight), n = 3..200
#   ER, IR   5000 simulated cores of n directions uniform on the sphere per n
#            (seed SEED + n), n = 3..200
# The interval is the central 95%: the quantiles at 0.025 and 0.975, as the
# article's baseline_table.csv (null_lower05, null_upper95).
#
# base R's parser is correctly rounded, so the numbers here agree to the bit
# with what the article's scar_baseline() reads.

CRF_ROOT <- Sys.getenv("CRF_ROOT", "H:/Core_reduction_framework")
src <- file.path(CRF_ROOT, "results", "baseline", "baseline_quantiles_fine.csv")
pub <- file.path(CRF_ROOT, "results", "baseline", "baseline_table.csv")
stopifnot(file.exists(src), file.exists(pub))

fine <- utils::read.csv(src, stringsAsFactors = FALSE)
fine <- fine[fine$null == "isotropic", ]
pick <- function(p) fine[abs(fine$prob - p) < 1e-9, c("parameter", "n", "q")]
lo <- pick(0.025)
hi <- pick(0.975)
key <- function(d) paste(d$parameter, d$n)
stopifnot(setequal(key(lo), key(hi)))

.random_interval_table <- data.frame(
  parameter = lo$parameter,
  n = as.integer(lo$n),
  lower = lo$q,
  upper = hi$q[match(key(lo), key(hi))],
  stringsAsFactors = FALSE
)
.random_interval_table <- .random_interval_table[
  order(match(.random_interval_table$parameter, c("SPI", "ER", "IR")),
        .random_interval_table$n), ]
rownames(.random_interval_table) <- NULL

# every parameter at every n = 3..200, nothing missing
stopifnot(!anyNA(.random_interval_table),
          all(table(.random_interval_table$parameter) == 198L),
          all(.random_interval_table$lower < .random_interval_table$upper))

# agrees with the article's published table (n = 3..60) to the bit
tab <- utils::read.csv(pub, stringsAsFactors = FALSE, comment.char = "#")
tab <- tab[tab$null == "isotropic", ]
m <- match(paste(tab$parameter, tab$n), key(.random_interval_table))
stopifnot(!anyNA(m),
          identical(.random_interval_table$lower[m], tab$null_lower05),
          identical(.random_interval_table$upper[m], tab$null_upper95))

save(.random_interval_table, file = file.path("R", "sysdata.rda"),
     compress = "xz", version = 2)
cat(sprintf("wrote R/sysdata.rda: %d rows (SPI, ER, IR x n = %d..%d)\n",
            nrow(.random_interval_table), min(.random_interval_table$n),
            max(.random_interval_table$n)))
