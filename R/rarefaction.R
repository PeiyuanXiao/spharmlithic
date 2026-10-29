# ==============================================================================
# rarefaction.R
# Raw and rarefied SPI, ER and IR of one set of unit directions (internal).
#
# Rarefied indices: every combination of m scars on a core is taken (m = 3 by
# default), SPI, ER and IR are computed on each, and each is averaged over the
# combinations. For scars drawn independently from one pattern the mean over
# all m-subsets is a U-statistic of order m (Hoeffding 1948), whose expectation
# does not depend on the scar count n. Above `max_comb` combinations, max_comb
# distinct combinations are drawn at random. n == m gives the raw values;
# n < m gives NA.
# ==============================================================================

# SPI of the unit rows of V; `w` optional weights (scar lengths). The
# arithmetic of compute_spi() before 0.2.0, kept to the bit.
.spi_raw <- function(V, w = NULL) {
  dx <- V[, 1]; dy <- V[, 2]; dz <- V[, 3]
  if (!is.null(w)) {
    dx <- dx * w
    dy <- dy * w
    dz <- dz * w
  }
  sqrt(sum(dx)^2 + sum(dy)^2 + sum(dz)^2) / sum(sqrt(dx^2 + dy^2 + dz^2))
}

# ER, IR and the eigenvalues of the orientation tensor of the unit rows of V.
.ei_raw <- function(V) {
  T_mat  <- (t(V) %*% V) / nrow(V)
  lambda <- pmax(sort(eigen(T_mat, symmetric = TRUE)$values, decreasing = TRUE), 0)
  ok <- lambda[1] > 1e-10
  c(ER = if (ok) 1 - lambda[2] / lambda[1] else NA_real_,
    IR = if (ok) lambda[3] / lambda[1] else NA_real_,
    lambda1 = lambda[1], lambda2 = lambda[2], lambda3 = lambda[3])
}

# Rarefied SPI, ER and IR of the unit rows of V. Returns list(values =
# c(SPI, ER, IR), n_comb, enumerated); values are NA when nrow(V) < m.
.rarefied_fabric <- function(V, m = 3L, max_comb = 2000L, seed = 2226L) {
  n <- nrow(V)
  if (n < m) {
    return(list(values = c(SPI = NA_real_, ER = NA_real_, IR = NA_real_),
                n_comb = 0L, enumerated = NA))
  }
  cmb <- .with_seed(seed, .rarefied_combinations(n, m, max_comb))
  val <- .rarefied_on_combinations(V, cmb$idx)
  out <- colMeans(val, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  list(values = out, n_comb = nrow(cmb$idx), enumerated = cmb$enumerated)
}

.check_rarefaction_args <- function(m, max_comb) {
  if (length(m) != 1L || !is.numeric(m) || is.na(m) || m < 3 || m != round(m)) {
    stop("`m` must be a single whole number >= 3.", call. = FALSE)
  }
  if (length(max_comb) != 1L || !is.numeric(max_comb) || is.na(max_comb) ||
      max_comb < 1 || max_comb != round(max_comb)) {
    stop("`max_comb` must be a single whole number >= 1.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---- which combinations ----------------------------------------------------------
# Returns idx, a K x m matrix of scar indices (one combination per row, sorted
# within the row), and whether it is the full enumeration.
#
# Three regimes, all uniform over combinations and all without repetition:
#   C(n, m) <= max_comb        every combination
#   C(n, m) <= 20 * max_comb   list them all, then keep max_comb of the rows
#   larger                     draw random m-subsets until max_comb distinct
# The middle regime exists because rejection sampling slows down when nearly
# every combination is wanted; the last because listing C(40, 5) = 658,008
# rows to keep 2,000 is wasteful.
.rarefied_combn_cache <- new.env(parent = emptyenv())

.rarefied_all_combinations <- function(n, m) {
  key <- paste(n, m)
  if (is.null(.rarefied_combn_cache[[key]]))
    .rarefied_combn_cache[[key]] <- t(utils::combn(n, m))
  .rarefied_combn_cache[[key]]
}

.rarefied_combinations <- function(n, m, max_comb) {
  total <- choose(n, m)
  if (total <= max_comb)
    return(list(idx = .rarefied_all_combinations(n, m), enumerated = TRUE))

  if (total <= 20 * max_comb) {
    all <- .rarefied_all_combinations(n, m)
    rows <- sort(sample.int(nrow(all), max_comb))
    return(list(idx = all[rows, , drop = FALSE], enumerated = FALSE))
  }

  idx <- matrix(integer(0), 0L, m)
  keys <- character(0)
  while (nrow(idx) < max_comb) {
    k <- max_comb - nrow(idx) + 16L
    draw <- .rarefied_random_subsets(k, n, m)
    dkey <- do.call(paste, c(as.data.frame(draw), sep = "."))
    fresh <- !duplicated(dkey) & !(dkey %in% keys)
    idx <- rbind(idx, draw[fresh, , drop = FALSE])
    keys <- c(keys, dkey[fresh])
  }
  list(idx = idx[seq_len(max_comb), , drop = FALSE], enumerated = FALSE)
}

# k independent uniform m-subsets of 1..n, one per row, sorted within the row:
# each row is the first m entries of a random permutation of 1..n.
.rarefied_random_subsets <- function(k, n, m) {
  row <- rep(seq_len(k), times = n)
  col <- rep(seq_len(n), each = k)
  perm <- matrix(col[order(row, stats::runif(k * n))], k, n, byrow = TRUE)
  sel <- as.vector(perm[, seq_len(m), drop = FALSE])
  matrix(sel[order(rep(seq_len(k), m), sel)], k, m, byrow = TRUE)
}

# ---- the three indices on every combination at once ----------------------------
# V     n x 3 unit directions
# idx   K x m scar indices
# Returns a K x 3 matrix, one row per combination.
.rarefied_on_combinations <- function(V, idx) {
  K <- nrow(idx)
  m <- ncol(idx)
  pick <- function(M, j) matrix(M[idx, j], K, m)

  # directions: SPI from the resultant, ER and IR from the orientation tensor
  x <- pick(V, 1L); y <- pick(V, 2L); z <- pick(V, 3L)
  spi <- sqrt(rowMeans(x)^2 + rowMeans(y)^2 + rowMeans(z)^2)
  lam <- .sym3_eigenvalues(rowMeans(x * x), rowMeans(y * y), rowMeans(z * z),
                           rowMeans(x * y), rowMeans(x * z), rowMeans(y * z))
  lam <- pmax(lam, 0)                      # as .ei_raw()
  ok1 <- lam[, 1] > 1e-10
  er <- ifelse(ok1, 1 - lam[, 2] / lam[, 1], NA_real_)
  ir <- ifelse(ok1, lam[, 3] / lam[, 1], NA_real_)

  cbind(SPI = spi, ER = er, IR = ir)
}

# Eigenvalues, descending, of many symmetric 3 x 3 matrices at once, one per
# element of the six input vectors (Smith 1961, the trigonometric solution of
# the characteristic cubic). The arccos loses digits when two eigenvalues
# nearly coincide (|det(B)/2| near 1), so those few matrices go to eigen()
# instead; everywhere else the two agree to rounding error. Returns a K x 3
# matrix.
.sym3_eigenvalues <- function(a11, a22, a33, a12, a13, a23) {
  q  <- (a11 + a22 + a33) / 3
  p2 <- (a11 - q)^2 + (a22 - q)^2 + (a33 - q)^2 +
        2 * (a12^2 + a13^2 + a23^2)
  p  <- sqrt(p2 / 6)
  flat <- !(p > 0)                          # all three eigenvalues equal to q
  p[flat] <- 1
  b11 <- (a11 - q) / p; b22 <- (a22 - q) / p; b33 <- (a33 - q) / p
  b12 <- a12 / p;       b13 <- a13 / p;       b23 <- a23 / p
  detb <- b11 * (b22 * b33 - b23 * b23) -
          b12 * (b12 * b33 - b23 * b13) +
          b13 * (b12 * b23 - b22 * b13)
  r <- pmin(pmax(detb / 2, -1), 1)
  phi <- acos(r) / 3
  l1 <- q + 2 * p * cos(phi)
  l3 <- q + 2 * p * cos(phi + 2 * pi / 3)
  l2 <- 3 * q - l1 - l3
  l1[flat] <- q[flat]; l2[flat] <- q[flat]; l3[flat] <- q[flat]
  out <- cbind(l1, l2, l3)

  near <- which(!flat & abs(r) > 1 - 1e-6)
  for (i in near) {
    A <- matrix(c(a11[i], a12[i], a13[i],
                  a12[i], a22[i], a23[i],
                  a13[i], a23[i], a33[i]), 3L, 3L)
    out[i, ] <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  }
  out
}
