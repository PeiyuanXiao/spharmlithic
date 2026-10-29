# Synthetic cores for the scar-pattern tests --------------------------------

# A scar table of several cores; core i has n[i] scars whose directions lean
# towards +X by lean[i].
make_cores <- function(n = c(12, 8, 20), lean = c(1, 0, 0.5), seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_along(n), function(i) {
    s <- matrix(rnorm(3 * n[i], sd = 10), n[i], 3)
    d <- matrix(rnorm(3 * n[i]), n[i], 3)
    d[, 1] <- d[, 1] + lean[i]
    data.frame(ID = paste0("core", i), Typology = c("A", "B", "A")[(i - 1) %% 3 + 1],
               Start_X = s[, 1], Start_Y = s[, 2], Start_Z = s[, 3],
               End_X = s[, 1] + 5 * d[, 1], End_Y = s[, 2] + 5 * d[, 2],
               End_Z = s[, 3] + 5 * d[, 3], stringsAsFactors = FALSE)
  }))
}

# A disc core of radius r in the z = 0 plane: k scars start on the rim and
# run at `cpi` degrees from the outward radius (180 = to the centre), dipping
# by `pra` degrees below the plane.
make_disc <- function(k = 12, r = 30, cpi = 180, pra = 20, id = "disc") {
  t <- seq(0, 2 * pi, length.out = k + 1)[-(k + 1)]
  start <- cbind(r * cos(t), r * sin(t), 0)
  a <- cpi * pi / 180
  radial <- cbind(cos(t), sin(t))
  inplane <- cbind(radial[, 1] * cos(a) - radial[, 2] * sin(a),
                   radial[, 1] * sin(a) + radial[, 2] * cos(a))
  p <- pra * pi / 180
  dir <- cbind(inplane * cos(p), -sin(p))
  end <- start + 10 * dir
  data.frame(ID = id, Start_X = start[, 1], Start_Y = start[, 2], Start_Z = start[, 3],
             End_X = end[, 1], End_Y = end[, 2], End_Z = end[, 3])
}

disc_outline <- function(r = 30, n = 200) {
  u <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  cbind(X = r * cos(u), Y = r * sin(u), Z = 0)
}

# Rotate and move every coordinate column of a scar table.
rigid_move <- function(d, R, shift) {
  for (p in list(c("Start_X", "Start_Y", "Start_Z"), c("End_X", "End_Y", "End_Z"))) {
    M <- as.matrix(d[, p]) %*% t(R)
    d[, p] <- sweep(M, 2, shift, "+")
  }
  d
}

random_rotation <- function(seed = 3) {
  set.seed(seed)
  q <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  if (det(q) < 0) q[, 1] <- -q[, 1]
  q
}

# A minimal IGES file whose outline is the closed polyline through `pts`
# (rows), written as degree-1 B-spline segments (entity 126) that split the
# loop at `cut` and store the second segment reversed.
write_iges_polyline <- function(path, pts, cut = ceiling(nrow(pts) / 2), extra_type = NULL) {
  pts <- rbind(pts, pts[1, ])
  segs <- list(pts[1:cut, , drop = FALSE], pts[nrow(pts):cut, , drop = FALSE])
  fmt <- function(x) sprintf("%.6f", x)
  d_lines <- character(0)
  p_lines <- character(0)
  pseq <- 0L
  for (j in seq_along(segs)) {
    s <- segs[[j]]
    K <- nrow(s) - 1L
    knots <- c(0, 0:K, K)
    par <- c("126", K, 1, 0, 0, 1, 0, fmt(knots), fmt(rep(1, K + 1)),
             fmt(as.vector(t(s))), fmt(0), fmt(K), "0", "0", "1")
    txt <- paste0(paste(par, collapse = ","), ";")
    de <- 2L * j - 1L
    chunks <- regmatches(txt, gregexpr(".{1,64}", txt))[[1]]
    first <- pseq + 1L
    for (ch in chunks) {
      pseq <- pseq + 1L
      p_lines <- c(p_lines, sprintf("%-64s%8d%s%7d", ch, de, "P", pseq))
    }
    d_lines <- c(d_lines,
                 sprintf("%8d%8d%-56s%s%7d", 126L, first, "", "D", de),
                 sprintf("%8d%-64s%s%7d", 126L, "", "D", de + 1L))
  }
  if (!is.null(extra_type)) {
    de <- length(d_lines) + 1L
    d_lines <- c(d_lines, sprintf("%8d%8d%-56s%s%7d", extra_type, 1L, "", "D", de),
                 sprintf("%8d%-64s%s%7d", extra_type, "", "D", de + 1L))
  }
  s_line <- sprintf("%-72s%s%7d", "test outline", "S", 1L)
  g_line <- sprintf("%-72s%s%7d", "1H,,1H;;", "G", 1L)
  writeLines(c(s_line, g_line, d_lines, p_lines), path)
  invisible(path)
}
