## ----setup, include = FALSE---------------------------------------------------
# This vignette uses several Suggested packages. If any are missing we still
# render the prose and show the code, but skip evaluation, so the vignette
# never fails to build.
can_run <- requireNamespace("vegan",        quietly = TRUE) &&
           requireNamespace("compositions", quietly = TRUE) &&
           requireNamespace("ggplot2",      quietly = TRUE) &&
           requireNamespace("tidyr",        quietly = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  warning  = FALSE,
  message  = FALSE,
  fig.width  = 6,
  fig.height = 4,
  eval     = can_run
)

## ----load-data----------------------------------------------------------------
# load the tables
scar  <- read.csv(
  system.file("extdata", "exp_cores_scar.csv",  package = "spharmlithic"),
  check.names = FALSE)
morph <- read.csv(
  system.file("extdata", "exp_cores_morph.csv", package = "spharmlithic"),
  check.names = FALSE)

# set metadata
strategy_levels <- c("Unidirectional", "Bidirectional",
                     "Levallois", "Discoid", "Multiplatform")
scar$Typology  <- factor(scar$Typology,  levels = strategy_levels)
morph$Typology <- factor(morph$Typology, levels = strategy_levels)

# prepare a palette for visualization
typology_colors <- c(
  "Unidirectional" = "#BA8530",
  "Bidirectional"  = "#788C4A",
  "Levallois"      = "#4A6E8A",
  "Discoid"        = "#802520",
  "Multiplatform"  = "#8A7A68"
)

# display summary table
knitr::kable(
  as.data.frame(table(Strategy = scar$Typology)),
  col.names = c("Reduction strategy", "n"),
  caption   = "Experimental cores per reduction strategy")

## ----preview------------------------------------------------------------------
knitr::kable(
  head(scar[, 1:6]),
  caption = "First few rows of the scar-direction power spectra")

# Confirm the spectra are closed (compositional)
range(rowSums(scar[, grep("^power_l", names(scar))]))

## ----feature-matrices---------------------------------------------------------
scar_cols   <- paste0("power_l", 1:6) # first 6 degrees of the power spectrum
morph_cols  <- paste0("power_l", 1:8) # first 8 degrees of the power spectrum
scar_feat   <- scar[,  scar_cols]
morph_feat  <- morph[, morph_cols]
groups      <- scar$Typology     

## ----ilr----------------------------------------------------------------------
library(spharmlithic)
library(dplyr)

# ILR coordinates: D parts -> (D - 1) real-valued coordinates per specimen.
scar_ilr  <- make_ilr(scar_feat)    # 6 degrees -> 5 ILR coordinates
morph_ilr <- make_ilr(morph_feat)   # 8 degrees -> 7 ILR coordinates

# first few rows
head(scar_ilr)

## ----degree-selection, fig.width = 7, fig.height = 3.4, fig.cap = "Per-degree diagnostics for the two power spectra (experimental cores, n = 58): across-specimen coefficient of variation (left) and cumulative power (right), each overlaying the two domains. Red dashed lines mark the retained truncations (degree 6 for scar direction, degree 8 for morphology); grey dashed lines mark CV = 100% and the 95% / 99% cumulative-power levels."----
library(ggplot2)

deg_df <- rbind(
  degree_diagnostics(scar,  "SP-SPHARM (scar direction)"),
  degree_diagnostics(morph, "M-SPHARM (morphology)"))
deg_df$descriptor <- factor(deg_df$descriptor,
  levels = c("M-SPHARM (morphology)", "SP-SPHARM (scar direction)"))

descriptor_colors <- c(
  "M-SPHARM (morphology)"      = "#4A6E8A",
  "SP-SPHARM (scar direction)" = "#BA8530")

# Long form so the two diagnostics share one faceted figure.
deg_long <- tidyr::pivot_longer(deg_df, c(cv_pct, cumul_pct),
                                names_to = "metric", values_to = "value")
deg_long$metric <- factor(deg_long$metric, c("cv_pct", "cumul_pct"),
  labels = c("Across-specimen CV (%)", "Cumulative power (%)"))

# Panel-specific reference lines: CV = 100% (noise threshold), and the 95 / 99%
# cumulative-power levels.
ref_lines <- data.frame(
  metric = factor(c("Across-specimen CV (%)",
                    "Cumulative power (%)", "Cumulative power (%)"),
                  levels = levels(deg_long$metric)),
  yint   = c(100, 95, 99))

ggplot(deg_long, aes(degree, value, colour = descriptor)) +
  geom_hline(
    data = ref_lines,
    aes(yintercept = yint),
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = c(6, 8),
    colour = "red",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = descriptor_colors) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  facet_wrap( ~ metric, scales = "free_y") +
  labs(x = "Spherical harmonic degree", y = NULL, colour = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "bottom")

## ----pca----------------------------------------------------------------------
pca_scar  <- prcomp(scar_ilr)
pca_morph <- prcomp(morph_ilr)

# Variance explained by the first few components, both domains
var_tab <- data.frame(
  PC = paste0("PC", 1:5),
  `Scar direction (%)` = round(100 * summary(pca_scar)$importance[2, 1:5], 1),
  `Morphology (%)`     = round(100 * summary(pca_morph)$importance[2, 1:5], 1),
  check.names = FALSE)


knitr::kable(var_tab, row.names = FALSE,
  caption = "Variance explained by the leading principal components, by domain")

## ----pca-plot, fig.cap = "Cores in the space of the first two principal components of the ILR coordinates, by domain, coloured by reduction strategy, with convex hulls outlining each group. Axis scales are free; variance explained is in the table above."----
library(ggplot2)

pca_df <- rbind(
  data.frame(Domain = "Scar direction",
             PC1 = pca_scar$x[, 1],  PC2 = pca_scar$x[, 2],  Strategy = groups),
  data.frame(Domain = "Morphology",
             PC1 = pca_morph$x[, 1], PC2 = pca_morph$x[, 2], Strategy = groups))

pca_df$Domain   <- factor(pca_df$Domain, levels = c("Scar direction", "Morphology"))
pca_df$Strategy <- factor(pca_df$Strategy, levels = strategy_levels)

# Convex hulls per strategy, matching the project's LDA / ternary scatterplots
hull_df <- pca_df %>%
  group_by(Domain, Strategy) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

ggplot(pca_df, aes(PC1, PC2, colour = Strategy)) +
  geom_hline(
    yintercept = 0,
    colour = "grey50",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    colour = "grey50",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_polygon(
    data = hull_df,
    aes(fill = Strategy, group = Strategy),
    alpha = 0.25,
    colour = NA
  ) +
  geom_point(size = 2,
             alpha = 0.85,
             shape = 16) +
  scale_colour_manual(values = typology_colors) +
  scale_fill_manual(values  = typology_colors) +
  facet_wrap( ~ Domain, scales = "free") +
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom")

## ----spectrum-plot, fig.cap = "Mean power spectrum by reduction strategy and domain (ribbon = +/- 1 SE). Values are closed proportions."----
library(tidyr)

# The two domains now keep different numbers of degrees, so pivot each to long
# form (a shared schema) before combining.
spec_long <- rbind(
  data.frame(
    Domain = "Scar direction",
    Strategy = groups,
    scar_feat,
    check.names = FALSE
  ) |>
    tidyr::pivot_longer(
      cols = all_of(scar_cols),
      names_to = "degree",
      values_to = "power"
    ),
  data.frame(
    Domain = "Morphology",
    Strategy = groups,
    morph_feat,
    check.names = FALSE
  ) |>
    tidyr::pivot_longer(
      cols = all_of(morph_cols),
      names_to = "degree",
      values_to = "power"
    )
)

spec_summ <- spec_long %>%
  mutate(degree =  readr::parse_number(degree),
        Domain =   factor(Domain, levels = c("Scar direction", "Morphology")),
        Strategy = factor(Strategy, levels = strategy_levels)) %>% 
  group_by(Domain, Strategy, degree) %>%
  summarise(
    power.m =  mean(power, na.rm = TRUE),
    power.se = sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power))),
    .groups = "drop"
  )

ggplot(spec_summ,
       aes(degree, power.m, colour = Strategy, fill = Strategy)) +
  geom_ribbon(
    aes(ymin = power.m - power.se, ymax = power.m + power.se),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = typology_colors) +
  scale_fill_manual(values  = typology_colors) +
  facet_wrap( ~ Domain, scales = "free") +
  labs(x = "Spherical harmonic degree",
       y = "Mean power (proportion)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

## ----permanova----------------------------------------------------------------
library(vegan)

grp_df <- data.frame(Strategy = groups)

set.seed(1)
pm_scar  <- adonis2(dist(scar_ilr)  ~ Strategy,
                    data = grp_df, permutations = 999)
pm_morph <- adonis2(dist(morph_ilr) ~ Strategy,
                    data = grp_df, permutations = 999)

## ----permanova-table, echo = FALSE--------------------------------------------
permanova_summary <- data.frame(
  Domain    = c("Scar direction", "Morphology"),
  R2        = c(pm_scar$R2[1],        pm_morph$R2[1]),
  `Pseudo-F`= c(pm_scar$F[1],         pm_morph$F[1]),
  `p-value` = c(pm_scar$`Pr(>F)`[1],  pm_morph$`Pr(>F)`[1]),
  check.names = FALSE)

knitr::kable(permanova_summary, digits = 3,
  caption = "PERMANOVA for differences among reduction strategies, by domain")

## ----pairwise-permanova-------------------------------------------------------
set.seed(1)

pw_scar  <- pairwise_permanova(scar_ilr,  groups)
pw_morph <- pairwise_permanova(morph_ilr, groups)

## ----pairwise-table, echo = FALSE---------------------------------------------
pw_all <- rbind(
  cbind(Domain = "Scar direction", pw_scar),
  cbind(Domain = "Morphology",     pw_morph))

knitr::kable(pw_all, digits = 3, row.names = FALSE,
  caption = "Pairwise PERMANOVA between reduction strategies. p_adj: Holm-Bonferroni correction within each domain.")

## ----pairwise-permanova-plot, fig.cap = "Adjusted p-values for pairwise PERMANOVA tests. We plot  −log₁₀(p_adj) instead of raw adjusted p-values to spread out tiny p-values and make them more interpretable and visually informative. Bigger bars mean stronger evidence for signficant difference. Vertical dashed line indicates p = 0.05"----

pw_all_log <- 
 pw_all %>%
  mutate(
    sig = p_adj < 0.05,
    score = -log10(p_adj)
  )

ggplot(
  pw_all_log,
  aes(
    forcats::fct_reorder(pair, score),
    score,
    fill = sig
  )
) +
  geom_col() +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(-log[10](adjusted~p))
  ) +
  scale_fill_manual(
  name = "Test result",
  values = c("FALSE" = "grey70", 
             "TRUE" = "steelblue"),
  labels = c("not significant", "significant")
) +
  theme_minimal() +
  facet_wrap( ~Domain)


## ----betadisper---------------------------------------------------------------
# type = "centroid" gives the original PERMDISP of Anderson (2006): distances are measured to the group centroid (betadisper's default is the spatial median).
bd_scar  <- betadisper(dist(scar_ilr),  groups, type = "centroid")
bd_morph <- betadisper(dist(morph_ilr), groups, type = "centroid")

set.seed(1)

pt_scar  <- permutest(bd_scar,  permutations = 999)
pt_morph <- permutest(bd_morph, permutations = 999)

## ----betadisper-table, echo = FALSE-------------------------------------------
disp_summary <- data.frame(
  Domain    = c("Scar direction", "Morphology"),
  `F`       = c(pt_scar$tab$F[1],         pt_morph$tab$F[1]),
  `p-value` = c(pt_scar$tab$`Pr(>F)`[1],  pt_morph$tab$`Pr(>F)`[1]),
  check.names = FALSE)

knitr::kable(disp_summary, digits = 3,
  caption = "Homogeneity of multivariate dispersion (PERMDISP) by domain")

## ----betadisper-plot, fig.cap = "Distance to group centroid by reduction strategy and domain. Lower and tighter = more standardised."----
disp_df <- rbind(
  data.frame(Domain = "Scar direction", Strategy = groups, distance = bd_scar$distances),
  data.frame(Domain = "Morphology",     Strategy = groups, distance = bd_morph$distances))

disp_df$Domain   <- factor(disp_df$Domain, levels = c("Scar direction", "Morphology"))
disp_df$Strategy <- factor(disp_df$Strategy, levels = strategy_levels)

ggplot(disp_df,
       aes(Strategy, distance, colour = Strategy, fill = Strategy)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.25,
               linewidth = 0.35) +
  geom_jitter(
    width = 0.15,
    size = 1.6,
    alpha = 0.7,
    shape = 16
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 16,
    size = 2.4,
    colour = "white"
  ) +
  scale_colour_manual(values = typology_colors) +
  scale_fill_manual(values  = typology_colors) +
  facet_wrap( ~ Domain, scales = "free_y") +
  labs(x = NULL, y = "Distance to group centroid") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "none"
  )

