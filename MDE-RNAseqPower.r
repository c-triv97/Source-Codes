############################################################
# MDE calculator for pseudobulk DE in snRNA-seq
# Based on Hart et al. (2013) / RNASeqPower (Negative Binomial)
# Two-group design, two-sided alpha.
#
# References:
# - RNASeqPower & Hart et al.: https://www.rdocumentation.org/packages/RNASeqPower/versions/1.12.0/topics/rnapower
#   and https://bioconductor.org/packages/.../RNASeqPower.pdf; power.hart notes: https://rdrr.io/cran/FDRsamplesize2/man/power.hart.html
# - Pseudobulk best practices: https://www.sc-best-practices.org/conditions/differential_gene_expression.html
#   and https://hbctraining.github.io/Pseudobulk-for-scRNAseq/lessons/03_pseudobulk_DESeq2.html
############################################################

# --------- 1) Utilities: normal quantiles ---------
# We'll use qnorm for standard normal quantiles.
z_alpha_two_sided <- function(alpha) {
  # Z_{1 - alpha/2}
  qnorm(1 - alpha/2)
}

z_power <- function(power) {
  # Z_{power}
  qnorm(power)
}

# --------- 2) Hart et al. MDE formula ----------
# Hart et al. sample size equation (per group) for detecting a given effect:
#   n = [ 2 * (Z_{1-a/2} + Z_power)^2 * (1/mu + CV^2) ] / [ (ln(effect))^2 ]
# Solve for effect (minimum detectable fold-change) for given n, alpha, power, mu, CV:
mde_fold_change <- function(n, alpha = 0.05, power = 0.8, mu, CV) {
  za <- z_alpha_two_sided(alpha)
  zp <- z_power(power)
  term <- 2 * (za + zp)^2 * (1/mu + CV^2) / n
  exp(sqrt(term)) # returns fold-change (FC)
}

# Convert FC to log2FC for convenience
log2fc <- function(fc) log(fc, base = 2)

# --------- 3) Example grid over mu and CV ----------
# This reproduces the numbers shown previously for n=3 per group.
alpha <- 0.05
powers <- c(0.8, 0.9)   # 80% and 90% power
n_per_group <- 4        # change this here for n per group (3, 4, 5, ...)
mu_grid <- c(5, 10, 20, 50, 100)
cv_grid <- c(0.4, 0.6, 0.8, 1.0)

grid_results <- list()

for (pwr in powers) {
  cat("\n==== Results for power =", pwr, "====\n")
  for (mu in mu_grid) {
    for (cv in cv_grid) {
      fc <- mde_fold_change(n = n_per_group, alpha = alpha, power = pwr, mu = mu, CV = cv)
      cat(sprintf("mu=%-4.0f, CV=%.1f -> MDE %5.2f× (log2FC=%.2f)\n", mu, cv, fc, log2fc(fc)))
      grid_results[[length(grid_results) + 1]] <- data.frame(
        power = pwr, mu = mu, CV = cv, MDE_FC = fc, MDE_log2FC = log2fc(fc)
      )
    }
  }
}

# Combine into a single data.frame for further use (plotting, saving)
grid_df <- do.call(rbind, grid_results)

# --------- 4) Connect μ to 100k nuclei via pseudobulk ----------
# Total nuclei: 100,000 across 6 samples => ~16,667 nuclei per sample.
# For a given cell type, the per-sample gene counts μ ≈ cells_per_sample * cell_type_prop * mean_UMI_per_cell_for_gene.
# NOTE: This is a simplification (UMI for the specific gene); for DE analyses, μ is per gene.
cells_total <- 100000
n_samples <- 8 # total samples (e.g., 4 vs 4)
cells_per_sample <- cells_total / n_samples  

# Define scenarios of cell-type abundance and mean UMI per cell for the gene of interest
scenarios <- expand.grid(
  cell_type_prop      = c(0.10, 0.05, 0.01),   # 10%, 5%, 1% cell type
  umi_per_cell_for_gene = c(0.05, 0.10)        # mean UMIs per cell for the gene
)

# Choose CV values to demonstrate
cv_values <- c(0.4, 0.6, 0.8)
power_target <- 0.8

cat(sprintf("\n==== Pseudobulk-derived μ scenarios (power 80%%, alpha 0.05, n=%d) ====\n", n_per_group))
pseudo_results <- list()
for (i in seq_len(nrow(scenarios))) {
  prop <- scenarios$cell_type_prop[i]
  umi_per_cell <- scenarios$umi_per_cell_for_gene[i]
  mu <- cells_per_sample * prop * umi_per_cell
  for (cv in cv_values) {
    fc <- mde_fold_change(n = n_per_group, alpha = alpha, power = power_target, mu = mu, CV = cv)
    cat(sprintf("prop=%4.2f, umi/cell=%4.2f, mu=%6.1f, CV=%.1f -> MDE %5.2f× (log2FC=%.2f)\n",
                prop, umi_per_cell, mu, cv, fc, log2fc(fc)))
    pseudo_results[[length(pseudo_results) + 1]] <- data.frame(
      prop = prop, umi_per_cell = umi_per_cell, mu = mu, CV = cv,
      MDE_FC = fc, MDE_log2FC = log2fc(fc),
      alpha = alpha, power = power_target, n_per_group = n_per_group
    )
  }
}
pseudo_df <- do.call(rbind, pseudo_results)

# --------- 5) Optional cross-check using RNASeqPower::rnapower ----------
# rnapower() can solve for sample size or power given effect, depth (mu), CV, etc.
# Here we invert: choose an effect and ask rnapower for resulting power.
# NOTE: rnapower's "depth" argument is the average coverage per gene (≈ μ here).
# If you don't have RNASeqPower installed, uncomment install line.

# install.packages("BiocManager")
# BiocManager::install("RNASeqPower")

suppressWarnings({
  has_rnapower <- requireNamespace("RNASeqPower", quietly = TRUE)
})

if (has_rnapower) {
  library(RNASeqPower)
  cat("\n==== Cross-check with RNASeqPower::rnapower (selected scenario) ====\n")
  # Pick scenario: mu=83.3, CV=0.6, n=3 per group, alpha=0.05
  mu_chk <- 83.3
  cv_chk <- 0.6
  n_chk  <- 3
  alpha_chk <- 0.05
  # Use the effect equal to our computed MDE at 80% power to see if rnapower returns ~0.8
  effect_mde <- mde_fold_change(n_chk, alpha_chk, 0.8, mu_chk, cv_chk)
  power_est  <- rnapower(depth = mu_chk, n = n_chk, cv = cv_chk, effect = effect_mde, alpha = alpha_chk)
  cat(sprintf("rnapower check -> mu=%.1f, CV=%.1f, n=%d, effect=%.2f×: power ≈ %.3f\n",
              mu_chk, cv_chk, n_chk, effect_mde, power_est))
} else {
  cat("\nNOTE: RNASeqPower not installed; skipping cross-check. See docs:\n",
      "https://www.rdocumentation.org/packages/RNASeqPower/versions/1.12.0/topics/rnapower\n")
}

# --------- 6) (Optional) Quick plot of MDE vs mu for given CV ----------
# This illustrates diminishing returns of increasing μ relative to adding replicates.
suppressWarnings({
  has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
})
if (has_ggplot2) {
  library(ggplot2)
  p <- ggplot(subset(grid_df, power == 0.8 & CV %in% c(0.4,0.6,0.8)),
              aes(x = mu, y = MDE_log2FC, color = factor(CV), group = CV)) +
    geom_line(size = 1.1) + geom_point(size = 2) +
    scale_color_brewer(palette = "Dark2", name = "CV") +
    labs(title = sprintf("Minimum detectable log2 fold-change vs. μ (n = %d, α = 0.05, power = 80%%)", n_per_group),
         x = expression(mu~"(mean counts per gene per sample)"),
         y = "MDE (log2 fold-change)") +
    theme_bw(base_size = 12)
  print(p)
}

############################################################
# End of script
############################################################