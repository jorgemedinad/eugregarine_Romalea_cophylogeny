#############################################################
# Title: PACo and ParaFit Analyses of Host-Parasite Cophylogenetic Signal
# Creator: Jorge Medina-Duran
# Date: Apr/27/2025
#
# Description:
# This script performs global-fit cophylogenetic analyses between Romalea grasshoppers and gregarine parasites.
# It includes:
# - Loading phylogenies and association matrix
# - Performing PACo analysis with jackknife resampling
# - Plotting host-parasite link contributions and tanglegrams
# - Performing subgroup analyses for Boliviana clade
# - Visualizing global m² significance test
# - Running ParaFit global and link-specific significance tests
#
# Output:
# - PACo link contributions and jackknife confidence intervals
# - Tanglegram weighted by PACo contributions
# - Boxplots comparing Boliviana vs. other interactions
# - Density plots of residuals
# - Null distribution of m² statistic
# - ParaFit link significance results
#
# Notes:
# - Set working directory and adjust `output` path for reproducibility.
# - Figures and results are saved into `/figures`, `/paco`, and `/parafit` subdirectories.
#############################################################

# Load libraries
library(ape)
library(paco)
library(phytools)
library(ggplot2)
library(tidyr)

#### Define paths ----
output <- "YOUR_PATH/output/"

# Load trees and association matrix
parasite_tree <- read.tree(paste0(output, "trimmed_tree/gregarine_tree_trim.tre"))
host_tree <- read.tree(paste0(output, "trimmed_tree/romalea_tree_trim.tre"))
int <- read.csv(paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny_renamed.csv"), row.names = 1)

# Compute patristic distance matrices
ptree <- cophenetic(parasite_tree)
htree <- cophenetic(host_tree)

#### PACo analysis ----
D <- prepare_paco_data(htree, ptree, int)
D <- add_pcoord(D, correction = "cailliez")
D <- PACo(D, nperm = 10000, seed = 13, method = "r0", symmetric = TRUE, shuffled = TRUE)
D <- paco_links(D)

# Residuals and jackknife contributions
res <- residuals_paco(D$proc)
res_sq <- res^2
res_df <- data.frame(link = rownames(res), res = res, res_sq = res_sq, jack = D$jackknife)

#### Jackknife Confidence Intervals ----
HP.ones <- which(as.matrix(int) > 0, arr.ind = TRUE)
NLinks <- nrow(HP.ones)

SQres.jackn <- matrix(NA, nrow = NLinks, ncol = NLinks)
colnames(SQres.jackn) <- paste(rownames(int)[HP.ones[,1]], colnames(int)[HP.ones[,2]], sep = "-")

for (i in 1:NLinks) {
  int_jack <- as.matrix(int)
  int_jack[HP.ones[i,1], HP.ones[i,2]] <- 0
  D_jack <- prepare_paco_data(htree, ptree, int_jack)
  D_jack <- add_pcoord(D_jack, correction = "cailliez")
  D_jack <- PACo(D_jack, nperm = 1000, method = "r0", symmetric = TRUE, shuffled = TRUE)
  res_jack <- residuals_paco(D_jack$proc)^2
  res_jack <- append(res_jack, NA, after = i-1)
  SQres.jackn[i, ] <- res_jack
}

SQres <- res_sq * NLinks
SQres.jackn <- t(apply(SQres.jackn * -(NLinks - 1), 1, "+", SQres))

phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE)
phi.SD <- apply(SQres.jackn, 2, sd, na.rm = TRUE)
t.critical <- qt(0.975, df = NLinks - 1)

jack_df <- data.frame(
  link = colnames(SQres.jackn),
  contribution = phi.mean,
  LCI = phi.mean - t.critical * phi.SD / sqrt(NLinks),
  UCI = phi.mean + t.critical * phi.SD / sqrt(NLinks)
)

# Merge with residuals
plot_df <- merge(res_df, jack_df, by = "link")
plot_df$UCI_rescaled <- plot_df$jack + (plot_df$UCI - plot_df$contribution)
plot_df$LCI_rescaled <- plot_df$jack + (plot_df$LCI - plot_df$contribution)

# Save PACo scores
dir.create(paste0(output, "paco/"), recursive = TRUE, showWarnings = FALSE)
write.csv(plot_df, paste0(output, "paco/paco_scores.csv"))

#### Plot link contributions ----
dir.create(paste0(output, "figures/"), recursive = TRUE, showWarnings = FALSE)

svg(paste0(output, "figures/link_contribution.svg"), width = 16, height = 8)
ggplot(plot_df, aes(x = reorder(link, jack), y = jack)) +
  geom_bar(stat = "identity", fill = "white", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = jack, ymax = UCI_rescaled), width = 0.2) +
  geom_hline(yintercept = median(plot_df$res, na.rm = TRUE), linetype = "dashed", color = "red") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  labs(x = "Host-parasite link", y = "Squared residuals (PACo contribution)")
dev.off()

#### Tanglegram weighted by PACo residuals ----
assoc <- data.frame(
  greg = rownames(int)[which(int == 1, arr.ind = TRUE)[,'row']],
  hos = colnames(int)[which(int == 1, arr.ind = TRUE)[,'col']]
)
weight <- (res^-2) / 70

a <- cophylo(host_tree, parasite_tree, assoc, rotate = TRUE)

svg(paste0(output, "figures/cophylogeny.svg"), width = 16, height = 17)
plot(a, link.type = "curved", link.lwd = weight, link.col = make.transparent("blue", 0.5))
dev.off()

#### Specific analysis: Boliviana vs others ----
cophy_int <- grep('Boliviana', names(res))
noncophy <- res[-cophy_int]
cophy <- res[cophy_int]

ttest <- t.test(cophy, noncophy)
print(ttest)

dat <- rbind(
  data.frame(cophy = cophy, level = 'Boliviana'),
  data.frame(cophy = noncophy, level = 'Other')
)

svg(paste0(output, "figures/BolivianaVsRest_boxplot.svg"), width = 16, height = 17)
ggplot(dat, aes(x = level, y = cophy, fill = level)) +
  geom_boxplot(alpha = 0.85) +
  theme_bw() +
  scale_fill_brewer(palette = 'Paired') +
  ylab('Procrustes residual') +
  theme(legend.position = "none")
dev.off()

#### Density plot: Boliviana residuals vs background ----
svg(paste0(output, "figures/procrustes_residuals_Boliviana.svg"), width = 17, height = 10)
ggplot(data.frame(res = noncophy), aes(x = res)) +
  geom_density(fill = 'grey70') +
  geom_vline(data = data.frame(res = cophy), aes(xintercept = res), color = 'darkorange1') +
  theme_bw() +
  xlab('Procrustes residuals') +
  ylab('Frequency')
dev.off()

#### PACo m² null distribution plot ----
null <- as.data.frame(D$shuffled)
m2 <- as.data.frame(D$gof$ss)

svg(paste0(output, "figures/procrustes_residuals_significance.svg"), width = 15, height = 9)
ggplot(null, aes(x = `(D[["shuffled"]])`)) +
  geom_density(fill = 'grey70') +
  geom_vline(data = m2, aes(xintercept = `D$gof$ss`), col = 'darkorange1', linewidth = 1) +
  theme_bw() +
  xlab('Procrustes sum of squared residuals') +
  ylab('Frequency')
dev.off()

#### ParaFit analysis ----
z <- parafit(htree, ptree, int, nperm = 999, test.links = TRUE, seed = 123, correction = "cailliez")

link_significance <- data.frame(z$link.table)
link_significance$host_name <- rownames(int)[link_significance$Host]
link_significance$parasite_name <- colnames(int)[link_significance$Parasite]

link_significance$signif1_0.05 <- link_significance$p.F1 <= 0.05
link_significance$signif2_0.05 <- link_significance$p.F2 <= 0.05

dir.create(paste0(output, "parafit/"), recursive = TRUE, showWarnings = FALSE)
write.csv(link_significance, paste0(output, "parafit/link_significance.csv"))
