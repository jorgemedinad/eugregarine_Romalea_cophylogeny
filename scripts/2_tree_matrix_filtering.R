#############################################################
# Title: Prepare Data for Cophylogenetic Analysis (Romalea × Gregarines)
# Creator: Jorge Medina
# Date: Sep/26/2024
#
# Description:
# This script prepares and trims data for cophylogenetic analysis between *Romalea* grasshoppers and their eugregarine parasites.
# Outputs occurrence matrix and trimmed trees ready for global-fit or event-based analyses.
#
# Output:
# - Updated occurrence matrix (RDS and CSV)
# - Trimmed Romalea tree (NEXUS)
# - Trimmed Gregarine tree (NEXUS)
#
# Note:
# - Change only `input` and `output` folder paths below.
#############################################################

# Load required libraries
library(ape)  # for reading, trimming, and saving phylogenetic trees

#### Define input and output directories ----
input <- "YOUR_PATH/input/"
output <- "YOUR_PATH/output/"

#### Load data ----
# Load and transpose occurrence matrix
occurrence <- readRDS(paste0(output, "occurrence_matrix/occurrence_matrix.rds"))
occurrence <- t(as.matrix(occurrence))

# Remove duplicated parasite sequence
occurrence <- occurrence[, colnames(occurrence) != "Amoebogregarina_4"]

# Load phylogenetic trees
gregarine.tree <- read.nexus(paste0(input, "trees/good/gregarines_dummy_dated.tre"))
romalea.tree <- read.nexus(paste0(input, "trees/good/romalea_dummy_dated.tre"))

# Optional: visualize trees
# par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
# plot(gregarine.tree); nodelabels()
# plot(romalea.tree)

#### Prepare and trim trees ----
# Identify tips not present in occurrence matrix
trim.romalea <- setdiff(romalea.tree$tip.label, rownames(occurrence))
trim.gregarines <- setdiff(gregarine.tree$tip.label, colnames(occurrence))

# Prune trees
romalea.tree.trim <- drop.tip(romalea.tree, trim.romalea)
gregarine.tree.trim <- drop.tip(gregarine.tree, trim.gregarines)

#### Rename host tree tips ----
new.names.h.tree <- c(
  "microptera_3", "microptera_2", "microptera_1", "microptera_4",
  "microptera_5", "microptera_6", "auricornis_1", "auricornis_2",
  "centurio", "eques", "obscura_1", "obscura_2", "picticornis",
  "reticulata", "tamaulipensis", "varipennis", "Xyleus"
)

romalea.tree.trim$tip.label <- new.names.h.tree

# Drop specific problematic sample
romalea.tree.trim <- drop.tip(romalea.tree.trim, "microptera_5")

# Adjust occurrence matrix accordingly
rownames(occurrence) <- new.names.h.tree
occurrence <- occurrence[rownames(occurrence) != "microptera_5", ]

# Remove wrong link
occurrence["reticulata", "Boliviana_1"] <- 0

#### Create output directories if needed ----
if (!dir.exists(paste0(output, "occurrence_matrix/"))) {
  dir.create(paste0(output, "occurrence_matrix/"), recursive = TRUE)
}
if (!dir.exists(paste0(output, "trimmed_tree/"))) {
  dir.create(paste0(output, "trimmed_tree/"), recursive = TRUE)
}

#### Save processed objects ----
# Occurrence matrix
saveRDS(occurrence, paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny.rds"))
write.csv(occurrence, paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny_renamed.csv"))

# Trimmed trees
write.tree(romalea.tree.trim, paste0(output, "trimmed_tree/romalea_tree_trim.tre"), tree.names = FALSE)
write.tree(gregarine.tree.trim, paste0(output, "trimmed_tree/gregarine_tree_trim.tre"), append = TRUE)

