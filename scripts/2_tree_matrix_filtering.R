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
# - Figures and results are saved into `/occurrence_matrix`, and `/trimmed_tree` subdirectories.
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

# Load phylogenetic trees
gregarine.tree <- read.nexus(paste0(input, "trees/gregarines_dummy_dated.tre"))
romalea.tree <- read.nexus(paste0(input, "trees/romalea_dummy_dated.tre"))

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
  "microptera_5", "auricornis_1", "auricornis_2",
  "centurio", "eques", "obscura_1", "obscura_2", "picticornis",
  "reticulata", "tamaulipensis", "varipennis", "Xyleus"
)

romalea.tree.trim$tip.label <- new.names.h.tree

# Adjust occurrence matrix accordingly
rownames(occurrence) <- new.names.h.tree

#### Create output directories if needed ----
if (!dir.exists(paste0(output, "occurrence_matrix/"))) {
  dir.create(paste0(output, "occurrence_matrix/"), recursive = TRUE)
}
if (!dir.exists(paste0(output, "trimmed_tree/"))) {
  dir.create(paste0(output, "trimmed_tree/"), recursive = TRUE)
}

#### Save processed objects ----
# Occurrence matrix
saveRDS(occurrence, paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny_renamed.rds"))
write.csv(occurrence, paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny_renamed.csv"))

# Trimmed trees
write.tree(romalea.tree.trim, paste0(output, "trimmed_tree/romalea_tree_trim.tre"), tree.names = FALSE)
write.tree(gregarine.tree.trim, paste0(output, "trimmed_tree/gregarine_tree_trim.tre"), append = TRUE)

