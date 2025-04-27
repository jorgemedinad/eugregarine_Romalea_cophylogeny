### Creator: Jorge Medina
### Date: Sep/26/2024


#### This code perform cophylogenetic analysis of gregarine 18S tree and Romalea grasshoppers

#### Install necessary packages
# install.packages("devtools")
library(devtools)
#devtools::install_github("mllaberia/Rtapas", force =T)
library(Rtapas)
#library(ape)


# Set working directory
input <- "C:/Users/beto2/Documents/cophylogeny_project/input/"

output <- "C:/Users/beto2/Documents/cophylogeny_project/output/"

#### Load data ----
# occurrence matrix
# occurrence <- read.csv("C:/Users/beto2/Google Drive/Research/DATA_dissertation/Chapter4_Taeniopoda/cophylogeny_project/
#                        output/ocurrence_matrix_cophylogeny_renamed.csv", header = T, row.names = 1)
occurrence <- readRDS(paste0(output, "occurrence_matrix/occurrence_matrix.rds"))
occurrence <- as.matrix(occurrence)
occurrence <- t(occurrence)
# Remove duplicated sequence
occurrence <- occurrence[,colnames(occurrence)!="Amoebogregarina_4"]
# gregarine phylogeny
par(mar=c(0,0,0,0), oma=c(0, 0, 0, 0))
gregarine.tree <- read.nexus(paste0(input, "trees/good/gregarines_dummy_dated.tre"))
plot(gregarine.tree)
nodelabels()

# Romalea tree
romalea.tree <- read.nexus(paste0(input, "trees/good/romalea_dummy_dated.tre"))
plot(romalea.tree)

#### Prepare data ----
# Identify tips to trim
trim.romalea <- dplyr::setdiff(romalea.tree$tip.label, rownames(occurrence))
trim.romalea

trim.gregarines <- dplyr::setdiff(gregarine.tree$tip.label, colnames(occurrence))
# trim.gregarines

# trim tips
romalea.tree.trim <- drop.tip(romalea.tree, trim.romalea)
gregarine.tree.trim <- drop.tip(gregarine.tree, trim.gregarines)


#### rename tree tips ----

new.names.h.tree <- c("microptera_3","microptera_2","microptera_1","microptera_4",
                      "microptera_5","microptera_6","auricornis_1","auricornis_2",
                      "centurio","eques","obscura_1","obscura_2","picticornis",
                      "reticulata","tamaulipensis","varipennis","Xyleus")

romalea.tree.trim$tip.label <- new.names.h.tree
romalea.tree.trim <- drop.tip(romalea.tree.trim, "microptera_5")


### Remove wrong links
rownames(occurrence) <- new.names.h.tree
occurrence <- occurrence[rownames(occurrence)!="microptera_5",]
occurrence["reticulata","Boliviana_1"] = 0



if(!dir.exists(paste0(output, "occurrence_matrix/"))){
  dir.create(paste0(output, "occurrence_matrix/"), recursive = T)
}

if(!dir.exists(paste0(output, "trimmed_tree/"))){
  dir.create(paste0(output, "trimmed_tree/"), recursive = T)
}



### Save objects for cophylogenetic analysis ----
# Save updated matrix
saveRDS(occurrence, paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny.rds"))
write.csv(occurrence, paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny_renamed.csv"))
# Romalea tree
write.tree(romalea.tree.trim, paste0(output, "trimmed_tree/romalea_tree_trim.tre"), tree.names = F)
# Gregarine tree
write.tree(gregarine.tree.trim, paste0(output, "trimmed_tree/gregarine_tree_trim.tre"), append = T)


