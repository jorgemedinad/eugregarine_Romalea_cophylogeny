#############################################################
# Title: Generate Subsampled One-to-One Associations for eMPRess
# Creator: Jorge Medina-Duran
# Date: Apr/27/2025
#
# Description:
# This script generates random one-to-one subsampled networks from Romalea–gregarine associations.
# It prepares 50 matrices and corresponding eMPRess-compatible mappings.
#
# Output:
# - 50 subsampled matrices (CSV)
# - 50 association mappings (eMPRess format)
#
# Note:
# - Change only `input` and `output` folder paths below.
# - Results are saved into `/subsampled_matrix` subdirectory.
#############################################################

# Load required libraries
library(dplyr)
library(digest)

#### Define input and output directories ----
input <- "YOUR_PATH/input/"
output <- "YOUR_PATH/output/"

#### Load occurrence matrix ----
network <- read.csv(paste0(output, "occurrence_matrix/occurrence_matrix_cophylogeny_renamed.csv"), row.names = 1)
network2 <- as.matrix(network)

#### Subsampling function ----
subsample_matrix <- function(matrix) {
  subsampled_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  rownames(subsampled_matrix) <- rownames(matrix)
  colnames(subsampled_matrix) <- colnames(matrix)
  
  for (parasite in colnames(matrix)) {
    associated_hosts <- rownames(matrix)[matrix[, parasite] > 0]
    if (length(associated_hosts) > 0) {
      chosen_host <- sample(associated_hosts, 1)
      subsampled_matrix[chosen_host, parasite] <- 1
    }
  }
  
  return(as.data.frame(subsampled_matrix))
}

#### Generate 50 subsampled matrices ----
set.seed(123)  # For reproducibility
subsampled_matrices <- lapply(1:50, function(i) subsample_matrix(network2))

#### Check for duplicate matrices ----
matrix_hashes <- sapply(subsampled_matrices, function(mat) digest(mat, algo = "md5"))
duplicate_indices <- which(duplicated(matrix_hashes))

if (length(duplicate_indices) > 0) {
  cat("Duplicate matrices detected at indices:", paste(duplicate_indices, collapse = ", "), "\n")
} else {
  cat("No duplicate matrices detected.\n")
}

#### Save subsampled matrices ----
subsample_dir <- paste0(output, "subsampled_matrix/")
if (!dir.exists(subsample_dir)) {
  dir.create(subsample_dir, recursive = TRUE)
}

for (i in seq_along(subsampled_matrices)) {
  write.csv(subsampled_matrices[[i]], paste0(subsample_dir, "subsampled_matrix_", i, ".csv"), row.names = TRUE)
}

#### Convert matrices to eMPRess format ----
convert_association_format <- function(matrix) {
  associations <- c()
  for (parasite in colnames(matrix)) {
    host <- rownames(matrix)[matrix[, parasite] > 0]
    if (length(host) > 0) {
      associations <- c(associations, paste0(parasite, ":", host))
    }
  }
  return(associations)
}

# Generate and save mappings
subsampled_associations <- lapply(subsampled_matrices, convert_association_format)

for (i in seq_along(subsampled_associations)) {
  writeLines(subsampled_associations[[i]], paste0(subsample_dir, "associations_", i, ".mapping"))
}

#### Print example outputs ----
cat("Example subsampled matrix:\n")
print(subsampled_matrices[[1]])

cat("\nExample associations mapping:\n")
print(subsampled_associations[[1]])


