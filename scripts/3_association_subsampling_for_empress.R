## This code performs subsamplings on gregarine-Romalea associations to prepare input to be used in emPREss

# Creator: Jorge Medina-Duran

library(data.table)
library(dplyr) 
library(digest)


# Set working directory
input <- "C:/Users/beto2/Documents/cophylogeny_project/input/"

output <- "C:/Users/beto2/Documents/cophylogeny_project/output/"




## generate one-to-one networks
network <- read.table(paste0(output, "ocurrence_matrix/ocurrence_matrix_cophylogeny_renamed.csv"), sep=",", stringsAsFactors = F)
network2 <- network
rownames(network2) <- network[,1]
colnames(network2) <- network[1,]
network2 <- network2[-1,-1]


# Function to subsample one host per parasite
subsample_matrix <- function(matrix) {
  # Create an empty matrix with the same dimensions but filled with zeros
  subsampled_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  rownames(subsampled_matrix) <- rownames(matrix)
  colnames(subsampled_matrix) <- colnames(matrix)
  
  for (parasite in colnames(matrix)) {
    associated_hosts <- rownames(matrix)[matrix[, parasite] > 0]
    
    if (length(associated_hosts) > 0) {  # If the parasite has at least one host
      chosen_host <- sample(associated_hosts, 1)  # Randomly select one host
      subsampled_matrix[chosen_host, parasite] <- 1
    }
  }
  
  return(as.data.frame(subsampled_matrix))
}

# Generate 50 subsampled matrices
set.seed(123)  # For reproducibility
subsampled_matrices <- lapply(1:50, function(i) subsample_matrix(network2))


# Detect duplicates using hashing
matrix_hashes <- sapply(subsampled_matrices, function(mat) digest(mat, algo = "md5"))
duplicate_indices <- which(duplicated(matrix_hashes))

# Print duplicate matrices
if (length(duplicate_indices) > 0) {
  cat("Duplicate matrices detected at indices:\n", paste(duplicate_indices, collapse = ", "), "\n")
} else {
  cat("No duplicate matrices detected.\n")
}





# Save each matrix as a CSV file
subsample_dir <- paste0(output, "/subsampled_matrix/")
if(!dir.exists(subsample_dir)){
dir.create(subsample_dir, showWarnings = FALSE)
}

for (i in 1:50) {
  write.csv(subsampled_matrices[[i]], paste0(subsample_dir, "subsampled_matrix_", i, ".csv"), row.names = TRUE)
}

# Print first subsampled matrix as example
print(subsampled_matrices[[1]])






### Convert format to emPress ####

# Function to convert matrix to parasite-host format
convert_association_format <- function(matrix) {
  associations <- c()  # Initialize empty vector
  
  for (parasite in colnames(matrix)) {
    host <- rownames(matrix)[matrix[, parasite] > 0]  # Find associated host
    
    if (length(host) > 0) {  # Ensure there's an association
      associations <- c(associations, paste0(parasite,":",host))
    }
  }
  
  return(associations)
}

# Convert all 50 matrices into this format
subsampled_associations <- lapply(subsampled_matrices, convert_association_format)

# Save each association list to a text file (Optional)
for (i in 1:50) {
  writeLines(subsampled_associations[[i]], paste0(subsample_dir, "associations_", i, ".mapping"))
}

# Print the first few associations of the first matrix
print(subsampled_associations[[1]])







