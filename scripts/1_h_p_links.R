#############################################################
# Title: FASTA Sequence Deduplication and Occurrence Matrix Construction
# Creator: Jorge Medina
# Date: Sep/26/2024
#
# Description:
# This script processes a FASTA file of gregarine sequences associated with *Romalea* grasshoppers.
# It performs:
# 1. Deduplication of identical FASTA sequences
# 2. Mapping of sequences to host collection data
# 3. Construction of an occurrence matrix linking unique parasite haplotypes to host tips
#
# Output:
# - A binary occurrence matrix (gregarine haplotypes × host tips) in RDS and CSV formats.
#
# Note:
# - Change `input` and `output` folder paths.
# - Figures and results are saved into `/occurrence_matrix` subdirectory.
#############################################################

# Load required libraries
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
library(dplyr)
library(tidyr)

#### Define input and output directories ----
input <- "YOUR_PATH/input/"
output <- "YOUR_PATH/output/"

#### Load and prepare databases ----
# Load FASTA sequences
fasta.file <- readDNAStringSet(paste0(input, "fasta_sequences_gregarines/gregarines_from_taeniopoda_updated_unaligned.fasta"))

# Load host tip names database
host.tip.names <- read.csv(paste0(input, "databases/host_tip_names.csv")) %>%
  select(host_species, county, host_tip_name)

# Load host collecting data and clean
romalea.collection <- read.csv(paste0(input, "databases/Romalea_collecting_data.csv")) %>%
  filter(hosts_species_match != "") %>%
  distinct(hosts_species_match, .keep_all = TRUE)

#### Deduplicate FASTA sequences ----
# Convert to data frame
sequences.df <- as.data.frame(fasta.file)

# Keep only unique sequences
unique.seq <- sequences.df %>% distinct(x)

# Preserve sequence names
unique.seq$seq_name <- rownames(unique.seq)

#### Find duplicated sequence names ----
sequences.df.dup <- as.data.frame(fasta.file)
sequences.df.dup$duplicated_names <- rownames(sequences.df.dup)

duplicated.groups <- sequences.df.dup %>%
  group_by(x) %>%
  summarize(duplicated_names = paste(duplicated_names, collapse = ", "), .groups = "drop")

#### Prepare merged FASTA metadata ----
merged.fasta <- merge(unique.seq, duplicated.groups, by = "x")

# Expand duplicate names into long format
merged.fasta.long <- merged.fasta %>%
  separate_rows(duplicated_names, sep = ",") %>%
  mutate(hosts_species_match = sub("^[^_]*_[^_]*_[^_]*_([^_]*).*$", "\\1", duplicated_names))

#### Match host species and collecting data ----
# Match romalea collection info
merged.fasta.long <- merged.fasta.long %>%
  left_join(romalea.collection, by = "hosts_species_match") %>%
  filter(!is.na(host_species))

# Match host tip names
merged.fasta.long <- merged.fasta.long %>%
  left_join(host.tip.names, by = c("host_species", "county"))

#### Create occurrence matrix ----
occurrence <- merged.fasta.long %>%
  select(seq_name, host_tip_name) %>%
  distinct()

# Create binary matrix
occurrence_matrix <- table(occurrence$seq_name, occurrence$host_tip_name)
occurrence_matrix_df <- as.data.frame.matrix(occurrence_matrix)

# Rename rows with standardized gregarine names
new.names.gregarines <- c(
  "Ataeniopoda_MK181531", "Amoebogregarina_7", "Amoebogregarina_8",
  "Amoebogregarina_3", "Amoebogregarina_1",
  "Amoebogregarina_6", "Amoebogregarina_5", "Amoebogregarina_2",
  "Boliviana_3", "Boliviana_1", "Boliviana_4", "Boliviana_2",
  "Coronoepimeritus_3", "Coronoepimeritus_2", "Coronoepimeritus_4",
  "Coronoepimeritus_5", "Coronoepimeritus_1", "Coronoepimeritus_6",
  "Coronoepimeritus_7", "Gregarina_2", "Gregarina_1", "Cmexicana_MK181532"
)
row.names(occurrence_matrix_df) <- new.names.gregarines

#### Save occurrence matrix ----
if (!dir.exists(paste0(output, "occurrence_matrix/"))) {
  dir.create(paste0(output, "occurrence_matrix/"), recursive = TRUE)
}

saveRDS(occurrence_matrix_df, paste0(output, "occurrence_matrix/occurrence_matrix.rds"))
write.csv(occurrence_matrix_df, paste0(output, "occurrence_matrix/occurrence_matrix.csv"))

