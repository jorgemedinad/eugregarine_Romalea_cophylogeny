### Creator: Jorge Medina
### Date: Sep/26/2024


#### This code perform fasta sequence deduplication

#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
library(dplyr)
library(tidyr)

# Set working directory
input <- "C:/Users/beto2/Documents/cophylogeny_project/input/"

output <- "C:/Users/beto2/Documents/cophylogeny_project/output/"

#### Load and prepare databases ----
# sequences fasta

fasta.file <- readDNAStringSet(paste0(input, "fasta_sequences_gregarines/gregarines_from_taeniopoda_updated_unaligned.fasta"))  # Replace with your FASTA file path

# hosts_tips_names database
host.tip.names <- read.csv("C:/Users/beto2/Documents/cophylogeny_project/input/databases/host_tip_names.csv")
host.tip.names.filtered <- host.tip.names %>% select(host_species, county, host_tip_name)

# host collecting data
romalea.collection <- read.csv("C:/Users/beto2/Documents/cophylogeny_project/input/databases/Romalea_collecting_data.csv")
# Filter empty and repeated values of columns to match (hosts_ID_simplified)
romalea.collection <- romalea.collection %>%
  filter(hosts_species_match != "") %>%
  distinct(hosts_species_match, .keep_all = TRUE)




#### Create dataframe with unique FASTA sequences and names  ----
# Convert fasta files in a data frame
sequences.df <- as.data.frame(fasta.file)

# Remove duplicates
unique.seq <- sequences.df %>% distinct(x)

# Convert row names in a column
unique.seq$seq_name <- rownames(unique.seq)




#### Create dataframe with all the names of duplicated FASTA sequences ----
# Convert fasta files in a data frame
sequences.df.dup <- as.data.frame(fasta.file)

# Convert row names in a column
sequences.df.dup$duplicated_names <- rownames(sequences.df.dup)

# Get all the names of duplicated sequences
duplicated.groups <- sequences.df.dup %>% group_by(x) %>% summarize(duplicated_names = paste(duplicated_names, collapse = ", "))




#### Prepare data frame with FASTA to create occurrence matrix  ----
# Combine the two data frames
### THIS TABLE CONTAIN THE RELATIONSHIP BETWEEN GREGARINE TIP NAMES AND THE DUPLICATED GENE NAMES ----
merged.fasta <- merge(unique.seq, duplicated.groups, by = "x")

# Parse and pivot long values of column "duplicated_names"
merged.fasta.long <- merged.fasta %>%
  separate_rows(duplicated_names, sep = ",")

# Extract the host species code from duplicated_names column
merged.fasta.long$hosts_species_match <- sub("^[^_]*_[^_]*_[^_]*_([^_]*).*$", "\\1", merged.fasta.long$duplicated_names)




#### Match and merge romalea.collection database in merged.fasta.long dataframe
# Match databases
merged.fasta.long$position <- match(merged.fasta.long$hosts_species_match, romalea.collection$hosts_species_match)
# Remove NA values
merged.fasta.long <- merged.fasta.long %>% filter(!is.na(position))


# Create columns in merged.fasta.long dataframe
merged.fasta.long$host_species <- NA
merged.fasta.long$county <- NA

# Fill columns based on match position
merged.fasta.long$host_species <- romalea.collection$host_species[merged.fasta.long$position]
merged.fasta.long$county <- romalea.collection$county[merged.fasta.long$position]


#### Join merged.fasta.long with host.tip.names ----
# Join the two dataframes based on 'host_species' and 'county' columns
merged.fasta.long <- merged.fasta.long %>%
  left_join(host.tip.names.filtered, by = c("host_species", "county"))


#### Create occurrence matrix ----
#Select columns to create occurrence matrix
occurrence <- merged.fasta.long %>% select(seq_name, host_tip_name)

# Remove rows with duplicated combinations
occurrence <- occurrence %>%
  distinct(seq_name, host_tip_name)


# Create the occurrence matrix using table()
occurrence_matrix <- table(occurrence$seq_name, occurrence$host_tip_name)

# Convert occurrence matrix to a dataframe
occurrence_matrix_df <- as.data.frame.matrix(occurrence_matrix)
#View(occurrence_matrix_df)

new.names.gregarines <- c("Ataeniopoda_MK181531","Amoebogregarina_7","Amoebogregarina_8",
                        "Amoebogregarina_4","Amoebogregarina_3","Amoebogregarina_1",
                        "Amoebogregarina_6","Amoebogregarina_5","Amoebogregarina_2",
                        "Boliviana_3","Boliviana_1","Boliviana_4","Boliviana_2",
                        "Coronoepimeritus_3","Coronoepimeritus_2","Coronoepimeritus_4",
                        "Coronoepimeritus_5","Coronoepimeritus_1","Coronoepimeritus_6",
                        "Coronoepimeritus_7","Gregarina_2","Gregarina_1",
                        "Cmexicana_MK181532")

row.names(occurrence_matrix_df) <- new.names.gregarines

# new.names.host <- c("microptera_3","microptera_2","microptera_1","microptera_4",
#                         "microptera_5","microptera_6","auricornis_1","auricornis_2",
#                         "centurio","eques","obscura_1","obscura_2","picticornis",
#                         "reticulata","tamaulipensis","varipennis","Xyleus")
# colnames(occurrence_matrix_df) <- new.names.host


if(!dir.exists(paste0(output, "occurrence_matrix/"))){
  dir.create(paste0(output, "occurrence_matrix/"), recursive = T)
}

# Save occurrence matrix
saveRDS(occurrence_matrix_df, paste0(output, "occurrence_matrix/occurrence_matrix.rds"))
write.csv(occurrence_matrix_df, paste0(output, "occurrence_matrix/occurrence_matrix.csv"))
 

