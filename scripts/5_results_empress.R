#############################################################
# Title: Analyze Empress Reconciliations and Generate Event Visualizations
# Creator: Jorge Medina-Duran
# Date: Apr/27/2025
#
# Description:
# This script summarizes the outputs of Empress reconciliations and p-value tests.
# It performs:
# - Extracts p-values and event types
# - Computes event counts, time inconsistencies, and ratios
# - Saves summary tables
# - Generates event and ratio visualizations
#############################################################

# Load required libraries
library(ape)
library(ggplot2)
library(tidyr)

#### Define paths ----
output <- "YOUR_PATH/output/"
results_dir <- paste0(output, "results_empress/")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
# setwd(results_dir)  # Optional: usually better to work without changing working directory

#### Define cost schemes and parameters ----
cost_schemes <- c("d1_t1_l1", "d1_t8_l1", "d2_t1_l2", "d2_t4_l1", "d4_t1_l1")
n_replicates <- 50

# Load trees
parasite_tree <- read.tree(paste0(output, "trimmed_tree/gregarine_tree_trim.tre"))
host_tree <- read.tree(paste0(output, "trimmed_tree/romalea_tree_trim.tre"))

# Compute node ages
node_age_hosts <- node.depth.edgelength(host_tree)
node_age_hosts <- max(node_age_hosts) - node_age_hosts
names(node_age_hosts) <- 1:length(node_age_hosts)

#### Initialize result data frames ----
res_stats <- data.frame()
all_pvalues <- data.frame()
event_data <- data.frame()

#### Process reconciliations and p-values ----
for (costs in cost_schemes) {
  list_pvalues <- c()
  print(paste("Processing cost scheme:", costs))
  
  for (rep in 1:n_replicates) {
    recon_file <- paste0(results_dir, "reconciliation_", costs, "_associations_", rep, ".csv")
    res_file <- paste0(results_dir, "res_", costs, "_associations_", rep, ".svg")
    
    if (file.exists(res_file) && file.exists(recon_file)) {
      # Extract p-value from SVG
      res <- read.table(res_file, comment.char = "", fill = TRUE, sep = ";")
      res_pvalue <- as.numeric(gsub(" -->", "", gsub("    <!-- p-value = ", "", res$V1[grep("p-value", res$V1)])))
      
      # Load reconciliation table
      recon <- read.csv(recon_file, stringsAsFactors = FALSE)
      colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
      
      # Sanity check on frequencies
      if (!all(recon$node_frequency >= recon$event_frequency)) {
        warning("Problem in event frequencies!")
      }
      
      # Rename parasite and host nodes
      recon$parasite <- sapply(recon$parasite, function(x) {
        if (grepl("^_p", x)) as.numeric(gsub("_p", "", x)) + Ntip(parasite_tree) + 1
        else which(parasite_tree$tip.label == x)
      })
      recon$host <- sapply(recon$host, function(x) {
        if (grepl("^_h", x)) as.numeric(gsub("_h", "", x)) + Ntip(host_tree) + 1
        else which(host_tree$tip.label == x)
      })
      
      # Compute time consistency
      recon$age_min <- recon$age_max <- NA
      for (i in seq_len(nrow(recon))) {
        if (recon$event[i] == "Contemporaneous") {
          recon$age_min[i] <- 0
          recon$age_max[i] <- 0
        }
        if (recon$event[i] == "Cospeciation") {
          recon$age_min[i] <- node_age_hosts[recon$host[i]]
          recon$age_max[i] <- node_age_hosts[recon$host[i]]
        }
        if (recon$event[i] %in% c("Transfer", "Loss", "Duplication")) {
          previous_node <- which(parasite_tree$edge[,2] == recon$parasite[i])
          recon$age_max[i] <- if (length(previous_node) == 1) {
            max(recon$age_max[which(recon$parasite == parasite_tree$edge[previous_node,1])])
          } else {
            max(node_age_hosts)
          }
          recon$age_min[i] <- node_age_hosts[recon$host[i]]
        }
      }
      recon$diff_age <- recon$age_max - recon$age_min
      
      # Count events
      all_events <- c("Contemporaneous", "Cospeciation", "Transfer", "Duplication", "Loss")
      event_counts <- table(factor(recon$event, levels = all_events))
      
      # Store event data
      list_pvalues <- rbind(list_pvalues, c(rep, res_pvalue, as.numeric(event_counts), sum(recon$diff_age < 0)))
      
      event_data <- rbind(event_data, data.frame(
        cost_scheme = costs,
        replicate = rep,
        pvalue = res_pvalue,
        cospeciation = event_counts["Cospeciation"],
        transfer = event_counts["Transfer"],
        duplication = event_counts["Duplication"],
        loss = event_counts["Loss"],
        time_inconsistency = sum(recon$diff_age < 0)
      ))
    }
  }
  
  # Save replicate-level results
  list_pvalues <- as.data.frame(list_pvalues)
  colnames(list_pvalues) <- c("replicate", "pvalue", "Contemporaneous", "Cospeciation", "Transfer", "Duplication", "Loss", "time_inconsistency")
  write.table(list_pvalues, paste0(results_dir, "results_emPress_", costs, ".csv"), sep = ",", row.names = FALSE)
  
  # Summarize overall
  res_stats <- rbind(res_stats, data.frame(
    cost_scheme = costs,
    mean_pvalue = mean(as.numeric(list_pvalues$pvalue), na.rm = TRUE),
    total_cospeciation = sum(as.numeric(list_pvalues$Cospeciation)),
    total_transfer = sum(as.numeric(list_pvalues$Transfer)),
    total_duplication = sum(as.numeric(list_pvalues$Duplication)),
    total_loss = sum(as.numeric(list_pvalues$Loss)),
    total_time_inconsistency = sum(as.numeric(list_pvalues$time_inconsistency))
  ))
  
  all_pvalues <- rbind(all_pvalues, data.frame(cost_scheme = costs, pvalue = list_pvalues$pvalue))
}

# Save overall stats
write.table(res_stats, paste0(results_dir, "results_stats_emPress.csv"), sep = ",", row.names = FALSE)
write.table(all_pvalues, paste0(results_dir, "pvalue_distribution.csv"), sep = ",", row.names = FALSE)

#### Prepare event data for plotting ----
event_data$signif <- event_data$pvalue <= 0.05
event_data$ratio <- event_data$transfer / event_data$cospeciation
event_data_long <- pivot_longer(event_data, cols = c(cospeciation, transfer, duplication, loss),
                                names_to = "event", values_to = "number_events")
event_data_long$event <- factor(event_data_long$event, levels = c("cospeciation", "transfer", "duplication", "loss"))

write.table(event_data, paste0(results_dir, "event_distribution.csv"), sep = ",", row.names = FALSE)
write.table(event_data_long, paste0(results_dir, "event_distribution_long.csv"), sep = ",", row.names = FALSE)

#### Generate figures ----
figures_dir <- paste0(output, "figures/")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Significant event distributions
svg(paste0(figures_dir, "empress_event_distributions_significant.svg"), width = 7, height = 5)
ggplot(event_data_long[event_data_long$signif == TRUE,], aes(x = event, y = number_events, group = interaction(cost_scheme, replicate), color = cost_scheme)) +
  geom_line(alpha = 0.6, size = 0.5) +
  geom_hline(yintercept = 0) +
  geom_boxplot(aes(group = event), fill = "gray80", color = "black", width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  xlab("") + ylab("Number of events") +
  scale_color_manual(values = c("#ba4a00", "#d68910", "#229954", "#2e86c1", "#7a1fa2"))
dev.off()

# Non-significant event distributions
svg(paste0(figures_dir, "empress_event_distributions_nonsignificant.svg"), width = 7, height = 5)
ggplot(event_data_long[event_data_long$signif == FALSE,], aes(x = event, y = number_events, group = interaction(cost_scheme, replicate), color = cost_scheme)) +
  geom_line(alpha = 0.6, size = 0.5) +
  geom_hline(yintercept = 0) +
  geom_boxplot(aes(group = event), fill = "gray80", color = "black", width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  xlab("") + ylab("Number of events") +
  scale_color_manual(values = c("#ba4a00", "#d68910", "#229954", "#2e86c1", "#7a1fa2"))
dev.off()

# Transfer/Cospeciation ratios
svg(paste0(figures_dir, "empress_ratio_transfer_cosp.svg"), width = 7, height = 5)
ggplot(event_data, aes(x = cost_scheme, y = ratio, fill = signif)) +
  geom_violin(trim = TRUE, draw_quantiles = 0.5, alpha = 0.6, width = 1.3) +
  geom_point(aes(color = signif), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 1.3), size = 0.8, alpha = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  scale_y_continuous(trans = "sqrt", breaks = c(0, 0.25, 1, 3, 5, 9)) +
  xlab("Cost scheme") +
  ylab("Ratio transfers vs cospeciations") +
  scale_fill_manual(values = c("TRUE" = "#1a5276", "FALSE" = "#b03a2e")) +
  theme(legend.title = element_blank(), legend.position = "right")
dev.off()



