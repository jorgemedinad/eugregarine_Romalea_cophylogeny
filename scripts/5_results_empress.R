
library(ape)
library(ggplot2)
library(tidyr)

# Set working directory
output <- "C:/Users/beto2/Documents/cophylogeny_project/output/"
results_dir <- paste0(output, "/results_empress")
setwd(results_dir)

# Define cost schemes and number of replicates
cost_schemes <- c("d1_t1_l1", "d1_t8_l1", "d2_t1_l2", "d2_t4_l1", "d4_t1_l1")
n_replicates <- 50  # Fixed at 50 replicates

# Load host and parasite trees (single phylogeny for all replicates)
parasite_tree <- read.tree(paste0(output, "trimmed_tree/gregarine_tree_trim.tre"))
host_tree <- read.tree(paste0(output, "trimmed_tree/romalea_tree_trim.tre"))

# Compute node ages
node_age_hosts <- node.depth.edgelength(host_tree)
node_age_hosts <- max(node_age_hosts) - node_age_hosts
names(node_age_hosts) <- 1:length(node_age_hosts)

res_stats <- data.frame()
all_pvalues <- data.frame()
event_data <- data.frame()

for (costs in cost_schemes) {
  list_pvalues <- c()
  print(costs)
  
  for (rep in 1:n_replicates) {
    recon_file <- paste0("reconciliation_", costs, "_associations_", rep, ".csv")
    res_file <- paste0("res_", costs, "_associations_", rep, ".svg")
    
    if (file.exists(res_file) && file.exists(recon_file)) {
      res <- read.table(res_file, comment.char = "", fill=TRUE, sep=";")
      
      res_pvalue <- res$V1[grep("p-value", res$V1)]
      res_pvalue <- gsub("    <!-- p-value = ", "", res_pvalue)
      res_pvalue <- as.numeric(gsub(" -->", "", res_pvalue))
      
      recon <- read.csv(recon_file, stringsAsFactors=FALSE)
      colnames(recon) <- c("parasite", "host", "event", "node_frequency", "event_frequency")
      
      if (!all(recon$node_frequency >= recon$event_frequency)) {
        print("problem in event frequencies")
      }
      
      # Convert node names (_pX, _hX) to numeric indices
      for (i in 1:nrow(recon)) {
        if (grepl("^_p", recon$parasite[i])) {
          recon$parasite[i] <- as.numeric(gsub("_p", "", recon$parasite[i])) + Ntip(parasite_tree) + 1
        } else {
          recon$parasite[i] <- which(parasite_tree$tip.label == recon$parasite[i])
        }
        
        if (grepl("^_h", recon$host[i])) {
          recon$host[i] <- as.numeric(gsub("_h", "", recon$host[i])) + Ntip(host_tree) + 1
        } else {
          recon$host[i] <- which(host_tree$tip.label == recon$host[i])
        }
      }
      
      # Time consistency checks
      recon$age_min <- NA
      recon$age_max <- NA
      
      for (i in 1:nrow(recon)) {
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
          if (length(previous_node) == 1) {
            recon$age_max[i] <- max(recon$age_max[which(recon$parasite == parasite_tree$edge[previous_node,1])])
          } else {
            recon$age_max[i] <- max(node_age_hosts)
          }
          recon$age_min[i] <- node_age_hosts[recon$host[i]]
        }
      }
      
      recon$diff_age <- recon$age_max - recon$age_min
      
      # Ensure all event types exist
      all_events <- c("Contemporaneous", "Cospeciation", "Transfer", "Duplication", "Loss")
      event_counts <- table(recon$event)
      event_counts <- as.numeric(event_counts[all_events])
      names(event_counts) <- all_events
      
      # Replace NA values with 0 (for missing event types)
      event_counts[is.na(event_counts)] <- 0
      
      # Store results correctly
      list_pvalues <- rbind(list_pvalues, c(rep, res_pvalue, event_counts, length(which(recon$diff_age < 0))))
      
      event_data <- rbind(event_data, data.frame(cost_scheme=costs, replicate=rep, pvalue=res_pvalue,
                                                 cospeciation=event_counts["Cospeciation"],
                                                 transfer=event_counts["Transfer"],
                                                 duplication=event_counts["Duplication"],
                                                 loss=event_counts["Loss"],
                                                 time_inconsistency=length(which(recon$diff_age < 0))))
    }
  }
  
  list_pvalues <- data.frame(list_pvalues)
  colnames(list_pvalues) <- c("replicate", "pvalue", "Contemporaneous", "Cospeciation", "Transfer", "Duplication", "Loss", "time_inconsistency")
  write.table(list_pvalues, paste0(results_dir, "/results_emPress_", costs, ".csv"), sep=",", row.names=FALSE)
  
  # Summarize statistics for this cost scheme
  res_stats <- rbind(res_stats, data.frame(cost_scheme=costs,
                                           mean_pvalue=mean(as.numeric(list_pvalues$pvalue), na.rm=TRUE),
                                           total_cospeciation=sum(as.numeric(list_pvalues$Cospeciation), na.rm=TRUE),
                                           total_transfer=sum(as.numeric(list_pvalues$Transfer), na.rm=TRUE),
                                           total_duplication=sum(as.numeric(list_pvalues$Duplication), na.rm=TRUE),
                                           total_loss=sum(as.numeric(list_pvalues$Loss), na.rm=TRUE),
                                           total_time_inconsistency=sum(as.numeric(list_pvalues$time_inconsistency), na.rm=TRUE)))
  
  all_pvalues <- rbind(all_pvalues, data.frame(cost_scheme=costs, pvalue=list_pvalues$pvalue))
}

colnames(res_stats) <- c("cost_scheme", "mean_pvalue", "total_cospeciation", "total_transfer", "total_duplication", "total_loss", "total_time_inconsistency")
write.table(res_stats, paste0(results_dir, "/results_stats_emPress.csv"), sep=",", row.names=FALSE)

# Save p-value distributions
write.table(all_pvalues, paste0(results_dir, "/pvalue_distribution.csv"), sep=",", row.names=FALSE)


# Add new columns
event_data$signif <- event_data$pvalue <= 0.05
event_data$time_inconsistency <- NULL
rownames(event_data) <- NULL
event_data$transfer[is.na(event_data$transfer)] <- 0
event_data$ratio <- event_data$transfer/event_data$cospeciation
# Save event distribution data for visualization
write.table(event_data, paste0(results_dir, "/event_distribution.csv"), sep=",", row.names=FALSE)



event_data_long <- pivot_longer(event_data, 
                              cols = c(cospeciation, transfer, duplication, loss), 
                              names_to = "event", 
                              values_to = "number_events")


# Reorder event factor levels
event_data_long$event <- factor(event_data_long$event, levels = c("cospeciation", "transfer", "duplication", "loss"))

# Save event distribution data for visualization
write.table(event_data_long, paste0(results_dir, "/event_distribution_long.csv"), sep=",", row.names=FALSE)


if(!dir.exists(paste0(output, "figures/"))){
  dir.create(paste0(output, "figures/"), recursive = T)
}


svg(paste0(output, "figures/empress_event_distributions_significant.svg"), width = 7, height = 5)
# significant reconciliations
ggplot(event_data_long[event_data_long$signif == TRUE,], 
       aes(x = event, y = number_events, 
           group = interaction(cost_scheme, replicate),  # Group by both cost_scheme & replicate
           color = cost_scheme)) +
  xlab("") +
  ylab("Number of events") +
  geom_line(alpha = 0.6, size = 0.5) +  # Transparency for better visibility
  geom_hline(yintercept = 0) +
  theme_minimal() +
  scale_color_manual(values = c("#ba4a00", "#d68910", "#229954", "#2e86c1", "#7a1fa2")) +  # Ensure distinct colors for cost_schemes
  theme(legend.position = "right") +
  geom_boxplot(aes(group = event), fill = "gray80", color = "black", width = 0.1, outlier.shape = NA) +  # Gray boxplots with black borders
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
dev.off()



svg(paste0(output, "figures/empress_event_distributions_nonsignificant.svg"), width=7, height=5)
# non-significant reconciliations
ggplot(event_data_long[event_data_long$signif == FALSE,], 
       aes(x = event, y = number_events, 
           group = interaction(cost_scheme, replicate),  # Group by both cost_scheme & replicate
           color = cost_scheme)) +
  xlab("") +
  ylab("Number of events") +
  geom_line(alpha = 0.6, size = 0.5) +  # Transparency for better visibility
  geom_hline(yintercept = 0) +
  theme_minimal() +
  scale_color_manual(values = c("#ba4a00", "#d68910", "#229954", "#2e86c1", "#7a1fa2")) +  # Ensure distinct colors for cost_schemes
  theme(legend.position = "right") +
  geom_boxplot(aes(group = event), fill = "gray80", color = "black", width = 0.1, outlier.shape = NA) +  # Gray boxplots with black borders
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
dev.off()


svg(paste0(output, "figures/empress_ratio_transfer_cosp.svg"), width=7, height=5)
# ratio transfer/cospeciation
ggplot(event_data, aes(y = ratio, x = cost_scheme, fill = signif)) +
  geom_hline(yintercept = 1) +
  geom_violin(trim = TRUE, draw_quantiles = 0.5, alpha = 0.6, width = 1.3) +  # Slight transparency to see points
  geom_point(aes(color = signif), position = position_jitterdodge(jitter.width = 0.12, 
                                                                  dodge.width = 1.3), size = 0.8, alpha = 0.5, show.legend = FALSE) +  
  theme_bw() +
  xlab("Cost scheme") +
  ylab("Ratio transfers vs cospeciations") +
  scale_y_continuous(trans = 'sqrt', breaks = c(0, 0.25, 1, 3, 5, 9)) +
  scale_fill_manual(
    values = c("TRUE" = "#1a5276", "FALSE" = "#b03a2e"), # Blue for significant, red for non-significant
    labels = c("FALSE" = "Non-Significant", "TRUE" = "Significant")
  ) +
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "black") # Match point colors to violin plots
  ) +
  theme(legend.title = element_blank(), legend.position = "right")
dev.off()




