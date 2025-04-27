
# This code performs cophylogenetic analysis using PACo and ParaFit
# 1. Loading phylogenies and association matrix
# 2. Perform PACo analysis
# 3. return interaction-specific cophylogenetic contributions
# 4. Performs Jackknife estimation of confidence intervals on cophylogenetic contributions
# 5. Plot link contribution to cophylogeny
# 6. Plot tanglegram with link contribution
# 7. Plot specific link contribution relative to other links (Bolviana against the rest)
#  7a. Welsh T-test Boliviana interactions against the rest
#  7b. Visualise residuals of Boliviana interactions against the rest
# 8. Plot m-squared significance value relative to null distribution
# 9. ParaFit significance test

library(ape)
library(paco)
library(phytools)
library(ggplot2)
library(tidyr)

# Set working directory
output <- "C:/Users/beto2/Documents/cophylogeny_project/output"

#### 1. Load host and parasite trees (single phylogeny for all replicates) ####
parasite_tree <- read.tree(paste0(output, "/trimmed_tree/gregarine_tree_trim.tre"))
plot(parasite_tree)

host_tree <- read.tree(paste0(output, "/trimmed_tree/romalea_tree_trim.tre"))
plot(host_tree)
# read links
int <- read.csv(paste0(output, "/occurrence_matrix/occurrence_matrix_cophylogeny_renamed.csv"), row.names=1)
# int <- int[rownames(int)!="microptera_5",]


# # Function to normalize patristic distance matrix
# normalize_max <- function(mat) {
#   return(mat / max(mat))
# }



# transform phylogenies to distance matrix
ptree <- cophenetic(parasite_tree)
htree <- cophenetic(host_tree)

# max_dist <- max(c(ptree, htree))
# 
# ptree <- ptree/max_dist
# htree <- htree/max_dist


# htree <- normalize_max(htree)
# ptree <- normalize_max(ptree)

  
#### PACo analysis ####
D <- prepare_paco_data(htree, ptree, int)
D <- add_pcoord(D, correction = "cailliez")
# perform cophylogenetic analysis
D <- PACo(D, nperm = 10000, seed = 13, method ="r0", symmetric = T, shuffled = T)
# Return interaction-specific cophylogenetic contributions based on a jackknifing procedure
D <- paco_links(D)

#### 3. return interaction-specific cophylogenetic contributions ####
res <- residuals_paco(D$proc)
res_df <- data.frame(res)
res_df$link <- row.names(res_df)
res_df$jack <- D$jackknife
res_sq <- res^2


# res_reordered <- res[order(res)]
# sum(res_reordered[1:13])
# mean(res_reordered)

#### 4. Jackknife estimation of confidence intervals ####
HP.ones <- which(as.matrix(int) > 0, arr.ind = TRUE)
NLinks <- nrow(HP.ones)

SQres.jackn <- matrix(NA, nrow = NLinks, ncol = NLinks)
colnames(SQres.jackn) <- paste(rownames(int)[HP.ones[,1]], colnames(int)[HP.ones[,2]], sep = "-")

for (i in 1:NLinks) {
  int_jack <- as.matrix(int)
  int_jack[HP.ones[i,1], HP.ones[i,2]] <- 0  # drop 1 link
  
  # Re-run PACo without the link
  D_jack <- prepare_paco_data(htree, ptree, int_jack)
  D_jack <- add_pcoord(D_jack, correction = "cailliez")
  D_jack <- PACo(D_jack, nperm = 1000, method = "r0", symmetric = TRUE, shuffled = TRUE)
  
  res_jack <- residuals_paco(D_jack$proc)^2
  res_jack <- append(res_jack, NA, after = i-1)
  SQres.jackn[i, ] <- res_jack
}


SQres <- res_sq * NLinks
SQres.jackn <- SQres.jackn * (-(NLinks - 1))
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres))

phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE)
phi.SD <- apply(SQres.jackn, 2, sd, na.rm = TRUE)
t.critical <- qt(0.975, df = NLinks - 1)

phi.UCI <- phi.mean + t.critical * phi.SD / sqrt(NLinks)
phi.LCI <- phi.mean - t.critical * phi.SD / sqrt(NLinks)

jack_df <- data.frame(
  link = colnames(SQres.jackn),
  contribution = phi.mean,
  LCI = phi.LCI,
  UCI = phi.UCI
)




#### 5. Plot the contributions and CI of individual links to Procrustean fit ####
plot_df <- merge(res_df, jack_df, by = "link")
res_sq <- res_sq[order(names(res_sq))]

# Add squared residuals to plot_df after the "res" column
plot_df <- plot_df %>%
  mutate(res_sq = as.numeric(res_sq)) %>%
  relocate(res_sq, .after = res)



plot_df$UCI_rescaled <- plot_df$jack + (plot_df$UCI - plot_df$contribution)
plot_df$LCI_rescaled <- plot_df$jack + (plot_df$LCI - plot_df$contribution)


svg(paste0(output, "/figures/link_contribution.svg"), width = 16, height = 8)
ggplot(plot_df, aes(x = reorder(link, jack), y = jack)) +
  geom_bar(stat = "identity", fill = "white", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = jack, ymax = UCI_rescaled), width = 0.2) +
  geom_hline(yintercept = median(plot_df$res, na.rm = TRUE), linetype = "dashed", color = "red") +
  coord_cartesian(clip = "off") +  # <<< key fix
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1, size = 13),
    axis.title = element_text(size = 14),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 45)
  ) +
  labs(
    x = "Host-parasite link",
    y = "Squared residuals (PACo contribution)"
  )

dev.off()


if(!dir.exists(paste0(output, "/paco"))){
  dir.create(paste0(output, "/paco"), recursive = T)
}

### Save PACo scores
write.csv(plot_df, paste0(output, "/paco/paco_scores.csv"))



#### 6. Plot tanglegram with link contribution using Ape function cophyloplot or similar ####
#(e.g., cophylo from phytools) weighted by interaction contribution.
# first we must make a list out of the interaction matrix
assoc <- data.frame(greg=rownames(int)[which(int==1, arr.ind=TRUE)[,'row']], hos=colnames(int)[which(int==1, arr.ind=TRUE)[,'col']])


# to weight the interactions we use the cophylogenetic contribution transformed to best show
# the differences graphically
weight <- (res^-2)/70
# weight <- -((-res^-2)*17)

# cophyloplot(host_tree, parasite_tree, assoc, show.tip.label=T, use.edge.length=T,
#             lwd=weight, col='steelblue', length.line=0, gap=27, space=180)
# png("c:/Users/beto2/Downloads/a.png", width = 1480, height = 1700, units = "px")
# cophyloplot(host_tree, parasite_tree, assoc, show.tip.label=T, use.edge.length=T,
#             lwd=weight, col='steelblue', length.line=0, gap=27, space=180)
# dev.off()

a <- cophylo(host_tree, parasite_tree, assoc, rotate=T)
plot(a, link.type="curved", link.lwd=weight,link.lty="solid",
     link.col=make.transparent("blue",0.5),cex=1,fsize=2)

svg(paste0(output, "/figures/cophylogeny.svg"), width = 16, height = 17)
plot(a, link.type="curved", link.lwd=weight,link.lty="solid",
     link.col=make.transparent("blue",0.5),cex=1,fsize=2)
dev.off()



#### 7. Plot specific link contribution relative to other links ####
# to analyse the links further we will split the interactions based on the plot
cophy_int <- c(grep('Boliviana', names(res)))
# grep('microptera', names(res)),

# remove the one interaction of Oxalis compacta that is not with the Lepidopterans
#cophy_int <- cophy_int[-grep('Anthidium', names(res[cophy_int]))]
noncophy <- res[-cophy_int]
cophy <- res[cophy_int]

##### 7a. Welch's t-test to test difference in cophylogenetic signal between different groups of interaction links #####
ttest <- t.test(cophy, noncophy)
ttest

# visualise the difference with a box and whisker plot
dat <- rbind(data.frame(cophy=cophy, level='high'), data.frame(cophy=noncophy, level='rest'))

svg(paste0(output, "/figures/BolivianaVsRest_boxplot.svg"), width = 16, height = 17)
ggplot(dat, aes(x=level, y=cophy, fill=level))+
  geom_boxplot(alpha=0.85)+
  scale_fill_brewer(palette='Paired')+
  scale_x_discrete(labels=c('Rm-B','Rest'))+
  ylab('Procrustes residual')+
  theme_bw()+
  theme(
    axis.title.x=element_blank(),
    panel.grid=element_blank(),
    legend.position='none',
    axis.text.x=element_text(size=14),
    axis.title.y=element_text(size=14)
  )
dev.off()


# test the influence of degree on cophylogenetic signal
special <- data.frame(cophy.sig=res, pla_deg=NA, row.names=names(res))
f <- function(x) length(grep(strsplit(x, '-')[[1]][2], names(res)))
special$pla_deg <- sapply(rownames(special), f)

greg_degree <- summary(lm(cophy.sig ~ pla_deg, data=special))

greg_degree


##### 7b. Visualise residuals of Boliviana interactions against the rest #####
noncophy_dat <- data.frame(res=noncophy)
cophy_dat <- data.frame(res=cophy)

svg(paste0(output, "/figures/procrustes_residuals_Boliviana.svg"), width = 17, height = 10)
ggplot(noncophy_dat, aes(x=res)) +
  geom_density(fill='grey70') +
  theme_bw()+
  geom_vline(data=cophy_dat, aes(xintercept=res), col='darkorange1') +
  theme(
    panel.grid=element_blank(),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    plot.margin = margin(50, 50, 50, 50, unit = "pt")
  ) +
  xlab('Procrustes residuals')+
  ylab('Frequency')+
  scale_x_continuous(limits=c(0, 0.6), expand=c(0.,0), breaks=c(0.25,0.5,0.75,1))+
  scale_y_continuous(expand=c(0,0), limits=c(0.00, 10))
dev.off()




#### 8. plot m-squared significance value relative to null distribution ####

null <- as.data.frame((D[["shuffled"]]))
m2 <- as.data.frame(D$gof$ss)

svg(paste0(output, "/figures/procrustes_residuals_significance.svg"), width = 15, height = 9)
ggplot(null, aes(x=`(D[["shuffled"]])`)) +
  geom_density(fill='grey70') +
  theme_bw() +
  geom_vline(data=m2, aes(xintercept=`D$gof$ss`), col='darkorange1', linewidth = 1) +
  theme(
    panel.grid=element_blank(),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    plot.margin = margin(50, 50, 50, 50, unit = "pt")
  ) +
  xlab('Procrustes sum of squared residuals')+
  ylab('Frequency')
dev.off()



#### 9. ParaFit significance test ####
z <- parafit(htree, ptree, int, nperm = 999, test.links = T, seed = 123, correction = "cailliez")

link_significance <- data.frame(z$link.table)
para_name <- colnames(int)
host_name <- row.names(int)
link_significance$host_name <- host_name[link_significance$Host]
link_significance$parasite_name <- para_name[link_significance$Parasite]


link_significance$signif1_0.05 <- link_significance$p.F1 <= 0.05
link_significance$signif1_0.02 <- link_significance$p.F1 <= 0.02

link_significance$signif2_0.05 <- link_significance$p.F2 <= 0.05
link_significance$signif2_0.02 <- link_significance$p.F2 <= 0.02


if(!dir.exists(paste0(output, "/parafit"))){
  dir.create(paste0(output, "/parafit"), recursive = T)
}


write.csv(link_significance, paste0(output, "/parafit/link_significance.csv"))

