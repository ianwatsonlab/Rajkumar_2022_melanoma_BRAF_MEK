# Code to reproduce Figure 1B from Rajkumar et al. 2022
# Author: Mathieu Lajoie (mathieu.lajoie2@mcgill.ca)
# February 18, 2022

# Load required packages ####
require(checkmate)
require(ggplot2)
require(data.table)
require(scales)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions ####
source("R/Functions_cooc_analysis.R")

# Load mutation data ####
file.mutation.table = "Data/mutation_table.tsv"
mt = data.table::fread(file.mutation.table)

# Figure 1B - Mutation probability as function of TMB ####

win = 201 # window size
pc = 1 # pseudocount

colnames(mt)
sel.mut = c("BRAF.Class.I","BRAF.Class.II","BRAF.Class.III","NRAS.G12.G13.Q61","NF1.LOF")

# Get mutation probabilities
for(sel in sel.mut){
  new.col = paste0(sel,".prob")
  mt[[new.col]] = mavg(mt[[sel]], mt$TMB, pc = pc, k = win)
  print(new.col)
}

# Melt for plot
sel.cols.prob = paste0(sel.mut,".prob")
sel.cols = c("TMB", sel.cols.prob)
dm = melt(mt[,..sel.cols], measure = list(sel.cols.prob), value.name = "prob")
dm$mutation = gsub(".prob$","",dm$variable) # mutation label
dm$mutation = factor(dm$mutation, levels = sort(unique(as.character(dm$mutation))))

colors = c("#00BF7D","#F8766D","#A3A500","#00B0F6","#E76BF3")

# Plot
gg.prob = ggplot(dm, aes(x=TMB, y=prob)) +
  geom_line(aes(color = mutation), lwd=0.75) + 
  coord_cartesian(x = c(1, 2000)) +
  ylab("Mutation probability") + xlab("Total number of mutations in tumor") + 
  scale_x_continuous(breaks = c(1, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = colors) + 
  theme_bw(base_size = 8) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
gg.prob

outfile.prob = paste0("Figures/Fig1B_mutation_probabilities_w", win, ".pdf"); outfile.prob
ggsave(filename = outfile.prob, gg.prob, width = 12, height = 8, units = "cm")
system(paste("open",outfile.prob))


