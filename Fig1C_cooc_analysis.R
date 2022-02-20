# Code to reproduce Figure 1C from Rajkumar et al. 2022
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

# Figure 1C - Co-occurrence analysis ####

# Mutation pairs to test
PT = data.frame(m1 = c("BRAF.Class.III","BRAF.Class.II","BRAF.Class.I","NRAS.G12.G13.Q61"), m2 = c("NF1.LOF","NF1.LOF", "NF1.LOF","BRAF.Class.III"),stringsAsFactors = FALSE)
PT$pair = paste(PT$m1,"and",PT$m2)

#nb_sim = 1e4; nb.str = "1e4" # Quick test
nb_sim = 1e6; nb.str = "1e6" # Number of simulations
win = 201 ; win.str = "w201" # Window size used to estimate mutation probability

lres = list()
set.seed(1)

for(r in 1:nrow(PT)){
  print(paste("Processing", PT$pair[r], "..."))
  m1 = mt[[PT$m1[r]]]; sum(m1) # Mutation status
  m2 = mt[[PT$m2[r]]]; sum(m2)
  p1 = mavg(m1, mt$TMB, k = win, pc = 1);sum(p1) # Mutation probability (given TMB)
  p2 = mavg(m2, mt$TMB, k = win, pc = 1);sum(p2)
  dfr = mc_sim_cooc(m1 = m1, m2 = m2, p1 = p1, p2 = p2, nb_sim = nb_sim) # Call to sim function
  dfr$pair = PT$pair[r]
  title(PT$pair[r])
  lres[[PT$pair[r]]] <- dfr
}

dfres = as.data.frame(do.call(rbind, lres))

# Update labels
dfres$label = "observed"
dfres$label[dfres$test=="sim"] <- "expected"
dfres$label = factor(dfres$label,levels = c("observed","expected"))

# Format p-values
dfres$pval.str = paste("p =", signif(dfres$pval.two.tailed, digits = 2))
dfres$pval.str[is.na(dfres$pval.two.tailed)] = NA
dfres$pair2 = gsub("and","\nand\n", dfres$pair)
newlab = ""
dfres = dfres[order(as.character(dfres$pair2)),]
dfres$pair2 = factor(dfres$pair2, levels = unique(as.character(dfres$pair2))) # Specify x label ordering for plot
dfres$pair2

# Plot
gg.cooc = ggplot(dfres, aes(y = nb.cooc, x = pair2, fill = label)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  xlab(newlab) + theme_bw() + labs(fill = newlab) +
  geom_errorbar(aes(ymin = ci95.low, ymax = ci95.high, width=.1), position = position_dodge(width = 0.75)) + 
  scale_fill_manual(values = c("grey36","lightgrey")) + 
  geom_text(label = dfres$pval.str, nudge_y = 11, nudge_x =  0, size=2) + 
  ylim(c(0,40)) + 
  theme_bw(base_size = 8) +
  theme(panel.grid.major.x = element_blank()) + 
  ylab("Number of co-occuring mutations") +
  theme(legend.key.size = unit(0.4, 'cm'))
gg.cooc

outfile.cooc = paste0("Figures/Fig1C_mutation_cooc_tests_r",nb.str,"_",win.str,".pdf"); outfile.cooc
ggsave(filename = outfile.cooc, gg.cooc, width = 12, height = 8, units = "cm")
system(paste("open",outfile.cooc))
