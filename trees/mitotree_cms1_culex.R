library(ape)
library(ggtree)

tree <- root(read.tree("seqs.fasta.raxml.bestTree"), 'CMS_012_RNA_A_S4')
metadata_all <- read.csv("../data/project-mosquito_sample-table_cleaned.csv")
metadata <- metadata_all[sapply(tree$tip.label, grep, metadata_all$sample_name), ]
sp_all <- read.csv("../../data/180806_project-mosquito_sample-table_pruned.csv", header=FALSE, stringsAsFactors = FALSE)
sp <- sp_all[sapply(tree$tip.label, grep, sp_all[, 1]), ]
df <- data.frame(seq=tree$tip.label, 
                 long=factor(round(metadata$long, 1)),
                 lat=factor(round(metadata$lat, 1)),
                 sp=factor(paste(sp[, 2], sp[, 3])))

tree_plot <- ggtree(tree) + theme_tree2() + xlim(0, .1)
tree_plot <- tree_plot %<+% df + geom_tiplab(aes(color=sp)) + 
  theme(legend.position="top", axis.line.x=element_line(size=1), axis.text.x=element_text(size=14),
        legend.text=element_text(size=14), legend.title=element_text(size=14))
tree_plot
ggsave("../figures/CMS001_culex_tree.pdf", tree_plot, width=5, height=4)

