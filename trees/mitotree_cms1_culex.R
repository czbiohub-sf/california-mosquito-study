library(ape)
library(ggtree)

# tree <- read.tree("seqs_COI.afa.raxml.bestTree")
# tree$tip.label <- gsub("_COI.fasta", "", tree$tip.label)
# tree <- root(tree, 'CMS_012_RNA_A_S4')
# 
# metadata_all <- read.csv("../data/project-mosquito_sample-table_cleaned.csv")
# metadata <- metadata_all[sapply(tree$tip.label, grep, metadata_all$sample_name), ]
# sp_all <- read.csv("../data/sample_genus_and_species.csv", header=FALSE, stringsAsFactors = FALSE)
# sp <- sp_all[sapply(tree$tip.label, grep, sp_all[, 1]), ]
# df <- data.frame(seq=tree$tip.label, 
#                  long=factor(round(metadata$long, 1)),
#                  lat=factor(round(metadata$lat, 1)),
#                  sp=factor(paste(sp[, 2], sp[, 3])))
# 
# tree_plot <- ggtree(tree) + theme_tree2() + xlim(0, .16)
# tree_plot <- tree_plot %<+% df + geom_tiplab(aes(color=sp)) + 
#   theme(legend.position="top", axis.line.x=element_line(size=1), axis.text.x=element_text(size=14),
#         legend.text=element_text(size=14), legend.title=element_text(size=14))
# tree_plot
# ggsave("../figures/CMS001_culex_COI_tree.pdf", tree_plot, width=10, height=8)
# 

tree <- read.tree("../../trees/culex_outgroup_COI_raxml.bestTree")
tree$tip.label <- gsub("_consensus", "", tree$tip.label)
tree <- root(tree, grep("Aedes", tree$tip.label, value=TRUE)) %>%
  drop.tip(., grep("Aedes", tree$tip.label, value=TRUE))

metadata_all <- read.csv("../../data/project-mosquito_sample-table_cleaned.csv")
metadata <- metadata_all[sapply(tree$tip.label, grep, metadata_all$sample_name), ]
sp_all <- read.csv("../../data/sample_genus_and_species.csv", header=FALSE, stringsAsFactors = FALSE)
sp <- sp_all[sapply(tree$tip.label, grep, sp_all[, 1]), ]



df <- data.frame(seq=tree$tip.label, 
                 long=factor(round(metadata$long, 1)),
                 lat=factor(round(metadata$lat, 1)),
                 Species=gsub("  ", " ", factor(paste(sp[, 2], sp[, 3]))),
                 location=c("NorCal", "SoCal")[(as.numeric(as.character(df$lat))>36)+1])
tree_plot <- ggtree(tree) + theme_tree2() 
tree_plot <- tree_plot %<+% df + geom_tiplab(aes(color=Species), align=TRUE) 
df2 <- select(df, location)
rownames(df2) <- df$seq
tree_plot <- gheatmap(tree_plot, df2, offset = .1, width=0.15, colnames_position="top", colnames_offset_y = 1) +
  theme(legend.position="top", axis.line.x=element_line(size=1), axis.text.x=element_text(size=12),
      legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  scale_fill_manual("Location", breaks=c("NorCal", "SoCal"), values=c("#3e50d5", "#578c00")) + 
  guides(color = guide_legend(order = 1, title.position="top", nrow=2), fill=guide_legend(order = 0, title.position="top"))
tree_plot


ggsave("../../figures/culex_outgroup_COI_tree.pdf", tree_plot, width=7, height=14)

