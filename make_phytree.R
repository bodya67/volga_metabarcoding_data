library(ggtree)
library(ggplot2)
library(phyloseq)
my_cols <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "red", "darkgreen", "lightskyblue", "deeppink", "gray76", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
             "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")


tree_contree <- read.tree("/media/bogdan/FF57-FFD2/projects/volga/18.04.2024/phytree/v3/iqtree2/v1/tree.contree")

tax_table <- read.table('/media/bogdan/FF57-FFD2/projects/volga/18.04.2024/phytree/v3/only_euk_taxonomy_classes_seqs_subset_pr2_asvs.csv', 
                        sep = ',',fill = T, col.names = c("id","Kingdom","Domain","Phylum","Class","Order","Family","Genus","Species"))
tax_table[tax_table == ''] <- 'Unknown' 
### fix names in tax table 
tax_table$id <- gsub('>','', tax_table$id)
### set new rownames 
rownames(tax_table) <- tax_table$id
count_table <- data.frame(cbind(tax_table[,0],1))
tax_table <- as.matrix(tax_table)


### import data as a phyloseq object 
OTU = otu_table(count_table, taxa_are_rows = T)
TAX = tax_table(tax_table)
phytree <- tree_contree
physeq = phyloseq(OTU,TAX,phytree)

tree_plot <- ggtree(physeq,aes(color = Domain),branch.length = 'none') + geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=2) +
  geom_tiplab(aes(label=paste0(Class,'_',id)), align=TRUE, linetype='dashed', linesize=.3, hjust=0, size = 1.8) +
  #geom_point(aes(x=x+hjust, color=SampleType, shape=Family, size=Abundance),na.rm=TRUE) +
  #scale_size_continuous(trans=log_trans(5)) +
  theme(legend.position="top") + #ggtitle("reproduce phyloseq by ggtree") + 
  scale_fill_manual(values=my_cols) + 
  #coord_cartesian(clip = 'off')+ 
  scale_color_manual(values=my_cols)# +
  #xlim(-0.5, 10)

ggsave(
  '/media/bogdan/FF57-FFD2/projects/volga/18.04.2024/finish_vis/phytree/v1/tree_contree.pdf',
  plot = tree_plot,
  device = 'pdf',
  path = NULL,
  #scale = 0.1,
  width = 15,
  height = 50,
  #units = "px",
  dpi = 300,
  limitsize = F
)

