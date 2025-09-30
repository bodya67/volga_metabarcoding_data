library("phyloseq")
library('dplyr')
library('ggplot2')
library('gridExtra')
library('ape')
library('vegan')
library('FSA')
library('plyr')
library('dunn.test')
library('ggsci')
library('EcolUtils')
library('cluster')
library(zoo)
library(ampvis2)
library("viridis") 
library(ggpmisc)
library(MASS)
library(ggcorrplot)
library(ggvenn)
library(VennDiagram)
library(ggpubr)
library(MuMIn)
library(MASS)

########################################################### Proces the otu table. Delete singletons/doubletons, remove taxons, fix taxonomy. 
###########################################################
###########################################################

# read all data
dir.create('finish_vis/',showWarnings = F, recursive = T)
mapping_file <- read.table('mapping.tsv', sep = '\t', dec = ',', quote = "", header = T)
mapping_with_HPP <- read.table('../mapping_DZ.txt', sep = '\t', dec = ',', quote = "", header = T)
mapping_file$HPP <- mapping_with_HPP$HPP
mapping_file <- mapping_file[-c(139,140),]
mapping_file$distance_from_source_of_the_river <- rev(mapping_file$distance_from_source_of_the_river)
### I will disable it for a while. It makes lm model work worse.
#mapping_file[,6:13] <- apply(mapping_file[,6:13],2, function(x){
#  ifelse(is.na(x),mean(x,na.rm = T),x)
#})

my_cols <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "red", "darkgreen", "lightskyblue", "deeppink", "gray76", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                                        "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")


### read otu table
full_table <- read.table('full_table.txt', sep = '\t', header = T,comment.char	= '')
rownames(full_table) <- full_table$Row.names
full_table <- subset(full_table, full_table$X139 == 0)
### Subset no domain ASVs
no_domain <- subset(full_table, full_table$Domain == 'k:Eukaryota' | full_table$Domain == '')
no_domain[,2:(ncol(no_domain) - 8)] <- apply(no_domain[,2:(ncol(no_domain) - 8)], 2, as.numeric)
no_domain <- subset(no_domain, apply(no_domain[,2:(ncol(no_domain) - 8)] , 1, sum) > 2)
### Subset asvs
library("seqinr")
fastafile <- read.fasta('ASVs.fa',as.string = TRUE, set.attributes = FALSE, forceDNAtolower = F)
no_domain_asvs <- fastafile[names(fastafile) %in% no_domain$Row.names]
no_domain_asvs_sums <- apply(no_domain[,2:(ncol(no_domain) - 8)] , 1, sum)
names(no_domain_asvs) <- paste0(no_domain$Row.names,';size=',no_domain_asvs_sums,';')
write.fasta(sequences = no_domain_asvs, names = names(no_domain_asvs), nbchar = 100000, file.out = "no_domain_asvs_for_usearch.fa")

#write.table(subset(full_table, full_table$Domain == ''), 'no_domain_asvs.tsv', col.names = T, sep = '\t', row.names = F)
full_table$Kingdom[full_table$Domain == 'k:Eukaryota'] <- 'k:Eukaryota'
full_table <- subset(full_table, full_table$Domain != 'k:Eukaryota' & full_table$Domain != '' & full_table$Kingdom != 'k:Bacteria')
full_table <- subset(full_table, full_table$Kingdom != "k:Eukaryota:nucl")
full_table <- full_table[c(colnames(full_table[,1:141]),'Kingdom', 'Domain', 'Phylum', 'Class','Order', 'Family', 'Genus', 'Species')]
#full_table <- subset(full_table, full_table$Domain != 'k:Eukaryota' & full_table$Domain != '' | full_table$Kingdom != '')
full_table[full_table == ''] <- NA
full_table <- data.frame(t(na.locf(t(full_table))))
full_table[,2:(ncol(full_table) - 8)] <- apply(full_table[,2:(ncol(full_table) - 8)], 2, as.numeric)
full_table <- subset(full_table, apply(full_table[,2:(ncol(full_table) - 8)] , 1, sum) > 2)
# Make fixed taxonomy
tax_levels <- c('Domain', 'Phylum', 'Class','Order', 'Family', 'Genus', 'Species')
patterns <- c('(k:.*)|(.*_X)|(.*_X_)|(.*_XX)|(.*_XX_)|(.*_XXX)|(.*_XXX_)|(.*_XXXX)|(.*_XXXX_)', 
              '(k:.*)|(d:.*)',
              '(k:.*)|(d:.*)|(p:.*)|(.*_X)|(.*_X_)|(.*_XX)|(.*_XX_)|(.*_XXX)|(.*_XXX_)|(.*_XXXX)|(.*_XXXX_)',
              '(k:.*)|(d:.*)|(p:.*)|(c:.*)|(.*_X)|(.*_X_)|(.*_XX)|(.*_XX_)|(.*_XXX)|(.*_XXX_)|(.*_XXXX)|(.*_XXXX_)', 
              '(k:.*)|(d:.*)|(p:.*)|(c:.*)|(o:.*)|(.*_X)|(.*_X_)|(.*_XX)|(.*_XX_)|(.*_XXX)|(.*_XXX_)|(.*_XXXX)|(.*_XXXX_))',
              '(k:.*)|(d:.*)|(p:.*)|(c:.*)|(o:.*)|(f:.*)|(.*_X)|(.*_X_)|(.*_XX)|(.*_XX_)|(.*_XXX)|(.*_XXX_)|(.*_XXXX)|(.*_XXXX_)',
              '(k:.*)|(d:.*)|(p:.*)|(c:.*)|(o:.*)|(f:.*)|(g:.*)|(.*_X)|(.*_X_.*)|(.*_XX)|(.*_XX_.*)|(.*_XXX)|(.*_XXX_.*)|(.*_XXXX)|(.*_XXXX_.*)|(.*_XXXXX)|(.*_XXXXX.*)|(.*_XXXXX )')

full_table_fixed <- full_table
count_var <- 1
for (fix_name in tax_levels) {
  full_table_fixed[,fix_name] <- gsub(patterns[count_var],'Unclassified',full_table_fixed[,fix_name])
  count_var <- count_var + 1
}
full_table_fixed$Genus <- gsub('UnclassifiedX','Unclassified',full_table_fixed$Genus)
full_table <- full_table_fixed
full_table <- subset(full_table, !(full_table$Phylum %in% c('p:Opisthokonta-Fungi','p:Opisthokonta-Metazoa','p:Streptophyta-Streptophyta_X')))
full_table <- subset(full_table, full_table$Domain != 'd:Cryptista:nucl')

### separate TSAR into subclades 
groups <- c("Alveolata", "Telonemia", "Stramenopiles", "Rhizaria")

for (g in groups) {
  # находим строки, где в Phylum есть подстрока g
  idx <- grepl(g, full_table$Phylum, fixed = TRUE)
  # присваиваем в domain новое значение
  full_table$Domain[idx] <- paste0("d:", g)
}

full_table$Domain <- gsub('d:TSAR','Unclassified', full_table$Domain)

### Change factors order 
mapping_file$Fraction_size <- factor(mapping_file$Fraction_size, levels = c('picofraction', 'microfraction'))

### write table with separated TSAR
write.table(full_table, 'full_table_separated_tsar.tsv', sep = '\t')


###########################################################
###########################################################
###########################################################


########################################################### get env functions 
###########################################################
###########################################################
library(EukFunc)
assign_check <- assign_func(full_table[,c(1,142:149)])
write.table(assign_check, 'assigned_eukfunc.tsv', sep = '\t')

assign_check$main_functional_class[is.na(assign_check$main_functional_class)] <- 'unknown'

pdf('finish_vis/euk_func_annot.pdf')
barplot(table(assign_check$main_functional_class))
dev.off()

subset <- merge(full_table, count, by = 'Row.names', all.y = T)

sum(apply(full_table[,2:141], 2, sum)) / sum(apply(subset[,2:141], 2, sum))


###########################################################
###########################################################
###########################################################


########################################################### analysis of env factors
###########################################################
###########################################################

### save summary 
capture.output(summary(mapping_file[,6:13]), file = "finish_vis/env_factors/factors_summary.txt")


### env factors against distance 
dir.create('finish_vis/env_factors',showWarnings = F, recursive = T)

mapping_with_HPP <- read.table('../mapping_DZ.txt', sep = '\t', dec = ',', quote = "", header = T)
mapping_with_HPP <- mapping_with_HPP[-c(139,140),]
mapping_with_HPP$site_after_HPP <- gsub('b_HPP','normal',mapping_with_HPP$site_after_HPP)
mapping_file_melted <- mapping_file
mapping_file_melted$site_after_HPP <- mapping_with_HPP$site_after_HPP
mapping_file_melted$site_after_HPP <- factor(mapping_file_melted$site_after_HPP, levels = c('normal', 'HPP'))

library(tidyr)
mapping_file_melted <- mapping_file_melted %>%
  gather(key = "variable", value = "value", c(6,8,9,10,11,12,13))

factor(mapping_file_melted$site_after_HPP, levels = c('normal', 'HPP'))
env_factors_dependency <- ggplot(mapping_file_melted, aes(x = mapping_file_melted$distance_from_source_of_the_river, y = value,label = Sample, color = Fast_Flowing_vs_Slow_Flowing, shape = site_after_HPP)) + 
  geom_point() + 
  theme_bw() +
  #geom_text() +
  theme(axis.text.x = element_text(vjust = 0.5, size = 15, angle = 45),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        legend.position = 'top',
        legend.text=element_text(size=13),
        legend.title = element_text(size=13)) + 
  labs(title = "", x ='distance_from_source_of_the_river') +
  scale_x_continuous(n.breaks = 10) + 
  facet_wrap(~ variable,ncol = 2, scales = "free_y")

ggsave(
  paste0('finish_vis/env_factors/all_factors.pdf'),
  plot = env_factors_dependency,
  device = NULL,
  path = NULL,
  width = 9.5,
  height = 7,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

ggsave(
  paste0('finish_vis/env_factors/all_factors.jpg'),
  plot = env_factors_dependency,
  device = NULL,
  path = NULL,
  width = 9.5,
  height = 7,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

### env factors against flowing/reach 

boxplot_flowing_reach <- ggplot(mapping_file_melted, aes(x = Fast_Flowing_vs_Slow_Flowing,y = value)) +
  geom_boxplot() +
  theme_bw() + 
  theme(plot.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 30),
        axis.text.y = element_text(size = 30),
        axis.title=element_text(size=30),
        axis.title.y = element_blank()) + 
  scale_fill_manual(values=my_cols) + 
  facet_wrap(~ variable, scales = "free_y")

ggsave(
  paste0('finish_vis/env_factors/boxplot_all_factors.pdf'),
  plot = boxplot_flowing_reach,
  device = 'pdf',
  path = NULL,
  scale = 1,
  width = 13,
  height = 13,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)


### Create corrplot 

# Вычисление корреляционной матрицы
corr_matrix <- cor(mapping_file[, 6:13], use = "complete.obs", method = "spearman")

# Вычисление матрицы p-значений
p_matrix <- cor_pmat(mapping_file[, 6:13], method = "spearman")

ggsave(
  'finish_vis/env_factors/spearman_cor_test_of_ind_var.pdf',
  plot = ggcorrplot(corr_matrix, p.mat = p_matrix, insig = "blank", lab = TRUE),
  width = 9.5,
  height = 7,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

ggsave(
  'finish_vis/env_factors/spearman_cor_test_of_ind_var.jpg',
  plot = ggcorrplot(corr_matrix, p.mat = p_matrix, insig = "blank", lab = TRUE),
  width = 9.5,
  height = 7,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)


###########################################################
###########################################################
###########################################################


########################################################### Write count tables
###########################################################
###########################################################

tax_levels <- c('Domain', 'Phylum', 'Class','Order','Family','Genus', 'Species')
dir.create('finish_vis/counts_and_percents/',showWarnings = F, recursive = T)

for (z in tax_levels) {
  aggregated_full_table <- aggregate(full_table[,2:(ncol(full_table) - 8)], by = list(full_table[,z]), FUN = sum)
  aggregated_full_table_with_procents <- as.data.frame(prop.table(as.matrix(aggregated_full_table[,-1]),2)*100)
  aggregated_full_table_with_procents[,z] <- aggregated_full_table[,1]
  write.table(aggregated_full_table ,paste0('finish_vis/counts_and_percents/counts_sums_per_each_',z,'.txt'), sep = "\t", quote = F, row.names = F)
  write.table(aggregated_full_table_with_procents ,paste0('finish_vis/counts_and_percents/percents_per_each_',z,'.txt'), sep = "\t", quote = F, row.names = F)

}

#################################  per each domain
for (domain_var in unique(full_table$Domain)) {
  dir.create(paste0('finish_vis/counts_and_percents/',gsub(':','_',domain_var)), showWarnings = F, recursive = T)
  domain_subset <- subset(full_table, full_table$Domain == domain_var)
  for (z in c('Phylum', 'Class','Order','Family','Genus', 'Species')) {
    aggregated_full_table_domain <- aggregate(domain_subset[,2:(ncol(domain_subset) - 8)], by = list(domain_subset[,z]), FUN = sum) 
    aggregated_full_table_domain_with_procents <- as.data.frame(prop.table(as.matrix(aggregated_full_table_domain[,-1]),2)*100)
    aggregated_full_table_domain_with_procents$taxa <- aggregated_full_table_domain[,1]
    # calculations of total percents
    total_percents <- as.data.frame(apply(aggregated_full_table_domain[,2:ncol(aggregated_full_table_domain)], 1, sum))
    total_percents <- total_percents[,1]/sum(total_percents[,1]) * 100
    total_percents <- cbind(aggregated_full_table_domain[,1], total_percents)
    # save data
    write.table(total_percents ,paste0('finish_vis/counts_and_percents/',gsub(':','_',domain_var),'/total_percents_for_each_',z,'.txt'), sep = "\t", quote = F, row.names = F)
    write.table(aggregated_full_table_domain ,paste0('finish_vis/counts_and_percents/',gsub(':','_',domain_var),'/counts_sums_per_each_',z,'.txt'), sep = "\t", quote = F, row.names = F)
    write.table(aggregated_full_table_domain_with_procents ,paste0('finish_vis/counts_and_percents/',gsub(':','_',domain_var),'/percents_per_each_',z,'.txt'), sep = "\t", quote = F, row.names = F)
  }
}

###########################################################
###########################################################
###########################################################


########################################################### prepare data for the phyloseq
###########################################################
###########################################################

count_table <- full_table[,2:(ncol(full_table) - 8)]
row.names(count_table) <- full_table[,1]
count_table <- count_table[,paste0('X',1:140)]
colnames(count_table) <- gsub('X', '', colnames(count_table))
count_table <- subset(count_table, select=-c(139:140))
  
tax_table <- full_table[,(ncol(full_table) - 7):(ncol(full_table))]
row.names(tax_table) <- full_table[,1]
tax_table <- as.matrix(tax_table)


################################# Import data as physeq artifacts
OTU = otu_table(count_table, taxa_are_rows = T)
TAX = tax_table(tax_table)
samples =sample_data(mapping_file)
row.names(samples) <- 1:nrow(samples)
phytree <- read_tree('ASVs_for_analysis.afa')
physeq = phyloseq(OTU, TAX,samples, phytree)

###########################################################
###########################################################
###########################################################

########################################################### alpha diversity calculation
###########################################################
###########################################################

pdf('finish_vis/rarecurve_plot.pdf')
rarecurve(t(otu_table(physeq)), step=50, cex=0.5)
dev.off()


ds <- amp_load(
  count_table,
  metadata = NULL,
  taxonomy = tax_table,
  fasta = NULL,
  tree = NULL,
  pruneSingletons = FALSE
)


pdf('finish_vis/octave_plot.pdf', width = 15, height = 15)
amp_octave(ds,
           #group_by = "SampleID",
           #tax_aggregate = "Genus",
           scales = "free_y",
           num_threads = 8
)
dev.off()

###########################################################
###########################################################
###########################################################


########################################################### Create alpha diversity and merge it with mapping file
###########################################################
###########################################################

alpha_div <- estimate_richness(physeq, split = TRUE, measures = c("Observed","ACE","Chao1", "Shannon", 'Fisher', 'Simpson'))
alpha_div_with_mapping <- cbind(alpha_div, mapping_file)
alpha_div_with_mapping$Sample <- as.character(alpha_div_with_mapping$Sample)


a_pico <- subset(alpha_div_with_mapping, alpha_div_with_mapping$Fraction_size == 'picofraction')
a_micro <- subset(alpha_div_with_mapping, alpha_div_with_mapping$Fraction_size == 'microfraction')


summary(lm(Shannon ~ distance_from_source_of_the_river, a_micro))
summary(lm(Shannon ~ distance_from_source_of_the_river, a_pico))


pdf('finish_vis/alpha_div_plot1.pdf')
plot_richness(physeq, measures=c("Observed","ACE","Chao1", "Shannon", 'Simpson')) + facet_grid(~ Fraction_size, scale="free_x", drop=TRUE) + scale_x_discrete(limits=mapping_file$Sample)
dev.off()

###########################################################
###########################################################
###########################################################

########################################################### alpha div indexes and normality test
###########################################################
########################################################### 

mapping_with_HPP <- read.table('../mapping_DZ.txt', sep = '\t', dec = ',', quote = "", header = T)
mapping_with_HPP <- mapping_with_HPP[-c(139,140),]



dir.create('finish_vis/alpha_div/normality_check/', showWarnings = F, recursive = T)
a_index <- c("Observed","ACE","Chao1", "Shannon", 'Simpson')
for (k in a_index) {
  alpha_div_cycle_plot <- plot_richness(physeq, measures=k,color = 'Fast_Flowing_vs_Slow_Flowing',shape = as.factor(gsub('b_HPP','normal',mapping_with_HPP$site_after_HPP))) +
    scale_x_discrete(limits=mapping_file$Sample) +
    geom_point(size = 2.5) +
    facet_wrap(~ Fraction_size, scales = "free_x", ncol = 1) + 
    scale_shape_manual(values=c(25, 19)) + 
    ylab(paste0(k, ' index value')) +
    theme(plot.title = element_blank(),
          legend.position="bottom",
          legend.text=element_text(size=30),
          legend.title = element_blank(),
          text = element_text(size=21),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90),
          strip.text.x = element_text(size = 30))
    

  ggsave(
    paste0('finish_vis/alpha_div/',k,'.pdf'),
    plot = alpha_div_cycle_plot,
    device = NULL,
    path = NULL,
    width = 15,
    height = 10,
    scale = 1,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  )
  
  ggsave(
    paste0('finish_vis/alpha_div/',k,'.png'),
    plot = alpha_div_cycle_plot,
    device = 'png',
    path = NULL,
    width = 15,
    height = 10,
    scale = 1,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  )
  
  pdf(paste0('finish_vis/alpha_div/normality_check/',k,'.pdf'))
  hist(alpha_div_with_mapping[,k])
  dev.off()
}



###########################################################
###########################################################
###########################################################

########################################################### alpha div indexes per each domain
###########################################################
###########################################################  

a_index <- c("Observed","ACE","Chao1", "Shannon", 'Simpson')
for (k in a_index) {
  for (d in unique(full_table_fixed$Domain)) {
    physeq_d <- subset_taxa(physeq, Domain == d)
    dir.create(paste0('finish_vis/alpha_div/by_domains/',gsub(':','_',d),'/') ,showWarnings = F, recursive = T)
    alpha_div_cycle_plot <- plot_richness(physeq_d, measures=k) +
      scale_x_discrete(limits=mapping_file$Sample) +
      facet_wrap(~ Fraction_size, scales = "free_x", ncol = 1) + 
      ggtitle(k) + theme(plot.title = element_text(hjust = 0.5))
    
    
    ggsave(
      paste0('finish_vis/alpha_div/by_domains/',gsub(':','_',d),'/',k,'.pdf'),
      plot = alpha_div_cycle_plot,
      device = NULL,
      path = NULL,
      width = 9.5,
      height = 5,
      scale = 1,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
  }
}

###########################################################
###########################################################
###########################################################

########################################################### boxplots of HPP
###########################################################
###########################################################

mapping_with_HPP <- read.table('../mapping_DZ.txt', sep = '\t', dec = ',', quote = "", header = T)
mapping_with_HPP <- mapping_with_HPP[-c(139,140),]
alpha_div_with_mapping$HPP <- mapping_with_HPP$HPP
alpha_div_with_mapping$site_after_HPP = mapping_with_HPP$site_after_HPP

a_index <- c('Shannon', 'ACE', 'Simpson')
for (i in a_index) {
  boxplot_HPP<- ggplot(alpha_div_with_mapping[alpha_div_with_mapping$site_after_HPP != 'normal' & 
                                                alpha_div_with_mapping$Fraction_size == 'microfraction',], aes(x = site_after_HPP,y = ACE)) +
    geom_boxplot() +
    stat_compare_means(aes(group = site_after_HPP), label = "p.format", method = 'kruskal.test', size = 6) +
    labs(y = paste0('ACE', ' index'), x = '') + 
    theme_bw() + 
    theme(plot.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title=element_text(size=20)) + 
    scale_fill_manual(values=my_cols) 
  
  ggsave(
    paste0('finish_vis/alpha_div/',i,'_HPP_boxplots.pdf'),
    plot = boxplot_HPP,
    device = 'pdf',
    path = NULL,
    scale = 1,
    width = 5,
    height = 5,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  )
  ggsave(
    paste0('finish_vis/alpha_div/',i,'_HPP_boxplots.png'),
    plot = boxplot_HPP,
    device = 'png',
    path = NULL,
    scale = 1,
    width = 5,
    height = 5,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  ) 
}


### For each HPP
boxplot_HPP<- ggplot(alpha_div_with_mapping[alpha_div_with_mapping$site_after_HPP != 'normal' & 
                                              alpha_div_with_mapping$Fraction_size == 'microfraction',], aes(x = HPP,y = Shannon)) +
  geom_boxplot() +
  stat_compare_means(aes(group = site_after_HPP), label = "p.format", method = 'kruskal.test', size = 6) +
  labs(y = 'Shannon index', x = '') + 
  theme_bw() + 
  theme(plot.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20)) + 
  scale_fill_manual(values=my_cols) 

ggsave(
  'finish_vis/alpha_div/each_HPP_boxplots.pdf',
  plot = boxplot_HPP,
  device = 'pdf',
  path = NULL,
  scale = 1,
  width = 5,
  height = 5,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)
ggsave(
  'finish_vis/alpha_div/each_HPP_boxplots.png',
  plot = boxplot_HPP,
  device = 'png',
  path = NULL,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

capture.output(summary(aov(Shannon ~ HPP,alpha_div_with_mapping)),
               file = 'finish_vis/alpha_div/HPP_boxplots.txt',quote = F)

capture.output(dunn.test(alpha_div_with_mapping$Shannon, alpha_div_with_mapping$HPP, method = "bonferroni", list = T), file =  'finish_vis/alpha_div/HPP_dunn_boxplots.txt')



###########################################################
###########################################################
###########################################################

########################################################### plot simple plots of linear regression
###########################################################
###########################################################


a_index <- c("Observed","ACE","Chao1", "Shannon", 'Simpson')
for (k in a_index) {
  dir.create(paste0('finish_vis/a_div_dependency/separated_filters/',k) ,showWarnings = F, recursive = T)
  dir.create(paste0('finish_vis/a_div_dependency/both_filters/',k) ,showWarnings = F, recursive = T)
  dir.create(paste0('finish_vis/a_div_dependency/regression_asumptions/separated_filters/',k) ,showWarnings = F, recursive = T)
  dir.create(paste0('finish_vis/a_div_dependency/regression_asumptions/both_filters/',k) ,showWarnings = F, recursive = T)
  dir.create('finish_vis/a_div_dependency/separated_filters/multiple_lm/' ,showWarnings = F, recursive = T)
  dir.create('finish_vis/a_div_dependency/both_filters/multiple_lm/' ,showWarnings = F, recursive = T)
  dir.create(paste0('finish_vis/a_div_dependency/both_filters_separated_filters/',k) ,showWarnings = F, recursive = T)
  dir.create(paste0('finish_vis/a_div_dependency/both_filters_separated_filters/multiple_lm') ,showWarnings = F, recursive = T)
  
  print(paste0('dirs created;', 'k=',k))
  for (i in colnames(mapping_file[,6:13])) {
    a_div_dependency <- ggplot(alpha_div_with_mapping, aes(x = alpha_div_with_mapping[,i], y = alpha_div_with_mapping[,k],label = Sample)) + 
      geom_point() + 
      #geom_text() +
      geom_smooth(method='lm', formula= y~x) +
      theme(axis.text.x = element_text(vjust = 0.5, size = 15),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 20)) + 
      labs(title = "", y = paste0(k, ' index value'), x ='distance_from_source_of_the_river') 
      #facet_wrap(~ Filter_size, scales = "free_x", ncol = 1)

    ggsave(
      paste0('finish_vis/a_div_dependency/both_filters/',k,'/',i,'.pdf'),
      plot = a_div_dependency,
      device = NULL,
      path = NULL,
      width = 9.5,
      height = 3.5,
      scale = 1,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
    
    a_div_dependency_separated <- ggplot(alpha_div_with_mapping, aes(x = alpha_div_with_mapping[,i], y = alpha_div_with_mapping[,k],label = Sample)) + 
      geom_point() +   
      #geom_text() +
      geom_smooth(method='lm', formula= y~x) +
      labs(title = "", y = paste0(k, ' index value'), x = i) +
      facet_wrap(~ Fraction_size, scales = "free_x", ncol = 1) + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            strip.text.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 15)) 
    
    ggsave(
      paste0('finish_vis/a_div_dependency/separated_filters/',k,'/',i,'.pdf'),
      plot = a_div_dependency_separated,
      device = NULL,
      path = NULL,
      width = 9.5,
      height = 9.5,
      scale = 1,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
    
    merged_regression <- grid.arrange(a_div_dependency_separated,a_div_dependency, ncol = 1,heights=c(2,1.2))
    ggsave(
      paste0('finish_vis/a_div_dependency/both_filters_separated_filters/',k,'/',i,'.pdf'),
      plot = merged_regression,
      device = NULL,
      path = NULL,
      width = 9.5,
      height = 10.5,
      scale = 1,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
  }
    
    print(paste0('dependensies plotted;', 'k=',k,' i=',i))
   
    main_formula <- paste(k, "~", paste(c('Fraction_size','Fast_Flowing_vs_Slow_Flowing','Max_depth','distance_from_source_of_the_river','T','DO','TDS','ntu','Transparancy','Color'), collapse = " + "))
    main_formula_for_fractions <- paste(k, "~", paste(c('Fast_Flowing_vs_Slow_Flowing','Max_depth','distance_from_source_of_the_river','T','DO','TDS','ntu','Transparancy','Color'), collapse = " + "))
    fit <- lm(main_formula, data = na.omit(alpha_div_with_mapping))
    step <- stepAIC(fit, direction="both")
    step$anova # display results 
    best_factors <- colnames(step$model)
    formula_string <- paste(best_factors[1], "~", paste(best_factors[-1], collapse = " + "))
    capture.output(summary(lm(formula_string, na.omit(alpha_div_with_mapping))),
                   file = paste0('finish_vis/a_div_dependency/both_filters_separated_filters/multiple_lm/',k,'.txt'),quote = F)
    print('lm has been writen')
    #print(paste0('mlm for both fractions has been created;', 'k=',k,))
    
    filter_size <- c('picofraction', 'microfraction')
    for (filter_size_cycle in filter_size) {
      separated_df <- subset(alpha_div_with_mapping, alpha_div_with_mapping$Fraction_size == filter_size_cycle)
      fit <- lm(main_formula_for_fractions, data = na.omit(separated_df))
      step <- stepAIC(fit, direction="both")
      best_factors <- colnames(step$model)
      formula_string <- paste(best_factors[1], "~", paste(best_factors[-1], collapse = " + "))
      capture.output(summary(lm(separated_df[,k] ~ separated_df[,i], separated_df)),
                     file = paste0('finish_vis/a_div_dependency/separated_filters/',k,'/',i,'_',filter_size_cycle,'.txt'),quote = F)
      capture.output(summary(lm(formula_string, separated_df)),
                     file = paste0('finish_vis/a_div_dependency/separated_filters/multiple_lm/',k,'_',filter_size_cycle,'.txt'),quote = F)


  }
}

############################################################### plot model selection summary 

library(olsrr)
k <- ols_step_all_possible(fit)
k$result
write.csv(k$result, "results.csv", row.names = FALSE) # Save the data frame as csv
plot(k)

############################################################### 

############################################################### linear regression per each domain both fractions

for (domain_var in unique(full_table$Domain)) {
  if(domain_var=="Unclassified") next
  domain_subset <- subset_taxa(physeq, Domain==domain_var)
  alpha_div_by_domain <- estimate_richness(domain_subset, split = TRUE, measures = c("Observed","ACE","Chao1", "Shannon", 'Simpson'))
  alpha_div_with_mapping_by_domain <- cbind(alpha_div_by_domain, mapping_file)
  alpha_div_with_mapping_by_domain$Sample <- as.character(alpha_div_with_mapping_by_domain$Sample)
  
  for (k in a_index) {
    dir.create(paste0('finish_vis/a_div_dependency/both_filters_separated_filters/by_domains/',gsub(':','_',domain_var),'/',k) ,showWarnings = F, recursive = T)
    dir.create(paste0('finish_vis/a_div_dependency/both_filters_separated_filters/by_domains/',gsub(':','_',domain_var),'/multiple_lm') ,showWarnings = F, recursive = T)
    for (i in colnames(mapping_file[,6:13])) {
      a_div_dependency_separated <- ggplot(alpha_div_with_mapping_by_domain, aes(x = alpha_div_with_mapping_by_domain[,i], y = alpha_div_with_mapping_by_domain[,k],label = Sample)) + 
        geom_point() +   
        #geom_text() +
        geom_smooth(method='lm', formula= y~x) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        labs(title = "", y = k, x = NULL) +
        facet_wrap(~ Fraction_size, scales = "free_x", ncol = 1)
      
      a_div_dependency <- ggplot(alpha_div_with_mapping_by_domain, aes(x = alpha_div_with_mapping_by_domain[,i], y = alpha_div_with_mapping_by_domain[,k],label = Sample)) + 
        geom_point() +   
        #geom_text() +
        geom_smooth(method='lm', formula= y~x) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        labs(title = "", y = k, x = NULL)
      
      merged_regression <- grid.arrange(a_div_dependency,a_div_dependency_separated, ncol = 1,heights=c(1,2))
      
      ggsave(
        paste0('finish_vis/a_div_dependency/both_filters_separated_filters/by_domains/',gsub(':','_',domain_var),'/',k,'/',i,'.pdf'),
        plot = merged_regression,
        device = NULL,
        path = NULL,
        width = 9.5,
        height = 7,
        scale = 1,
        #units = c("in", "cm", "mm"),
        dpi = 300,
        limitsize = F
      )
      
      main_formula <- paste(k, "~", paste(c('Fraction_size','Fast_Flowing_vs_Slow_Flowing','Max_depth','distance_from_source_of_the_river','T','DO','TDS','ntu','Transparancy','Color'), collapse = " + "))
      fit <- lm(main_formula, data = na.omit(alpha_div_with_mapping_by_domain))
      step <- stepAIC(fit, direction="both")
      best_factors <- colnames(step$model)
      formula_string <- paste(best_factors[1], "~", paste(best_factors[-1], collapse = " + "))
      capture.output(summary(lm(formula_string, alpha_div_with_mapping_by_domain)),
                     file = paste0('finish_vis/a_div_dependency/both_filters_separated_filters/by_domains/',gsub(':','_',domain_var),'/multiple_lm/',k,'.txt'),quote = F)
      
      
    }
  }
  
}

###########################################################
###########################################################
###########################################################


########################################################### Correlation for all independent variables
###########################################################
########################################################### 

ggsave(
  'finish_vis/env_factors/spearman_cor_test_of_ind_var.pdf',
  plot = ggcorrplot(cor(subset(mapping_file[,6:13]), use = "complete.obs", method = 'spearman')),
  device = NULL,
  path = NULL,
  width = 9.5,
  height = 7,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)
library("PerformanceAnalytics")
pdf('finish_vis/env_factors/spearman_cor_test_of_ind_var_full.pdf')
chart.Correlation(cor(subset(mapping_file[,6:13]), use = "complete.obs", method = 'spearman'), histogram=TRUE, pch=19)
dev.off()


write.table(format(round(cor(subset(mapping_file[,6:13]), use = "complete.obs", method = 'spearman'), 2), nsmall = 2), 'finish_vis/env_factors/spearman_cor_test_of_ind_var.txt', sep = '\t')

###########################################################
###########################################################
###########################################################

########################################################### Boxplots of alpha diversity indexes
###########################################################
########################################################### 

# By domains
dir.create('finish_vis/alpha_div/boxplot/by_domains/',showWarnings = F, recursive = T)
alpha_div_by_domains <- data.frame()
for (i in unique(full_table$Domain)) {
  alpha_div_by_domains_inter_var <- estimate_richness(subset_taxa(physeq, Domain==i), split = TRUE, measures = c("Observed","ACE","Chao1", "Shannon", 'Simpson'))
  print('1')
  alpha_div_by_domains_inter_var$Domain <- i
  print('2')
  alpha_div_by_domains_inter_var <- cbind(alpha_div_by_domains_inter_var, mapping_file)
  print('3')
  alpha_div_by_domains <- rbind(alpha_div_by_domains,alpha_div_by_domains_inter_var)
  print('4')
  alpha_div_by_domains_inter_var <- NULL
  print(i)
  alpha_div_by_domains_plus_both_filters <- alpha_div_by_domains
  alpha_div_by_domains_plus_both_filters$Fraction_size <- 'picofraction+microfraction'
  alpha_div_by_domains <- rbind(alpha_div_by_domains,alpha_div_by_domains_plus_both_filters)
}
alpha_div_by_domains$Domain <- gsub('d:','', alpha_div_by_domains$Domain)

for (i in c("Observed","Chao1", "Shannon", 'Simpson')) {

  boxplot_by_domains <- ggplot(alpha_div_by_domains, aes(x = Domain,y = alpha_div_by_domains[,i], fill = Fast_Flowing_vs_Slow_Flowing)) +
    geom_boxplot() +
    stat_compare_means(aes(group = Fast_Flowing_vs_Slow_Flowing), label = "p.format", method = 'kruskal.test', size = 3) +
    theme(text = element_text(size=20)) + 
    labs(y = i) + 
    theme_bw() + 
    theme(legend.position="top", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    scale_fill_manual(values=my_cols) +
    facet_wrap(~Fraction_size, ncol = 1) + theme(legend.key.size = unit(2, 'cm'), #change legend key size
                                     legend.key.height = unit(2, 'cm'), #change legend key height
                                     legend.key.width = unit(2, 'cm'), #change legend key width
                                     legend.title = element_text(size=20), #change legend title font size
                                     legend.text = element_text(size=16))
  
  ggsave(
    paste0('finish_vis/alpha_div/boxplot/by_domains/',i,'.pdf'),
    plot = boxplot_by_domains,
    device = NULL,
    path = NULL,
    width = 15,
    height = 10,
    scale = 1,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  )
    
}



# by full data
boxplot_cycle = 'Shannon'
for (boxplot_cycle in a_index) {

    alpha_div_inter_var <- estimate_richness(subset_taxa(physeq), split = TRUE, measures = c("Observed","ACE","Chao1", "Shannon", 'Simpson'))
    print('1')
    print('2')
    alpha_div_inter_var <- cbind(alpha_div_inter_var, mapping_file)
    print('3')
    alpha_div <- alpha_div_inter_var
    print('4')
    print(i)
    alpha_div_plus_both_filters <- alpha_div
    alpha_div_plus_both_filters$Fraction_size <- 'picofraction+\nmicrofraction'
    alpha_div <- rbind(alpha_div,alpha_div_plus_both_filters)
  
  
  
  boxplot_filter_size <- ggboxplot(alpha_div, x = "Fraction_size", y = boxplot_cycle,
              color = "Fast_Flowing_vs_Slow_Flowing", palette = my_cols,
              add = "jitter") + 
              stat_compare_means(aes(group = Fast_Flowing_vs_Slow_Flowing), label = "p.format", method = 'kruskal.test') +
              theme(text = element_text(size=20))+
              theme(axis.title.x = element_blank()) +
              ylab(paste0(boxplot_cycle,' index value'))
              
    
    
    ggsave(
      paste0('finish_vis/alpha_div/boxplot/',boxplot_cycle,'_Type_of_section.pdf'),
      plot = boxplot_filter_size,
      device = NULL,
      path = NULL,
      width = 9.5,
      height = 7,
      scale = 1,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
    
    
    ggsave(
      paste0('finish_vis/alpha_div/boxplot/',boxplot_cycle,'_Type_of_section.jpg'),
      plot = boxplot_filter_size,
      device = NULL,
      path = NULL,
      width = 9.5,
      height = 7,
      scale = 1,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
    

  capture.output(kruskal.test(alpha_div_with_mapping[,boxplot_cycle] ~ Fraction_size, alpha_div_with_mapping), file =  paste0('finish_vis/alpha_div/boxplot/',boxplot_cycle,'.txt'))
  alpha_div_with_mapping_0.2 <- subset(alpha_div_with_mapping,alpha_div_with_mapping$Fraction_size == 'picofraction')
  alpha_div_with_mapping_3 <- subset(alpha_div_with_mapping,alpha_div_with_mapping$Fraction_size == 'microfraction')
  capture.output(kruskal.test(alpha_div_with_mapping_0.2[,boxplot_cycle] ~ Fast_Flowing_vs_Slow_Flowing, alpha_div_with_mapping_0.2), 
                 file =  paste0('finish_vis/alpha_div/boxplot/',boxplot_cycle,'_type_of_section_pico.txt'))
  capture.output(kruskal.test(alpha_div_with_mapping_3[,boxplot_cycle] ~ Fast_Flowing_vs_Slow_Flowing, alpha_div_with_mapping_3), 
                 file =  paste0('finish_vis/alpha_div/boxplot/',boxplot_cycle,'_type_of_section_micro.txt'))
}

###########################################################
###########################################################
###########################################################


########################################################### Venn diag of shared ASVs and shared ASV distibution 
###########################################################
###########################################################

# Data prep
dir.create('finish_vis/check_filter', recursive = T)
ven_diag <- apply(count_table, 1, function(x){
  
  c(sum(x[1:69]), sum(x[70:130]))
  
})
ven_diag_TF <- t(ifelse(ven_diag > 0, 1,0))
check_filter <- data.frame(ifelse((ven_diag_TF[,1] + ven_diag_TF[,2]) == 2, 'Both', 'One'))
ven_diag_TF <- cbind(ven_diag_TF, full_table$Domain)
colnames(ven_diag_TF) <- c('3um','0.2um', 'Domain')
ven_diag_TF <- as.data.frame(ven_diag_TF)

for (z in levels(as.factor(ven_diag_TF$Domain))) {
  ven_diag_TF_cycle <- subset(ven_diag_TF, ven_diag_TF[,3] == z)
  subset_for_venn_3um <- rownames(subset(ven_diag_TF_cycle,ven_diag_TF_cycle[,1] == 1))
  subset_for_venn_0.2um <- rownames(subset(ven_diag_TF_cycle,ven_diag_TF_cycle[,2] == 1))
  
  
  
  
  x <- list(
    '3um' = subset_for_venn_3um, 
    '0.2um' = subset_for_venn_0.2um
  )
  venn_diag <- ggvenn(
    x, 
    fill_color = my_cols,
    fill_alpha = 0.75,
    stroke_size = 0.5, set_name_size = 4
  ) + ggtitle(z) + theme(plot.title = element_text(size = 12, hjust = 0.5),
                         text = element_text(size=20))
                         
                         
  
  ggsave(
    paste0('finish_vis/check_filter/',gsub(':','_',z),'.pdf'),
    plot = venn_diag,
    device = NULL,
    path = NULL,
    #width = 9.5,
    #height = 7,
    scale = 1,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  )
}

### ven diagramm for all domains 

subset_for_venn_3um <- rownames(subset(ven_diag_TF,ven_diag_TF[,1] == 1))
subset_for_venn_0.2um <- rownames(subset(ven_diag_TF,ven_diag_TF[,2] == 1))
x <- list(
  'Microfraction' = subset_for_venn_3um, 
  'Picofraction' = subset_for_venn_0.2um
)
venn_diag <- ggvenn(
  x, 
  fill_color = my_cols,
  fill_alpha = 0.75,
  stroke_size = 0.5, 
  set_name_size = 7.5,
  text_size = 7.5
) 

ggsave(
  paste0('finish_vis/check_filter/full_data.pdf'),
  plot = venn_diag,
  device = NULL,
  path = NULL,
  #width = 9.5,
  #height = 7,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

###

check_filter$sum <- apply(count_table, 1, sum)
check_filter <- check_filter[rev(order(check_filter$sum)),]
check_filter$value <- 1
check_filter$number <- 1:nrow(check_filter)
colnames(check_filter) <- c('filter', 'suma', 'value', 'number')


dist_filter <- ggplot(check_filter, aes(x = number, y = value, fill = as.factor(filter))) +
  geom_tile() + 
  geom_line(aes(x = number, y = log(suma)/7)) + 
  scale_fill_manual(values=my_cols)

ggsave(
  paste0('finish_vis/check_filter/dist.pdf'),
  plot = dist_filter,
  device = NULL,
  path = NULL,
  #width = 9.5,
  #height = 7,
  scale = 1,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

###########################################################
###########################################################
###########################################################

########################################################### Check shared ASVs between fractions.
###########################################################
###########################################################

### read the data 
rarefied_asvs <- data.frame(otu_table(rarefied_otu_table))

### remove low impact
get_contribution <- data.frame(prop.table(apply(rarefied_asvs[,1:138], 1, sum))) * 100
keep_asvs <- rownames(subset(get_contribution, get_contribution[,1] > 0.1))
count_table_imp <- rarefied_asvs[rownames(rarefied_asvs) %in% keep_asvs,]


sum_fraction_cols <- apply(count_table_imp, 1, function(x){
  
  c(sum(x[1:69]), sum(x[70:138]))
  
})

shared_df_tax <- merge(t(sum_fraction_cols), tax_table, by = 0)
colnames(shared_df_tax)[2:3] <- c('microfraction','picofraction')

write.table(shared_df_tax, 'finish_vis/abundant_shared_asvs.tsv', sep = '\t',row.names = F)

###########################################################
###########################################################
###########################################################


########################################################### Pie charts of taxons across all samples
###########################################################
###########################################################

dir.create('finish_vis/taxanomy/pie/' ,showWarnings = F, recursive = T)
tax_levels_pie_chart <- c('Domain')
pc_df_agr_fin <- data.frame()
pc_df <- data.frame(apply(full_table[,2:(ncol(full_table) - 8)], 1, sum))
pc_df <- cbind(full_table[,1],pc_df, full_table[,(ncol(full_table) - 7):ncol(full_table)])

pc_df <- data.frame(apply(count_table[,1:69], 1, sum))
pc_df <- cbind(full_table[,1],pc_df, full_table[,(ncol(full_table) - 7):ncol(full_table)])



for (pie_chart_count in tax_levels_pie_chart) {
  pc_df_agr <- aggregate(pc_df[,2], by = list(pc_df[,pie_chart_count]), FUN = sum)
  pc_df_agr$Taxa <- pie_chart_count
  colnames(pc_df_agr) <- c('taxa','value','rank') 
  pc_df_agr$value <- pc_df_agr$value/sum(pc_df_agr$value)*100
  format(pc_df_agr$value, nsmall = 1)
  
  
  
  pie_chart <- ggplot(pc_df_agr, aes(x="", y=value, fill=paste0(taxa, ' ', round(value,2), '%'))) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    #theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    theme(plot.title = element_text(hjust = 0, size = 25), legend.text=element_text(size=20), legend.title=element_text(size=20)) + 
    labs(title = 'Microfraction', fill = 'Taxa')  + 
    scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "red", "darkgreen", "lightskyblue", "deeppink", "gray76", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                                 "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey"))
  

    #geom_text(aes(label = round(pc_df_agr[,2], digits = 2)), color = "white", position = position_stack(vjust = 0.5))
  
  #pc_df_agr_fin <- rbind(pc_df_agr_fin, pc_df_agr) - Эта строчка для создания массива чтобы можно было делать facet_wrap

  ggsave(
    paste0('finish_vis/taxanomy/pie/microfraction_taxonomy_pie_chart_both_fractions',pie_chart_count,'.pdf'),
    plot = pie_chart,
    device = 'pdf',
    path = NULL,
    scale = 1,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  )
  
  
}


###########################################################
###########################################################
###########################################################


########################################################### rarefaction
###########################################################
########################################################### 

sample_sums(physeq)
hist(sample_sums(physeq))
rarefied_otu_table <- rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),
                                        rngseed = T, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


# phyloseq object with percents instead of counts


physeq_percents <- transform_sample_counts(physeq, function(x) x / sum(x))

###########################################################
###########################################################
###########################################################


########################################################### Plot and analyse beta-div
###########################################################
########################################################### 

dir.create(paste0('finish_vis/beta_div/') ,showWarnings = F, recursive = T)
dist_methods <- c('bray', 'jaccard')
for (dist in dist_methods) {
  dist_mx <- vegdist(t(as.data.frame(otu_table(rarefied_otu_table))), method= dist)
  ord <- metaMDS(dist_mx, k = 2)
  (fit <- envfit(ord, mapping_file[,6:13], perm = 999, na.rm = T))
  capture.output(fit, 
                 file = paste0('finish_vis/beta_div/vectors_',dist,'.txt'))
  
  new_fit <-as.data.frame(cbind(scores(fit, display = "vectors"), 'r^2' = fit$vectors$r, 'p.value' = fit$vectors$pvals))
  new_fit$Species <- rownames(new_fit)
  
  #new_fit$Species <- c("Максимальная глубина","Км. от устья","T","DO","Соленость","Мутность","Прозрачность","Цветность" )
  
  df <- data.frame(ord$points)
  colnames(df) <- c('NMDS1', 'NMDS2')
  df <- cbind(df, mapping_file)
  #df$Type_of_section <- replace(df$Type_of_section, df$Type_of_section == 'flowing', 'Проточный участок')
  #df$Type_of_section <- replace(df$Type_of_section, df$Type_of_section == 'reach', 'Плесовый участок')
  
  beta_div <- ggplot(df, aes(x = NMDS1, y = NMDS2))+
    geom_point(aes(color = Fraction_size)) + 
    geom_segment(data = new_fit,
                 aes(x = 0, xend = NMDS1*3, y = 0, yend = NMDS2*3),
                 arrow = arrow(length = unit(0.1, "cm")), colour = ifelse(new_fit$p.value < 0.05, 'red', 'black')) +
    geom_label(data = new_fit, aes(x = NMDS1*3, y = NMDS2*3, label = Species),
              size = 2, colour = "black", fontface = "bold") + 
    theme_bw() + 
    scale_color_manual(values=my_cols) + 
    annotate("text", x = Inf, y = Inf, label = paste0("stress = ", round(ord$stress, 2)), hjust = 1, vjust = 1)
    
  
  
  
  ggsave(
    paste0('finish_vis/beta_div/NMDS_vec_',dist,'.pdf'),
    plot = beta_div,
    device = NULL,
    path = NULL,
    scale = 1,
    #units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = F
  )
  
}

#### Perform permanova.
for (dist in dist_methods) {
dist_mx <- vegdist(t(as.data.frame(otu_table(rarefied_otu_table))), method= dist)
mapping_file_permanova <- mapping_file 
mapping_file_permanova[,6:13] <- apply(mapping_file_permanova[,6:13],2, function(x){
  ifelse(is.na(x),mean(x,na.rm = T),x)
})
capture.output(adonis2(dist_mx ~ Fraction_size + Fast_Flowing_vs_Slow_Flowing + Max_depth + distance_from_source_of_the_river + T + DO + TDS + ntu + Transparancy + Color + HPP, data = mapping_file_permanova, permutations = 999), 
               file = paste0('finish_vis/beta_div/adonis_test_',dist,'.txt'))
}


#### Perform permanova separatly for fractions 
for (dist in dist_methods) {
  dist_mx <- vegdist(t(as.data.frame(otu_table(rarefied_otu_table))), method= dist)
  mapping_file_permanova <- mapping_file 
  mapping_file_permanova[,6:13] <- apply(mapping_file_permanova[,6:13],2, function(x){
    ifelse(is.na(x),mean(x,na.rm = T),x)
  })
  capture.output(adonis2(dist_mx ~ Fraction_size + Fast_Flowing_vs_Slow_Flowing + Max_depth + distance_from_source_of_the_river + T + DO + TDS + ntu + Transparancy + Color + HPP, data = mapping_file_permanova, permutations = 999), 
                 file = paste0('finish_vis/beta_div/adonis_test_',dist,'.txt'))
}


for (dist in dist_methods) {
  dist_mx <- vegdist(t(as.data.frame(otu_table(rarefied_otu_table))[,1:69]), method= dist)
  mapping_file_permanova <- mapping_file[1:69,] 
  mapping_file_permanova[,6:13] <- apply(mapping_file_permanova[,6:13],2, function(x){
    ifelse(is.na(x),mean(x,na.rm = T),x)
  })
  capture.output(adonis2(dist_mx ~ Fast_Flowing_vs_Slow_Flowing + Max_depth + distance_from_source_of_the_river + T + DO + TDS + ntu + Transparancy + Color + HPP, data = mapping_file_permanova, permutations = 999, by = 'margin'), 
                 file = paste0('finish_vis/beta_div/adonis_test_',dist,'_micro.txt'))
}

for (dist in dist_methods) {
  dist_mx <- vegdist(t(as.data.frame(otu_table(rarefied_otu_table))[,70:138]), method= dist)
  mapping_file_permanova <- mapping_file[70:138,]
  mapping_file_permanova[,6:13] <- apply(mapping_file_permanova[,6:13],2, function(x){
    ifelse(is.na(x),mean(x,na.rm = T),x)
  })
  capture.output(adonis2(dist_mx ~ Fast_Flowing_vs_Slow_Flowing + Max_depth + distance_from_source_of_the_river + T + DO + TDS + ntu + Transparancy + Color + HPP, data = mapping_file_permanova, permutations = 999, by = 'margin'), 
                 file = paste0('finish_vis/beta_div/adonis_test_',dist,'_pico.txt'))
}

###########################################################
###########################################################
###########################################################

########################################################### Plot taxonomy by ggplot 
###########################################################
###########################################################

df_inter <- data.frame()
prep <- cbind(count_table, full_table[,143:149])


for (i in unique(full_table$Domain)) {
  print(paste0('on domain: ',i))
  for (sub_taxa in c('Phylum', 'Class','Order', 'Family', 'Genus', 'Species')) {
    print(paste0('on taxa level: ',sub_taxa))
    dir.create(paste0('finish_vis/taxanomy/by_domains/',gsub(':','_',i)) ,showWarnings = F, recursive = T)
    for_taxa <- subset(prep, prep$Domain == i)
    for_taxa <- subset(for_taxa, for_taxa[,sub_taxa] != 'Unclassified')
    if(nrow(for_taxa) == 0) next # skip iteration 
    for_taxa <- aggregate(for_taxa[,1:(ncol(for_taxa) - 7)], by = list(for_taxa[,sub_taxa]), FUN = sum)
    ### for both filters 
    for_taxa_both_fractions = as.data.frame(cbind(for_taxa[,1],as.matrix(for_taxa[,2:70]) + as.matrix(for_taxa[,71:139])))
    colnames(for_taxa_both_fractions) = c("taxa",paste0(1:69,'/',70:138))
    for_taxa_both_fractions[,2:70] = apply(for_taxa_both_fractions[,2:70], 2, as.numeric)
    ###
    for_taxa[,2:139] <- prop.table(as.matrix(for_taxa[,2:139]), margin = 2)*100
    for_taxa[for_taxa == 'NaN'] <- 0
    ### for both filters 
    for_taxa_both_fractions[,2:70] <- prop.table(as.matrix(for_taxa_both_fractions[,2:70]), margin = 2)*100
    for_taxa_both_fractions[for_taxa_both_fractions == 'NaN'] <- 0
    ###
    # Calculation global percent by sum all samples
    #for_taxa$auf <- prop.table(apply(for_taxa[,2:139], 1, sum))*100
    #for_taxa[,1] <- apply(for_taxa, 1, function(x){
    
    # ifelse(sum(as.numeric(x[140])) < 0.5, 'Other', x[1])
    
    #})
    
  
    ### calculation for normal data
    for_taxa[,1] <- apply(for_taxa, 1, function(x){
  
      ifelse(sum(as.numeric(x[2:139])) < 40, 'Other', x[1])
      
    })
    
    df_inter <- data.frame(Values=unlist(for_taxa[,2:139]))
    df_inter$Taxa <- rep_len(for_taxa$Group.1, length.out=nrow(for_taxa))
    ### 
    ### both filters
    for_taxa_both_fractions[,1] <- apply(for_taxa_both_fractions, 1, function(x){
      
      ifelse(sum(as.numeric(x[2:70])) < 80, 'Other', x[1])
      
    })
    
    df_inter_both_fractions <- data.frame(Values=unlist(for_taxa_both_fractions[,2:70]))
    df_inter_both_fractions$Taxa <- rep_len(for_taxa_both_fractions$taxa, length.out=nrow(for_taxa_both_fractions))
    ###
    
    ### data frame changing for normal data
    v <- mapping_file$Sample
    n <- nrow(for_taxa)
    ffillv <- function(i) rep(v[i], n)
    df_inter$Sample <- c(sapply(seq_len(length(v)), ffillv))
    df_inter$Sample <- as.character(df_inter$Sample)
    v <- mapping_file$River_part
    df_inter$River_part <- c(sapply(seq_len(length(v)), ffillv))
    df_inter$River_part <- as.character(df_inter$River_part)
    v <- mapping_file$Fraction_size
    df_inter$Fraction_size <- c(sapply(seq_len(length(v)), ffillv))
    df_inter$Fraction_size <- as.character(df_inter$Fraction_size)
    ### 
    ### data frame changing for both filters
    v <- paste0(1:69,'/',70:138)
    n <- nrow(for_taxa_both_fractions)
    ffillv <- function(i) rep(v[i], n)
    df_inter_both_fractions$Sample <- c(sapply(seq_len(length(v)), ffillv))
    df_inter_both_fractions$Sample <- as.character(df_inter_both_fractions$Sample)
    v <- mapping_file[1:69,'River_part']
    df_inter_both_fractions$River_part <- c(sapply(seq_len(length(v)), ffillv))
    df_inter_both_fractions$River_part <- as.character(df_inter_both_fractions$River_part)
    v <- mapping_file[1:69,'Fraction_size']
    df_inter_both_fractions$Fraction_size <- 'picofraction+\nmicrofraction'

    ### 
    
    df_inter = rbind(df_inter,df_inter_both_fractions)
    
    df_inter$River_part <- factor(df_inter$River_part, levels = rev(levels(factor(df_inter$River_part))))
    
    taxa_barplot_my <- ggplot(df_inter, aes(x=factor(df_inter$Sample, levels = unique(df_inter$Sample)), y=Values, fill = Taxa)) +
      geom_bar(stat="identity") + 
      scale_fill_manual(values = my_cols) + 
      facet_wrap(~factor(Fraction_size, levels = c("picofraction", "microfraction", "picofraction+microfraction")) + River_part, scales = 'free_x', ncol = 3) +
      theme(plot.title = element_text(hjust = 0.5) ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.text=element_text(size=20),
            legend.position="bottom") + 
      ggtitle(i) 
    
    ggsave(
      paste0('finish_vis/taxanomy/by_domains/',gsub(':','_',i),'/',sub_taxa,'.pdf'),
      plot = taxa_barplot_my,
      device = 'pdf',
      path = NULL,
      scale = 1,
      width = 20,
      height = 10,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
    ggsave(
      paste0('finish_vis/taxanomy/by_domains/',gsub(':','_',i),'/',sub_taxa,'.png'),
      plot = taxa_barplot_my,
      device = 'png',
      path = NULL,
      scale = 1,
      width = 20,
      height = 10,
      #units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = F
    )
 }
}

################ plot only domains 

for_taxa <- subset(prep)
for_taxa <- subset(for_taxa, for_taxa$Domain != 'Unclassified')
for_taxa <- aggregate(for_taxa[,1:(ncol(for_taxa) - 7)], by = list(for_taxa[,'Domain']), FUN = sum)
### for both filters 
for_taxa_both_fractions = as.data.frame(cbind(for_taxa[,1],as.matrix(for_taxa[,2:70]) + as.matrix(for_taxa[,71:139])))
colnames(for_taxa_both_fractions) = c("taxa",paste0(1:69,'/',70:138))
for_taxa_both_fractions[,2:70] = apply(for_taxa_both_fractions[,2:70], 2, as.numeric)
###
for_taxa[,2:139] <- prop.table(as.matrix(for_taxa[,2:139]), margin = 2)*100
for_taxa[for_taxa == 'NaN'] <- 0
### for both filters 
for_taxa_both_fractions[,2:70] <- prop.table(as.matrix(for_taxa_both_fractions[,2:70]), margin = 2)*100
for_taxa_both_fractions[for_taxa_both_fractions == 'NaN'] <- 0
###

df_inter <- data.frame(Values=unlist(for_taxa[,2:139]))
df_inter$Taxa <- rep_len(for_taxa$Group.1, length.out=nrow(for_taxa))
### 
### both filters
df_inter_both_fractions <- data.frame(Values=unlist(for_taxa_both_fractions[,2:70]))
df_inter_both_fractions$Taxa <- rep_len(for_taxa_both_fractions$taxa, length.out=nrow(for_taxa_both_fractions))
###

### data frame changing for sep filters
v <- mapping_file$Sample
n <- as.numeric(nrow(for_taxa))
ffillv <- function(i) rep(v[i], n)
df_inter$Sample <- c(sapply(seq_len(length(v)), ffillv))
df_inter$Sample <- as.character(df_inter$Sample)
v <- mapping_file$River_part
df_inter$River_part <- c(sapply(seq_len(length(v)), ffillv))
df_inter$River_part <- as.character(df_inter$River_part)
v <- mapping_file$Fraction_size
df_inter$Fraction_size <- c(sapply(seq_len(length(v)), ffillv))
df_inter$Fraction_size <- as.character(df_inter$Fraction_size)
### 
### data frame changing for both filters
v <- paste0(1:69,'/',70:138)
n <- nrow(for_taxa_both_fractions)
ffillv <- function(i) rep(v[i], n)
df_inter_both_fractions$Sample <- c(sapply(seq_len(length(v)), ffillv))
df_inter_both_fractions$Sample <- as.character(df_inter_both_fractions$Sample)
v <- mapping_file[1:69,'River_part']
df_inter_both_fractions$River_part <- c(sapply(seq_len(length(v)), ffillv))
df_inter_both_fractions$River_part <- as.character(df_inter_both_fractions$River_part)
v <- mapping_file[1:69,'Fraction_size']
df_inter_both_fractions$Fraction_size <- 'picofraction+microfraction'
### 
df_inter = rbind(df_inter,df_inter_both_fractions)
df_inter$River_part <- factor(df_inter$River_part, levels = rev(levels(factor(df_inter$River_part))))


taxa_barplot_my <- ggplot(df_inter, aes(x=factor(df_inter$Sample, levels = unique(df_inter$Sample)), y=Values, fill = Taxa)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = my_cols) + 
  facet_wrap(~factor(Fraction_size, levels = c("picofraction", "microfraction", "picofraction+microfraction")) + River_part, scales = 'free_x', ncol = 3) +
  theme(plot.title = element_text(hjust = 0.5) ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=30),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.text.x = element_text(size=30),
        text = element_text(size=25))

####################################### without dividing by river part
taxa_barplot_my <- ggplot(df_inter, aes(x=factor(df_inter$Sample, levels = unique(df_inter$Sample)), y=Values, fill = Taxa)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = my_cols) + 
  facet_wrap(~factor(Fraction_size, levels = c("picofraction", "microfraction", "picofraction+microfraction")), scales = 'free_x', ncol = 1) +
  theme(plot.title = element_text(hjust = 0.5) ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=40),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.text.x = element_text(size=40),
        text = element_text(size=30))
#######################################
df_inter$Fraction_size = gsub('picofraction+microfraction', 'picofraction+microfraction', df_inter$Fraction_size)
rev_df_inter <- df_inter[rev(1:nrow(df_inter)), ]
taxa_barplot_my <- ggplot(rev_df_inter, aes(y=factor(rev_df_inter$Sample, levels = unique(rev_df_inter$Sample)), x=Values, fill = Taxa)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = my_cols) + 
  facet_wrap(~factor(Fraction_size, levels = c("picofraction", "microfraction", "picofraction+microfraction")), scales = 'free_y', ncol = 3) +
  theme(plot.title = element_text(hjust = 0.5) ,axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=35),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.text.x = element_text(size=35),
        text = element_text(size=30))

#######################################
  
ggsave(
  paste0('finish_vis/taxanomy/my_domains_percents.pdf'),
  plot = taxa_barplot_my,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 30,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

ggsave(
  paste0('finish_vis/taxanomy/my_domains_percents.png'),
  plot = taxa_barplot_my,
  device = 'png',
  path = NULL,
  width = 20,
  height = 30,
  #units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = F
)

###########################################################
###########################################################
###########################################################




