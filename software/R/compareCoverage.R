library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = F)

mutation_samples <- read_delim("software/data/mutation_data/mutation_samples.txt", delim = "\t", skip = 1)
mutation_samples <- mutation_samples[rowSums(is.na(mutation_samples)) < 5,]
mutation_samples %<>% tbl_df

metabolite_samples <- read_delim("software/data/metabolomics_data/metabolite_samples.txt", delim = "\t") %>% tbl_df
metabolite_samples %<>% select(Sample_ID = Patient, Status) %>% mutate(Sample_ID_Status = paste0(Sample_ID, Status)) %>% unique

mutation_samples$`Alternate ID`[mutation_samples$`Alternate ID` != "-"] %>% length # 29 mutation samples
length(unique(metabolite_samples$Sample_ID)) # 49 mutation samples
print(length(intersect(mutation_samples$`Alternate ID`, metabolite_samples$Sample_ID))) # 22 samples in overlap
shared_IDs <- intersect(mutation_samples$`Alternate ID`, metabolite_samples$Sample_ID)
###

mutation_mutationlist <- read_delim("software/data/mutation_data/mutation_mutationlist.txt", delim = "\t", skip = 1)
mutation_mutationlist <- mutation_mutationlist[,colSums(is.na(mutation_mutationlist)) < 1000]
mutation_mutationlist %<>% tbl_df

mutation_mutationlist %>% count(`Gene Symbol`) %>% arrange(desc(n))


### 

metabolite_metaboliteQuant <- read_delim("software/data/metabolomics_data/metabolite_relative_abundance.txt", delim = "\t") %>%
  gather("Sample_ID_Status", "log2RA", -Metabolite) %>%
  mutate(Sample_ID = substr(Sample_ID_Status, 1, 10),
         Status = substr(Sample_ID_Status, 11, 100))

metabolite_RA <- metabolite_metaboliteQuant %>%
  group_by(Sample_ID, Metabolite) %>%
  summarize(TN_diff = log2RA[Status == "Tumor"] - log2RA[Status == "Normal"]) %>%
  filter(Sample_ID %in% shared_IDs) %>% # filter for shared IDs
  group_by(Metabolite) %>%
  filter(sum(is.na(TN_diff)) < n()/2) # remove poorly measured metabolites

# instead use PCs

#metabolite_RA <- metabolite_metaboliteQuant %>%
#  group_by(Sample_ID, Metabolite) %>%
#  summarize(TN_diff = log2RA[Status == "Tumor"] - log2RA[Status == "Normal"]) %>%
#  group_by(Metabolite) %>%
#  filter(sum(is.na(TN_diff)) < n()/2) %>%
#  spread(Sample_ID, TN_diff) 

#rownames(metabolite_RA) <- metabolite_RA$Metabolite
#metabolite_RA %<>% ungroup %>% select(-Metabolite) %>% as.matrix %>% apply(c(1,2), as.numeric)

#cluster::clara(metabolite_RA, k = 10)$medoids



# load previously identified metabolite which differ b/w tumor and benign

TNdiscoveries <- read.delim("software/data/metabolomics_data/TNdiscoveriesRenamed.tsv")
TNdiscovery_names <- read.delim("software/data/metabolomics_data/disc_colors.txt")
TNdiscoveries %<>% select(Discoveries) %>%
  left_join(TNdiscovery_names %>% select(Metabolite = Discoveries, Discoveries = Palatable.Name), by = "Discoveries")

# only look at metabolites that are up or down overall in tumors (rename to cleaner names from these discoveries)

metabolite_RA %<>% left_join(TNdiscoveries, by = "Metabolite") %>%
  filter(!is.na(Discoveries)) %>%
  ungroup %>%
  select(Sample_ID, Metabolite = Discoveries, TN_diff)

# recenter TN differences ?

#metabolite_RA %<>% group_by(Metabolite) %>%
#  mutate(TN_diff = TN_diff - mean(TN_diff, na.rm = T)) %>% ungroup


# Convert to heatmap form (cluster and rearrange by dendrogram)

metabolite_HM <- metabolite_RA %>% filter(!is.na(TN_diff)) %>%
  spread(key = Sample_ID, value = TN_diff)
rownames(metabolite_HM) <- metabolite_HM$Metabolite; metabolite_HM %<>% ungroup %>% select(-Metabolite)

#metabolite_dissimilarity <- 1-cor(t(as.matrix(metabolite_HM)), use="pairwise.complete.obs")
#cluster_features <- hclust(as.dist(metabolite_dissimilarity)) # clustering features
cluster_features <- hclust(dist(as.matrix(metabolite_HM)), method = "ward.D2") # clustering features
cluster_samples <- hclust(dist(t(as.matrix(metabolite_HM))), method = "ward.D2") # clustering samples
  
clustered_features <- data.frame(Metabolite = cluster_features$label[cluster_features$order], stringsAsFactors = F)  
clustered_samples <- data.frame(Sample_ID = cluster_samples$label[cluster_samples$order])  

metabolite_RA %<>% ungroup %>% mutate(Metabolite = factor(Metabolite, levels = clustered_features$Metabolite),
                                      Sample_ID = factor(Sample_ID, levels = clustered_samples$Sample_ID)) %>%
  mutate(TN_diff = ifelse(TN_diff < -3, -3, TN_diff),
         TN_diff = ifelse(TN_diff > 3, 3, TN_diff))

hm_theme <- theme_minimal() + theme(axis.text = element_text(size = 20), panel.grid = element_blank(),
                                    axis.title = element_blank(), legend.key.size = unit(1, "inches"),
                                    legend.text = element_text(size = 25), legend.title = element_text(size = 25),
                                    axis.text.x = element_text(angle = 90, hjust = 1),
                                    legend.position = "none")

metHM <- ggplot(metabolite_RA, aes(x = Sample_ID, y = Metabolite, fill = TN_diff)) + geom_tile() +
  scale_fill_gradient2("log2 concentration\nin tumor relative\nto benign", low = "green2", mid = "black", high = "firebrick1",
                       midpoint = 0, guide = "colourbar", breaks = seq(-3, 3, by = 1), labels = c("< -3", "-2", "-1", "0", "1", "2", "> 3")) + hm_theme +
  theme(axis.text.x = element_blank())

# Add in mutations

mutation_cutoff = 4
genes_of_interest <- c("MLL3", "MLL2", "MLL", "ARID1A", "CDKN2A", "KRAS", "SMAD4", "TGFBR2", "TP53")

#mutation_mutationlist %>% count(`Gene Symbol`) %>% filter(`Gene Symbol` %in% genes_of_interest)

mutation_mutationlist_subset <- mutation_mutationlist %>%
  select(Sample_ID = `Alternate ID`, Gene = `Gene Symbol`, Frac = `Mutant Tags (%)`) %>%
  filter(Sample_ID %in% shared_IDs) %>% unique %>%
  filter(Gene %in% genes_of_interest) %>%
  #mutate(Gene = ifelse(Gene %in% c("MLL3", "MLL2", "MLL"), "MLL", Gene)) %>%
  #ungroup %>% mutate(Gene = factor(Gene, levels = c("TP53", "TGFBR2", "SMAD4", "KRAS", "CDKN2A", "ARID1A", "MLL"))) %>%
  ungroup %>% mutate(Gene = factor(Gene, levels = rev(genes_of_interest))) %>%
  #group_by(Gene) %>% filter(n() >= mutation_cutoff) %>%
  #ungroup %>% mutate(Gene = factor(Gene, levels = Gene %>% table %>% data.frame %>%
  #                                   arrange(Freq) %>% select(1) %>% unlist)) %>%
  mutate(mutated = "mutant")
  
# add on categories where mutation is missing

all_mutation_gene_pairs <- expand.grid(Sample_ID = clustered_samples %>% unlist,
                                       Gene = unique(mutation_mutationlist_subset$Gene), stringsAsFactors = F) %>%
  mutate(Gene = factor(Gene, levels = levels(mutation_mutationlist_subset$Gene))) %>%
  anti_join(mutation_mutationlist_subset, by = c("Sample_ID", "Gene")) %>%
  mutate(mutated = "WT") %>% mutate(Frac = NA)

all_mutation_gene_pairs <- rbind(mutation_mutationlist_subset, all_mutation_gene_pairs) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = clustered_samples$Sample_ID))

mutHM <- ggplot(all_mutation_gene_pairs, aes(x = Sample_ID, y = Gene, fill = Frac)) + geom_tile(color = "black") + hm_theme +
  theme(axis.text.x = element_blank())


### Make a big heatmap

theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank()
)

feature_dendro <- ggplot(ggdendro::segment(ggdendro::dendro_data(as.dendrogram(cluster_features)))) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_none + coord_flip()

sample_dendro <- ggplot(ggdendro::segment(ggdendro::dendro_data(as.dendrogram(cluster_samples)))) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_none

library(grid)

pdf("analysis/metabolite_mutation_overlap/combinedHM.pdf", height = 20, width = 20)
grid::grid.newpage()
print(metHM, vp=viewport(0.6, 0.7, x=0.4, y=0.35))
print(mutHM, vp=viewport(0.4, 0.2, x=0.48, y=0.8))
print(sample_dendro, vp=viewport(0.52, 0.1, x=0.45, y=0.95))
print(feature_dendro, vp=viewport(0.1, 0.7, x=0.8, y=0.35))
dev.off()


### Stratify by mutations

all_met_genes <- metabolite_RA %>% left_join(all_mutation_gene_pairs, by = "Sample_ID")

ggplot(all_met_genes, aes(x = Sample_ID, y = Metabolite, fill = TN_diff)) + facet_grid(Gene ~ mutated) + geom_tile() +
  scale_fill_gradient2("log2 concentration\nin tumor relative\nto benign", low = "green2", mid = "black", high = "firebrick1",
                       midpoint = 0, guide = "colourbar", breaks = seq(-3, 3, by = 1), labels = c("< -3", "-2", "-1", "0", "1", "2", "> 3")) + hm_theme
