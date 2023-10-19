library(dplyr)
library(vegan)
library(tibble)
library(ggplot2)
library(cowplot)
library(qiime2R)
library(phyloseq)
library(microbiome)


### Import Data ###
#=================#

## Already normalized data ##

# 16S
feature_table_16S.norm = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_normalised_feature_table.rds")
# AD
feature_table_AD.norm = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_normalised_OBU_table-filtered.rds")


## Meatdata ##

metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  mutate(sample_type = case_when(caste_phase %in% c("YC", "MC", "OC", "C") ~ "Comb",
                                 caste_phase %in% c("SW", "LW", "SS", "LS", "S") ~ "Gut"))


## Removing Microtermes ##
metadata_noMicro = metadata %>%
  filter(genus != "Microtermes")

feature_table_16S.norm_noMicro = feature_table_16S.norm %>%
  filter(rownames(.) %in% metadata_noMicro$sample.id)

feature_table_AD.norm_noMicro = feature_table_AD.norm %>%
  filter(rownames(.) %in% metadata_noMicro$sample.id)

# Phyloseq for genus level plot
taxa_table = read_qza("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/16S_profiles_redo_merged_seqs-taxonamy.qza")$data %>%
  parse_taxonomy() %>%
  filter(rownames(.) %in% colnames(feature_table_16S.norm)) %>%
  as.matrix() 

metadata_phylo = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  column_to_rownames(var = "sample.id")

## Make phyloseq object ##
ASV_16S = otu_table(t(feature_table_16S.norm), taxa_are_rows = TRUE)
TAX_16S = tax_table(taxa_table)
samples_16S = sample_data(metadata_phylo)

physeq_16S = phyloseq(ASV_16S, TAX_16S, samples_16S)


physeq_genus = aggregate_taxa(physeq_16S, "Genus") %>% # collapse to genus rank
  subset_samples(., genus != "Microtermes") %>% # remove microtermes
  subset_samples(., species != "NA") # remove samples without species asignment

## Extract genus based feature table ##
feature_table_genus = as.data.frame(t(otu_table(physeq_genus)))

### Bray curtis dissimilarity ###
#===============================#

## 16S ##
bray_dis_16S.norm_noMicro = vegdist(feature_table_16S.norm_noMicro, method = "bray")
bray_dismat_16S.norm_noMicro = as.matrix(bray_dis_16S.norm_noMicro) # as a matrix

bray_dis_16S_genus = vegdist(feature_table_genus, method = "bray")
bray_dismat_16S_genus = as.matrix(bray_dis_16S_genus) # as a matrix

## AD ##
bray_dis_AD.norm_noMicro = vegdist(feature_table_AD.norm_noMicro, method = "bray")
bray_dismat_AD.norm_noMicro = as.matrix(bray_dis_AD.norm_noMicro) # as a matrix

# factor levels
levels(metadata_noMicro$caste_phase) = c(levels(metadata_noMicro$caste_phase), "W")

metadata_noMicro_oneS_W = metadata_noMicro %>%
  mutate(caste_phase2 = caste_phase) %>%
  mutate(caste_phase2 = replace(caste_phase2, caste_phase2 %in% c("S", "SS", "LS"), "S"),
         caste_phase2 = replace(caste_phase2, caste_phase2 %in% c("SW", "LW"), "W"))

shapes = c(C = 4, YC = 17, MC = 3, OC = 18, W = 15 , S = 16)

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")

#### Pannel A ####
#================#

### AVS plot ###
nmds_16S = vegan::metaMDS(bray_dis_16S.norm_noMicro, k = 4, maxit = 999, trymax = 500) # NMDS of the normalized feature table in bray-curtis space
nmds_16S.scores = as.data.frame(nmds_16S$points)

nmds_16S.scores.meta = nmds_16S.scores %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata_noMicro_oneS_W, by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn")
nmds_16S.scores.meta$caste_phase = factor(nmds_16S.scores.meta$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))
nmds_16S.scores.meta$caste_phase2 = factor(nmds_16S.scores.meta$caste_phase2, levels = c("C", "YC", "MC", "OC", "W", "S"))


NMDS_1_2.species.genus.cent = aggregate(cbind(MDS1, MDS2) ~ species, data = nmds_16S.scores.meta, FUN = mean)

NMDS_1_2.species.genus.segs <- merge(nmds_16S.scores.meta, setNames(NMDS_1_2.species.genus.cent, c('species','oNMDS1','oNMDS2')),
                                     by = 'species', sort = FALSE)

NMDS_1_2.species.genus = ggplot(nmds_16S.scores.meta, aes(x = MDS1, y = MDS2))+
  geom_point(size = 4, aes(colour = species, shape = caste_phase2))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  scale_shape_manual(values = shapes)+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = NMDS_1_2.species.genus.segs, mapping = aes(xend = oNMDS1, yend = oNMDS2),
               alpha = 0.2)+
  geom_point(data = NMDS_1_2.species.genus.cent, size = 1)+
  geom_text(data = NMDS_1_2.species.genus.cent, aes(label = species), size = 12)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/pannel A.pdf", NMDS_1_2.species.genus, width = 12, height = 12)

#### Pannel A end #### 

#### Pannel B ####

nmds_16S_genus = vegan::metaMDS(bray_dis_16S_genus, k = 4, maxit = 999, trymax = 500) # NMDS of the normalized feature table in bray-curtis space
nmds_16S_genus.scores = as.data.frame(nmds_16S_genus$points)

nmds_16S_genus.scores.meta = nmds_16S_genus.scores %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata_noMicro_oneS_W, by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn")

nmds_16S_genus.scores.meta$caste_phase = factor(nmds_16S_genus.scores.meta$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))
nmds_16S_genus.scores.meta$caste_phase2 = factor(nmds_16S_genus.scores.meta$caste_phase2, levels = c("C", "YC", "MC", "OC", "W", "S"))

NMDS_1_2.caste_phase.phase.cent = aggregate(cbind(MDS1, MDS2) ~ caste_phase2, data = nmds_16S_genus.scores.meta, FUN = mean)

NMDS_1_2.caste_phase.phase.segs <- merge(nmds_16S_genus.scores.meta, setNames(NMDS_1_2.caste_phase.phase.cent, c('caste_phase2','oNMDS1','oNMDS2')),
                                         by = 'caste_phase2', sort = FALSE)

NMDS_1_2.caste_phase.phase = ggplot(nmds_16S_genus.scores.meta, aes(x = MDS1, y = MDS2))+
  geom_point(size = 4, aes(colour = species, shape = caste_phase2))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  scale_shape_manual(values = shapes)+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = NMDS_1_2.caste_phase.phase.segs, mapping = aes(xend = oNMDS1, yend = oNMDS2),
               alpha = 0.2)+
  geom_point(data = NMDS_1_2.caste_phase.phase.cent, size = 1)+
  geom_text(data = NMDS_1_2.caste_phase.phase.cent, aes(label = caste_phase2), size = 12)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/pannel B.pdf", NMDS_1_2.caste_phase.phase, width = 12, height = 12)
#### Pannel B end #### 

#### Pannel C ####
nmds_AD = vegan::metaMDS(bray_dis_AD.norm_noMicro, k = 4, maxit = 999, trymax = 500) # NMDS of the normalized feature table in bray-curtis space
nmds_AD.scores = as.data.frame(nmds_AD$points)

nmds_AD.scores.meta = nmds_AD.scores %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata_noMicro_oneS_W, by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn")
nmds_AD.scores.meta$caste_phase = factor(nmds_AD.scores.meta$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))
nmds_AD.scores.meta$caste_phase2 = factor(nmds_AD.scores.meta$caste_phase2, levels = c("C", "YC", "MC", "OC", "W", "S"))

NMDS_1_2_AD.species.genus.cent = aggregate(cbind(MDS1, MDS2) ~ species, data = nmds_AD.scores.meta, FUN = mean)

NMDS_1_2_AD.species.genus.segs <- merge(nmds_AD.scores.meta, setNames(NMDS_1_2_AD.species.genus.cent, c('species','oNMDS1','oNMDS2')),
                                     by = 'species', sort = FALSE)

NMDS_1_2_AD.species.genus = ggplot(nmds_AD.scores.meta, aes(x = MDS1, y = MDS2))+
  geom_point(size = 4, aes(colour = species, shape = caste_phase2))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  scale_shape_manual(values = shapes)+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = NMDS_1_2_AD.species.genus.segs, mapping = aes(xend = oNMDS1, yend = oNMDS2),
               alpha = 0.2)+
  geom_point(data = NMDS_1_2_AD.species.genus.cent, size = 1)+
  geom_text(data = NMDS_1_2_AD.species.genus.cent, aes(label = species), size = 12)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/pannel C.pdf", NMDS_1_2_AD.species.genus, width = 12, height = 12)
#### Pannel C end #### 

#### Pannel D ####
proc.res = procrustes(nmds_16S, nmds_AD, symmetric = TRUE)
prot.res = protest(nmds_16S, nmds_AD, permutations = 999)

ctest <- data.frame(rda1=proc.res$Yrot[,1],
                    rda2=proc.res$Yrot[,2],
                    xrda1=proc.res$X[,1],
                    xrda2=proc.res$X[,2],each=8) %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata_noMicro, by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn")

procrustes_plot = ggplot(ctest) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2, colour = sample_type), linewidth = 1) +
  geom_point(aes(x=rda1, y=rda2, colour = species, shape = genus), size = 4, shape = 16) +
  geom_point(aes(x=xrda1, y=xrda2, colour = species, shape = genus), size = 4, shape = 17) +
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

ggsave(filename = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/pannel D.pdf", procrustes_plot, width = 12, height = 12)

#### Pannel D end #### 