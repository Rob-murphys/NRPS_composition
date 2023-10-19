library(dplyr)
library(vegan)
library(tibble)

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
# Sort metadata rows to match feature table
metadata.sorted = metadata[match(rownames(feature_table_AD.norm), metadata$sample.id),]

## Removing Microtermes ##
metadata.sorted_noMicro = metadata.sorted %>%
  filter(genus != "Microtermes")

feature_table_16S.norm_noMicro = feature_table_16S.norm %>%
  filter(rownames(.) %in% metadata.sorted_noMicro$sample.id)

feature_table_AD.norm_noMicro = feature_table_AD.norm %>%
  filter(rownames(.) %in% metadata.sorted_noMicro$sample.id)

### Bray curtis dissimilarity ###
#===============================#

## 16S ##
bray_dis_16S.norm_noMicro = vegdist(feature_table_16S.norm_noMicro, method = "bray")
bray_dismat_16S.norm_noMicro = as.matrix(bray_dis_16S.norm_noMicro) # as a matrix

## AD ##
bray_dis_AD.norm_noMicro = vegdist(feature_table_AD.norm_noMicro, method = "bray")
bray_dismat_AD.norm_noMicro = as.matrix(bray_dis_AD.norm_noMicro) # as a matrix

### nMDS ordination ###
#=====================#
nmds_16S = vegan::metaMDS(bray_dis_16S.norm_noMicro, k = 4, maxit = 999, trymax = 500)

nmds_AD = vegan::metaMDS(bray_dis_AD.norm_noMicro, k = 4, maxit = 999, trymax = 500)

### Procrustes analysis ###
proc.res = procrustes(nmds_16S, nmds_AD, symmetric = TRUE)
prot.res = protest(nmds_16S, nmds_AD, permutations = 999)

residual_df = data.frame(res = residuals(proc.res)) %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata.sorted_noMicro, by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn")

# Testing if sample type and genus effect the residual error
residual.aov = aov(res ~ caste_phase * species, data = residual_df)
summary(residual.aov)
TukeyHSD(residual.aov,which = "species")
