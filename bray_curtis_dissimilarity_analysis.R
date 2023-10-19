library(dplyr)
library(vegan)
library(pairwiseAdonis)
library(tibble)
library(microbiome)
library(phyloseq)
library(qiime2R)
library(MicEco)

#### Import Data ####
#===================#

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

## Phyloseq for genus level analysis ##
#=====================================#
taxa_table = read_qza("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/16S_profiles_redo_merged_seqs-taxonamy.qza")$data %>%
  parse_taxonomy() %>%
  filter(rownames(.) %in% colnames(feature_table_16S.norm)) %>%
  as.matrix() 

metadata_phylo = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  column_to_rownames(var = "sample.id")

## Make phyloseq object 
ASV_16S = otu_table(t(feature_table_16S.norm), taxa_are_rows = TRUE)
TAX_16S = tax_table(taxa_table)
samples_16S = sample_data(metadata_phylo)

physeq_16S = phyloseq(ASV_16S, TAX_16S, samples_16S)


physeq_genus = aggregate_taxa(physeq_16S, "Genus") %>% # collapse to genus rank
  subset_samples(., genus != "Microtermes") %>% # remove microtermes
  subset_samples(., species != "NA") # remove samples without species asignment

## Extract genus based feature table 
feature_table_genus = as.data.frame(t(otu_table(physeq_genus)))

#### Bray curtis dissimilarity ####
#=================================#

## 16S ##
bray_dis_16S.norm_noMicro = vegdist(feature_table_16S.norm_noMicro, method = "bray")
bray_dismat_16S.norm_noMicro = as.matrix(bray_dis_16S.norm_noMicro) # as a matrix

bray_dis_genus = vegdist(feature_table_genus, method = "bray")
bray_dismat_genus = as.matrix(bray_dis_genus) # as a matrix

## AD ##
bray_dis_AD.norm_noMicro = vegdist(feature_table_AD.norm_noMicro, method = "bray")
bray_dismat_AD.norm_noMicro = as.matrix(bray_dis_AD.norm_noMicro) # as a matrix


#### Permanovas ####
#==================#

#Take output of pairwise adonis2 and makes it into a dataframe becuase staring at huge lists is not fun
adonis2df = function(adonis_res){
  res_df = NULL
  for(pair in 2:length(adonis_res)){
    new_row = cbind(from = strsplit(names(adonis_res[pair]), "_vs_")[[1]][1],
                    to = strsplit(names(adonis_res[pair]), "_vs_")[[1]][2],
                    R2 = adonis_res[pair][[1]][3][1,1], 
                    fval = adonis_res[pair][[1]][4][1,1], 
                    pval = adonis_res[pair][[1]][5][1,1])
    res_df = rbind(res_df, new_row)
  }
  return(res_df)
}

filter_feature_table = function(ft_table, meta_df, samp_types){
  new_df = ft_table %>%
  rownames_to_column() %>%
  left_join(meta_df[colnames(meta_df) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 %in% samp_types) %>%
  column_to_rownames() %>%
  select(!c("caste_phase2"))
  return(new_df)
}

# Removing rows where species == NA
metadata.sorted_noMicro_noNA = metadata.sorted_noMicro %>%
  filter(!is.na(species))

# Change to only one soldier class - Merging all soldier castes
metadata.sorted_noMicro_noNA_oneS = metadata.sorted_noMicro_noNA %>%
  mutate(caste_phase2 = caste_phase) %>%
  mutate(caste_phase2 = replace(caste_phase2, caste_phase2 %in% c("S", "SS", "LS"), "S"))

  ## 16S ##
  #=======#

# Removing samples with NA species from feature table
feature_table_16S.norm_noMicro_noNA = feature_table_16S.norm_noMicro %>% 
  filter(rownames(.) %in% metadata.sorted_noMicro_noNA$sample.id)

adonis_16S.res = adonis2(feature_table_16S.norm_noMicro_noNA ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS, method = "bray", strata = metadata.sorted_noMicro_noNA_oneS$genus)
adonis_OmegaSq(adonis_16S.res)

    
  # Pairwise adonis #


adonis_pairwise_16S.res = pairwise.adonis2(feature_table_16S.norm_noMicro_noNA ~ caste_phase2 * species , data = metadata.sorted_noMicro_noNA_oneS, method = "bray", nperm = 9999, strata = metadata.sorted_noMicro_noNA_oneS$genus)

res1_df_16S = adonis2df(adonis_pairwise_16S.res)

# Effect size calculation for specific comparisons #
YC_SW_df = filter_feature_table(feature_table_16S.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "SW"))
YC_SW_res = adonis2(YC_SW_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "SW"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "SW"),]$genus)
adonis_OmegaSq(YC_SW_res)

YC_LW_df = filter_feature_table(feature_table_16S.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "LW"))
YC_LW_res = adonis2(YC_LW_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "LW"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "LW"),]$genus)
adonis_OmegaSq(YC_LW_res)

YC_S_df = filter_feature_table(feature_table_16S.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "S"))
YC_S_res = adonis2(YC_S_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "S"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "S"),]$genus)
adonis_OmegaSq(YC_S_res)

YC_MC_df = filter_feature_table(feature_table_16S.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "MC"))
YC_MC_res = adonis2(YC_MC_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "MC"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "MC"),]$genus)
adonis_OmegaSq(YC_MC_res)

YC_OC_df = filter_feature_table(feature_table_16S.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "OC"))
YC_OC_res = adonis2(YC_OC_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "OC"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "OC"),]$genus)
adonis_OmegaSq(YC_OC_res)

YC_C_df = filter_feature_table(feature_table_16S.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "C"))
YC_C_res = adonis2(YC_C_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "C"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "C"),]$genus)
adonis_OmegaSq(YC_C_res)

  ## AD ##
  #======#

# Removing samples with NA species from feature table
feature_table_AD.norm_noMicro_noNA = feature_table_AD.norm_noMicro %>%
  filter(rownames(.) %in% metadata.sorted_noMicro_noNA$sample.id)

adonis_AD.res = adonis2(feature_table_AD.norm_noMicro_noNA ~ caste_phase2 * species * genus, by = "terms", data = metadata.sorted_noMicro_noNA_oneS, method = "bray", strata = metadata.sorted_noMicro_noNA_oneS$genus)
adonis_OmegaSq(adonis_AD.res)

    # Pairwise adonis #


adonis_pairwise_AD.res = pairwise.adonis2(feature_table_AD.norm_noMicro_noNA ~ caste_phase2 * species , data = metadata.sorted_noMicro_noNA_oneS, nperm = 9999, strata = metadata.sorted_noMicro_noNA_oneS$genus)
res1_df_AD = adonis2df(adonis_pairwise_AD.res)

# Effect size calculation for specific comparisons #
YC_SW_df = filter_feature_table(feature_table_AD.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "SW"))
YC_SW_res = adonis2(YC_SW_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "SW"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "SW"),]$genus)
adonis_OmegaSq(YC_SW_res)

YC_LW_df = filter_feature_table(feature_table_AD.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "LW"))
YC_LW_res = adonis2(YC_LW_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "LW"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "LW"),]$genus)
adonis_OmegaSq(YC_LW_res)

YC_S_df = filter_feature_table(feature_table_AD.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "S"))
YC_S_res = adonis2(YC_S_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "S"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "S"),]$genus)
adonis_OmegaSq(YC_S_res)

YC_MC_df = filter_feature_table(feature_table_AD.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "MC"))
YC_MC_res = adonis2(YC_MC_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "MC"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "MC"),]$genus)
adonis_OmegaSq(YC_MC_res)

YC_OC_df = filter_feature_table(feature_table_AD.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "OC"))
YC_OC_res = adonis2(YC_OC_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "OC"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "OC"),]$genus)
adonis_OmegaSq(YC_OC_res)

YC_C_df = filter_feature_table(feature_table_AD.norm_noMicro_noNA, metadata.sorted_noMicro_noNA_oneS, c("YC", "C"))
YC_C_res = adonis2(YC_C_df ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "C"),], method = "bray", strata = metadata.sorted_noMicro_noNA_oneS[metadata.sorted_noMicro_noNA_oneS$caste_phase2 %in% c("YC", "C"),]$genus)
adonis_OmegaSq(YC_C_res)

  ## Genus ##
  #=========#

# Removing samples with NA species from feature table
feature_table_genus_noNA = feature_table_genus %>% 
  filter(rownames(.) %in% metadata.sorted_noMicro_noNA$sample.id)

adonis_16S.res = adonis2(feature_table_genus_noNA ~ caste_phase2* species, by = "terms", data = metadata.sorted_noMicro_noNA_oneS, method = "bray", strata = metadata.sorted_noMicro_noNA_oneS$genus)

# Pairwise adonis #

# Filtering distance matrix for pairwise adonis
bray_dismat_genus_noNA = bray_dismat_genus %>%
  as.data.frame() %>%
  filter(rownames(.) %in% metadata.sorted_noMicro_noNA$sample.id) %>%
  select(metadata.sorted_noMicro_noNA$sample.id) %>%
  as.matrix()

adonis_pairwise_16S.res = pairwise.adonis2(bray_dismat_genus_noNA ~ caste_phase2 * species , data = metadata.sorted_noMicro_noNA_oneS, nperm = 9999, strata = metadata.sorted_noMicro_noNA_oneS$genus)

res1_df_16S = adonis2df(adonis_pairwise_16S.res)