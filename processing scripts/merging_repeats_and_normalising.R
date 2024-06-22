library(qiime2R)
library(vegan)
library(dplyr)
library(tibble)
library(stringr)
library(compositions)

data_handle = function(infile, outfile1, outfile2){
  # Import data
  feature_table = as.data.frame(t(read_qza(infile)$data))
  
  # Merging technical replicates
  feature_table.merged = feature_table %>%
    rownames_to_column(var = "rn") %>%
    mutate(across("rn", str_replace, "-R1|-R2|-R3", "")) %>%
    dplyr::group_by(rn) %>%
    summarise(across(everything(), ~ mean(.))) %>%
    column_to_rownames(var = "rn")
  
  saveRDS(feature_table.merged, file = outfile1)
  # Data normalisation
  feature_table.merged.norm = feature_table.merged %>%
    sqrt(.) %>% # Squareoot standardization
    wisconsin(.) # Wisconson normalisation
  
  saveRDS(feature_table.merged.norm, file = outfile2)
}

clr_transform = function(infile, outfile1, outfile2, outfile3){
  
  # Import data
  feature_table = as.data.frame(t(read_qza(infile)$data))
  
  # Merging technical replicates
  feature_table.merged = feature_table %>%
    rownames_to_column(var = "rn") %>%
    mutate(across("rn", str_replace, "-R1|-R2|-R3", "")) %>%
    dplyr::group_by(rn) %>%
    summarise(across(everything(), ~ mean(.))) %>%
    column_to_rownames(var = "rn")
  
  saveRDS(feature_table.merged, file = outfile3)
  
  # Normalization
  merged_df = readRDS(merged_path)
  merged_df_clr = as.data.frame(clr(merged_df))
  saveRDS(merged_df_clr, file = outfile2)
  write.csv(file = outfile3, merged_df_clr, row.names = TRUE)
}

#### OLD #### 
#16S merged
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_mean_normalised_feature_table.rds")

# Comb 16S
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_comb_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_comb_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_comb_mean_normalised_feature_table.rds")

# Gut 16S
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_gut_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_gut_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_gut_mean_normalised_feature_table.rds")
#AD merged
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_normalised_feature_table.rds")

# Comb AD
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_comb_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_comb_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_comb_mean_normalised_feature_table.rds")

# Gut AD
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_gut_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_gut_mean_feature_table.rds",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_gut_mean_normalised_feature_table.rds")


metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  filter(!nest %in% c("negative-control", "positive-control", "environmental")) %>%
  mutate(across("sample.id", str_replace, "-R1|-R2|-R3", "")) %>%
  distinct(sample.id, .keep_all = TRUE)
# Remove NAs from metadata
metadata[metadata == ""] = NA

write.table(metadata, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", row.names = FALSE)


#====================#
## clr Transforming ##
#====================#

clr_transform = function(merged_path, outfile1, outfile2){
  merged_df = readRDS(merged_path)
  merged_df_clr = as.data.frame(clr(merged_df))
  saveRDS(merged_df_clr, file = outfile1)
  write.csv(file = outfile2, merged_df_clr, row.names = TRUE)
}

# 16S
clr_transform(merged_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_mean_feature_table.rds",
              outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_mean_clr_feature_table.rds",
              outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/16S_profiles_mean_clr_feature_table.csv")

# AD
clr_transform(merged_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_feature_table.rds",
              outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_clr_feature_table.rds",
              outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_clr_feature_table.csv")

#### ####


#### Redo ####

# 16S
## Normal ##
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_normalised_feature_table.rds")

## 99 ##
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_99_ASV_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_99_ASV_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_99_ASV_mean_normalised_feature_table.rds")

## 97 ##
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_97_ASV_merged_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_97_ASV_mean_feature_table.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_97_ASV_mean_normalised_feature_table.rds")

# AD
data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_merged_OBU_table-filtered.qza",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_OBU_table-filtered.rds",
            outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_normalised_OBU_table-filtered.rds")


#====================#
## clr Transforming ##
#====================#
# 16S
clr_transform(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_merged_table-filtered.qza",
              outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_mean_feature_table.rds",
              outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_mean_clr_feature_table.rds",
              outfile3 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_mean_clr_feature_table.csv")

# AD
clr_transform(merged_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_feature_table.rds",
              outfile2 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_clr_feature_table.rds",
              outfile3 = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/AD_profiles_paired_OBU_mean_clr_feature_table.csv")
