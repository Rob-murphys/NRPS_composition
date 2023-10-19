library(qiime2R)
library(tibble)
library(dplyr)
### Import Data ###
#=================#

## Un-normalized data ##
pre_feature_table_AD = as.data.frame(t(read_qza("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/pre_filt_feature_tables/AD_profiles_paired_redo_merged_OBU_table.qza")$data))
feature_table_AD = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_OBU_table-filtered.rds") 

pre_feature_table_16S = as.data.frame(t(read_qza("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/pre_filt_feature_tables/16S_profiles_redo_merged_table-dada2.qza")$data))
feature_table_16S = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_feature_table.rds") 

## Meatdata ##

metadata_AD = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_decontam_AD.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

metadata_16S = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_decontam_16S.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)



#### General stats about read counts etc ####



## Total OBU read count 

# Before
sum(pre_feature_table_AD) # 6763249 reads

df_AD = pre_feature_table_AD %>%
  rownames_to_column() %>%
  left_join(metadata_AD[colnames(metadata_AD) %in% c("sample.id", "Sample_or_Control")], by = c("rowname" = "sample.id")) %>%
  filter(!Sample_or_Control == "Control Sample") %>%
  select(!Sample_or_Control) %>%
  column_to_rownames()

rowsum_AD = rowSums(df_AD)
mean_AD = mean(rowsum_AD) # 15105.28
se_AD = sd(rowsum_AD)/sqrt(length(rowsum_AD)) # 691.2423

# After decontam, removing low abundance OBUs and averaging across triplicates
ncol(feature_table_AD) # 2407 OBUs
sum(feature_table_AD) # 1723428 reads

## Total ASV read count 

# Before
sum(pre_feature_table_16S) # 74445712 reads

df_16S = pre_feature_table_16S %>%
  rownames_to_column() %>%
  left_join(metadata_16S[colnames(metadata_16S) %in% c("sample.id", "Sample_or_Control")], by = c("rowname" = "sample.id")) %>%
  filter(!Sample_or_Control == "Control Sample") %>%
  select(!Sample_or_Control) %>%
  column_to_rownames()

rowsum_16S = rowSums(df_16S)
mean_16S = mean(rowsum_16S) # 159592
se_16S= sd(rowsum_16S)/sqrt(length(rowsum_16S)) # 3407.598

# After decontam, removing low abundance OBUs and averaging across triplicates
ncol(feature_table_16S) # 15960 ASVs
sum(feature_table_16S) # 19967957 reads



