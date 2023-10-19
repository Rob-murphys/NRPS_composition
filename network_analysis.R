library(tibble)
library(dplyr)
library(SpiecEasi)
library(igraph)
library(intergraph)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(GGally)
library(sna)
library(ggridges)
library(effectsize)

### Import Data ###
#=================#

# 16S
feature_table_16S = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_feature_table.rds")

# AD
feature_table_AD = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_OBU_table-filtered_2.rds") 

metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  mutate(sample_type = case_when(caste_phase %in% c("YC", "MC", "OC", "C") ~ "Comb",
                                 caste_phase %in% c("SW", "LW", "SS", "LS", "S") ~ "Gut"))
## Removing Microtermes ##
metadata_noMicro = metadata %>%
  filter(genus != "Microtermes")

feature_table_16S_noMicro = feature_table_16S %>%
  filter(rownames(.) %in% metadata_noMicro$sample.id)

feature_table_AD_noMicro = feature_table_AD %>%
  filter(rownames(.) %in% metadata_noMicro$sample.id)                                


#### Partitioning 16S data ####
#=============================#

  ## Guts ##

## Macrotermes ##
gut_16S_beli = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "bellicosus", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

gut_16S_sub = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "subhyalinus", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Ancistrotermes ##
gut_16S_cavi = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "cavithorax", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.60))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing columns that have colsums of 0 but in a stupid dplyr way

gut_16S_guin = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "guinensis", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.60))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing columns that have colsums of 0 but in a stupid dplyr way

## Odontotermes ##

gut_16S_odoA = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species A", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

gut_16S_odoB = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species B", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Pseudocanthotermes ##
gut_16S_mili = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "militaris-spiniger", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

  ## Comb ##

## Macrotermes ##
comb_16S_beli = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "bellicosus", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

comb_16S_sub = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "subhyalinus", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Ancistrotermes ##
comb_16S_cavi = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "cavithorax", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.85))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

comb_16S_guin = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "guinensis", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.85))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Odontotermes ##

comb_16S_odoA = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species A", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

comb_16S_odoB = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species B", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Pseudocanthotermes ##
comb_16S_mili = feature_table_16S_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "militaris-spiniger", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
  #select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

#### 16S data end ####

#### Partitioning AD data ####
#============================#

## Guts ##

## Macrotermes ##
gut_AD_beli = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "bellicosus", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

gut_AD_sub = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "subhyalinus", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Ancistrotermes ##
gut_AD_cavi = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "cavithorax", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing columns that have colsums of 0 but in a stupid dplyr way

gut_AD_guin = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "guinensis", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing columns that have colsums of 0 but in a stupid dplyr way

## Odontotermes ##

gut_AD_odoA = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species A", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

gut_AD_odoB = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species B", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Pseudocanthotermes ##
gut_AD_mili = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "militaris-spiniger", sample_type == "Gut") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Comb ##

## Macrotermes ##
comb_AD_beli = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "bellicosus", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

comb_AD_sub = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "subhyalinus", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Ancistrotermes ##
comb_AD_cavi = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "cavithorax", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.60))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

comb_AD_guin = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "guinensis", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.60))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Odontotermes ##

comb_AD_odoA = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species A", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

comb_AD_odoB = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species B", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

## Pseudocanthotermes ##
comb_AD_mili = feature_table_AD_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_noMicro[colnames(metadata_noMicro) %in% c("sample.id", "sample_type", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "militaris-spiniger", sample_type == "Comb") %>%
  column_to_rownames() %>%
  select(!c(sample_type, species)) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.50))
#select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) # Removing colums that have colsums of 0 but in a stupid dplyr way

#### AD data end ####

#### spiec.easi 16S networks with "mb" method ####
#================================================#

  ## Gut ##

# Macrotermes
mb_res_16S_belli_gut = spiec.easi(data = as.matrix(gut_16S_beli), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_belli_gut)
saveRDS(mb_res_16S_belli_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_belli_gut.rds")

mb_res_16S_sub_gut = spiec.easi(data = as.matrix(gut_16S_sub), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_sub_gut)
saveRDS(mb_res_16S_sub_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_sub_gut.rds")

#Ancistrotermes
mb_res_16S_cavi_gut = spiec.easi(data = as.matrix(gut_16S_cavi), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_cavi_gut)
saveRDS(mb_res_16S_cavi_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_cavi_gut.rds")

mb_res_16S_guin_gut = spiec.easi(data = as.matrix(gut_16S_guin), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_guin_gut)
saveRDS(mb_res_16S_guin_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_guin_gut.rds")

#Odontotermes
mb_res_16S_odoA_gut = spiec.easi(data = as.matrix(gut_16S_odoA), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_odoA_gut)
saveRDS(mb_res_16S_odoA_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoA_gut.rds")

mb_res_16S_odoB_gut = spiec.easi(data = as.matrix(gut_16S_odoB), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_odoB_gut)
saveRDS(mb_res_16S_odoB_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoB_gut.rds")

#Pseudocanthotermes
mb_res_16S_mili_gut = spiec.easi(data = as.matrix(gut_16S_mili), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_mili_gut)
saveRDS(mb_res_16S_mili_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_mili_gut.rds")


  ## Comb ##

# Macrotermes
mb_res_16S_belli_comb = spiec.easi(data = as.matrix(comb_16S_beli), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_belli_comb)
saveRDS(mb_res_16S_belli_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_belli_comb.rds")

mb_res_16S_sub_comb = spiec.easi(data = as.matrix(comb_16S_sub), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_sub_comb)
saveRDS(mb_res_16S_sub_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_sub_comb.rds")

#Ancistrotermes
mb_res_16S_cavi_comb = spiec.easi(data = as.matrix(comb_16S_cavi), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_cavi_comb)
saveRDS(mb_res_16S_cavi_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_cavi_comb.rds")

mb_res_16S_guin_comb = spiec.easi(data = as.matrix(comb_16S_guin), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), nlambda = 100)                       
getStability(mb_res_16S_guin_comb)
saveRDS(mb_res_16S_guin_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_guin_comb.rds")

#Odontotermes
mb_res_16S_odoA_comb = spiec.easi(data = as.matrix(comb_16S_odoA), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_odoA_comb)
saveRDS(mb_res_16S_odoA_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoA_comb.rds")

mb_res_16S_odoB_comb = spiec.easi(data = as.matrix(comb_16S_odoB), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_odoB_comb)
saveRDS(mb_res_16S_odoB_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoB_comb.rds")

#Pseudocanthotermes
mb_res_16S_mili_comb = spiec.easi(data = as.matrix(comb_16S_mili), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500))                       
getStability(mb_res_16S_mili_comb)
saveRDS(mb_res_16S_mili_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_mili_comb.rds")
#### 16S network end ####

#### spiec.easi AD networks with "mb" method ####
#===============================================#

  ## Guts ##

# Macrotermes
mb_res_AD_belli_gut = spiec.easi(data = as.matrix(gut_AD_beli), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3)                       
getStability(mb_res_AD_belli_gut)
saveRDS(mb_res_AD_belli_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_belli_gut.rds")

mb_res_AD_sub_gut = spiec.easi(data = as.matrix(gut_AD_sub), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-4, nlambda = 100)                       
getStability(mb_res_AD_sub_gut)
saveRDS(mb_res_AD_sub_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_sub_gut.rds")

#Ancistrotermes
mb_res_AD_cavi_gut = spiec.easi(data = as.matrix(gut_AD_cavi), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_cavi_gut)
saveRDS(mb_res_AD_cavi_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_cavi_gut.rds")

mb_res_AD_guin_gut = spiec.easi(data = as.matrix(gut_AD_guin), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_guin_gut)
saveRDS(mb_res_AD_guin_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_guin_gut.rds")

#Odontotermes
mb_res_AD_odoA_gut = spiec.easi(data = as.matrix(gut_AD_odoA), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_odoA_gut)
saveRDS(mb_res_AD_odoA_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoA_gut.rds")

mb_res_AD_odoB_gut = spiec.easi(data = as.matrix(gut_AD_odoB), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_odoB_gut)
saveRDS(mb_res_AD_odoB_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoB_gut.rds")

#Pseudocanthotermes
mb_res_AD_mili_gut = spiec.easi(data = as.matrix(gut_AD_mili), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_mili_gut)
saveRDS(mb_res_AD_mili_gut, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_mili_gut.rds")


## Comb ##

# Macrotermes
mb_res_AD_belli_comb = spiec.easi(data = as.matrix(comb_AD_beli), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_belli_comb)
saveRDS(mb_res_AD_belli_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_belli_comb.rds")

mb_res_AD_sub_comb = spiec.easi(data = as.matrix(comb_AD_sub), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_sub_comb)
saveRDS(mb_res_AD_sub_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_sub_comb.rds")

#Ancistrotermes
mb_res_AD_cavi_comb = spiec.easi(data = as.matrix(comb_AD_cavi), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 150)                       
getStability(mb_res_AD_cavi_comb)
saveRDS(mb_res_AD_cavi_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_cavi_comb.rds")

mb_res_AD_guin_comb = spiec.easi(data = as.matrix(comb_AD_guin), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 300)                       
getStability(mb_res_AD_guin_comb)
saveRDS(mb_res_AD_guin_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_guin_comb.rds")

#Odontotermes
mb_res_AD_odoA_comb = spiec.easi(data = as.matrix(comb_AD_odoA), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_odoA_comb)
saveRDS(mb_res_AD_odoA_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoA_comb.rds")

mb_res_AD_odoB_comb = spiec.easi(data = as.matrix(comb_AD_odoB), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_odoB_comb)
saveRDS(mb_res_AD_odoB_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoB_comb.rds")

#Pseudocanthotermes
mb_res_AD_mili_comb = spiec.easi(data = as.matrix(comb_AD_mili), method = "mb", pulsar.select = TRUE, sel.criterion = "stars", icov.select.params=list(rep.num=500), lambda.min.ratio = 1e-3, nlambda = 100)                       
getStability(mb_res_AD_mili_comb)
saveRDS(mb_res_AD_mili_comb, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_mili_comb.rds")

#### AD network end ####


#### Read in all networks ####
## 16S ##
  # Guts #
mb_res_16S_belli_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_belli_gut.rds")
mb_res_16S_sub_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_sub_gut.rds")

mb_res_16S_cavi_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_cavi_gut.rds")
mb_res_16S_guin_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_guin_gut.rds")

mb_res_16S_odoA_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoA_gut.rds")
mb_res_16S_odoB_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoB_gut.rds")

mb_res_16S_mili_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_mili_gut.rds")
  
  # Comb #
mb_res_16S_belli_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_belli_comb.rds")
mb_res_16S_sub_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_sub_comb.rds")

mb_res_16S_cavi_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_cavi_comb.rds")
mb_res_16S_guin_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_guin_comb.rds")

mb_res_16S_odoA_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoA_comb.rds")
mb_res_16S_odoB_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_odoB_comb.rds")

mb_res_16S_mili_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/mb_res_16S_mili_comb.rds")

## AD ##
  # Guts #
mb_res_AD_belli_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_belli_gut.rds")
mb_res_AD_sub_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_sub_gut.rds")

mb_res_AD_cavi_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_cavi_gut.rds")
mb_res_AD_guin_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_guin_gut.rds")

mb_res_AD_odoA_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoA_gut.rds")
mb_res_AD_odoB_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoB_gut.rds")

mb_res_AD_mili_gut = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_mili_gut.rds")

  # Comb #
mb_res_AD_belli_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_belli_comb.rds")
mb_res_AD_sub_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_sub_comb.rds")

mb_res_AD_cavi_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_cavi_comb.rds")
mb_res_AD_guin_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_guin_comb.rds")

mb_res_AD_odoA_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoA_comb.rds")
mb_res_AD_odoB_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_odoB_comb.rds")

mb_res_AD_mili_comb = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/mb_res_AD_mili_comb.rds")
#### read networks end ####


#### Network metrics ####
#=======================#

#### Degree ####

## 16S ##

  # Gut #

# Mactorermes
refit_16S_belli_gut = getRefit(mb_res_16S_belli_gut)
deg_dist_16S_belli_gut = igraph::degree(adj2igraph(refit_16S_belli_gut))
deg_dist_16S_belli_gut_df = data.frame(degrees = deg_dist_16S_belli_gut, species = "bellicosus")
mean(igraph::degree(adj2igraph(refit_16S_belli_gut)))
sd(igraph::degree(adj2igraph(refit_16S_belli_gut)))

refit_16S_sub_gut = getRefit(mb_res_16S_sub_gut)
deg_dist_16S_sub_gut = igraph::degree(adj2igraph(refit_16S_sub_gut))
deg_dist_16S_sub_gut_df = data.frame(degrees = deg_dist_16S_sub_gut, species = "subhyalinus")
mean(igraph::degree(adj2igraph(refit_16S_sub_gut)))
sd(igraph::degree(adj2igraph(refit_16S_sub_gut)))

# Ancistrotermes
refit_16S_cavi_gut = getRefit(mb_res_16S_cavi_gut)
deg_dist_16S_cavi_gut = igraph::degree(adj2igraph(refit_16S_cavi_gut))
deg_dist_16S_cavi_gut_df = data.frame(degrees = deg_dist_16S_cavi_gut, species = "cavithorax")
mean(igraph::degree(adj2igraph(refit_16S_cavi_gut)))
sd(igraph::degree(adj2igraph(refit_16S_cavi_gut)))

refit_16S_guin_gut = getRefit(mb_res_16S_guin_gut)
deg_dist_16S_guin_gut = igraph::degree(adj2igraph(refit_16S_guin_gut))
deg_dist_16S_guin_gut_df = data.frame(degrees = deg_dist_16S_guin_gut, species = "guinensis")
mean(igraph::degree(adj2igraph(refit_16S_guin_gut)))
sd(igraph::degree(adj2igraph(refit_16S_guin_gut)))

# Odontotermes
refit_16S_odoA_gut = getRefit(mb_res_16S_odoA_gut)
deg_dist_16S_odoA_gut = igraph::degree(adj2igraph(refit_16S_odoA_gut))
deg_dist_16S_odoA_gut_df = data.frame(degrees = deg_dist_16S_odoA_gut, species = "Species A")
mean(igraph::degree(adj2igraph(refit_16S_odoA_gut)))
sd(igraph::degree(adj2igraph(refit_16S_odoA_gut)))

refit_16S_odoB_gut = getRefit(mb_res_16S_odoB_gut)
deg_dist_16S_odoB_gut = igraph::degree(adj2igraph(refit_16S_odoB_gut))
deg_dist_16S_odoB_gut_df = data.frame(degrees = deg_dist_16S_odoB_gut, species = "Species B")
mean(igraph::degree(adj2igraph(refit_16S_odoB_gut)))
sd(igraph::degree(adj2igraph(refit_16S_odoB_gut)))

# Pseudocanthotermes
refit_16S_mili_gut = getRefit(mb_res_16S_mili_gut)
deg_dist_16S_mili_gut = igraph::degree(adj2igraph(refit_16S_mili_gut))
deg_dist_16S_mili_gut_df = data.frame(degrees = deg_dist_16S_mili_gut, species = "militaris-spiniger")
mean(igraph::degree(adj2igraph(refit_16S_mili_gut)))
sd(igraph::degree(adj2igraph(refit_16S_mili_gut)))

def_16S_gut_df = rbind(deg_dist_16S_belli_gut_df, deg_dist_16S_sub_gut_df, 
                       deg_dist_16S_cavi_gut_df, deg_dist_16S_guin_gut_df, 
                       deg_dist_16S_odoA_gut_df, deg_dist_16S_odoB_gut_df,
                       deg_dist_16S_mili_gut_df)
def_16S_gut_df$species = factor(def_16S_gut_df$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")

gut_16S_degree = ggplot(def_16S_gut_df, aes(x = species, y = degrees, fill = species))+
  geom_violin()+
  theme_pubr()+
  stat_summary(fun.y=mean, geom="point", size=3, colour = "black")+
  scale_fill_manual(values = cols)
  xlab("Species")+
  ylab("Degree")+
  ggtitle("16S gut")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/gut_16S_degree.pdf", gut_16S_degree)

  # Comb #

# Mactorermes
refit_16S_belli_comb = getRefit(mb_res_16S_belli_comb)
deg_dist_16S_belli_comb = igraph::degree(adj2igraph(refit_16S_belli_comb))
deg_dist_16S_belli_comb_df = data.frame(degrees = deg_dist_16S_belli_comb, species = "bellicosus")
mean(igraph::degree(adj2igraph(refit_16S_belli_comb)))
sd(igraph::degree(adj2igraph(refit_16S_belli_comb)))

refit_16S_sub_comb = getRefit(mb_res_16S_sub_comb)
deg_dist_16S_sub_comb = igraph::degree(adj2igraph(refit_16S_sub_comb))
deg_dist_16S_sub_comb_df = data.frame(degrees = deg_dist_16S_sub_comb, species = "subhyalinus")
mean(igraph::degree(adj2igraph(refit_16S_sub_comb)))
sd(igraph::degree(adj2igraph(refit_16S_sub_comb)))

# Ancistrotermes
refit_16S_cavi_comb = getRefit(mb_res_16S_cavi_comb)
deg_dist_16S_cavi_comb = igraph::degree(adj2igraph(refit_16S_cavi_comb))
deg_dist_16S_cavi_comb_df = data.frame(degrees = deg_dist_16S_cavi_comb, species = "cavithorax")
mean(igraph::degree(adj2igraph(refit_16S_cavi_comb)))
sd(igraph::degree(adj2igraph(refit_16S_cavi_comb)))

refit_16S_guin_comb = getRefit(mb_res_16S_guin_comb)
deg_dist_16S_guin_comb = igraph::degree(adj2igraph(refit_16S_guin_comb))
deg_dist_16S_guin_comb_df = data.frame(degrees = deg_dist_16S_guin_comb, species = "guinensis")
mean(igraph::degree(adj2igraph(refit_16S_guin_comb)))
sd(igraph::degree(adj2igraph(refit_16S_guin_comb)))

# Odontotermes
refit_16S_odoA_comb = getRefit(mb_res_16S_odoA_comb)
deg_dist_16S_odoA_comb = igraph::degree(adj2igraph(refit_16S_odoA_comb))
deg_dist_16S_odoA_comb_df = data.frame(degrees = deg_dist_16S_odoA_comb, species = "Species A")
mean(igraph::degree(adj2igraph(refit_16S_odoA_comb)))
sd(igraph::degree(adj2igraph(refit_16S_odoA_comb)))

refit_16S_odoB_comb = getRefit(mb_res_16S_odoB_comb)
deg_dist_16S_odoB_comb = igraph::degree(adj2igraph(refit_16S_odoB_comb))
deg_dist_16S_odoB_comb_df = data.frame(degrees = deg_dist_16S_odoB_comb, species = "Species B")
mean(igraph::degree(adj2igraph(refit_16S_odoB_comb)))
sd(igraph::degree(adj2igraph(refit_16S_odoB_comb)))

# Pseudocanthotermes
refit_16S_mili_comb = getRefit(mb_res_16S_mili_comb)
deg_dist_16S_mili_comb = igraph::degree(adj2igraph(refit_16S_mili_comb))
deg_dist_16S_mili_comb_df = data.frame(degrees = deg_dist_16S_mili_comb, species = "militaris-spiniger")
mean(igraph::degree(adj2igraph(refit_16S_mili_comb)))
sd(igraph::degree(adj2igraph(refit_16S_mili_comb)))


def_16S_comb_df = rbind(deg_dist_16S_belli_comb_df, deg_dist_16S_sub_comb_df, 
                       deg_dist_16S_cavi_comb_df, deg_dist_16S_guin_comb_df, 
                       deg_dist_16S_odoA_comb_df, deg_dist_16S_odoB_comb_df,
                       deg_dist_16S_mili_comb_df)

def_16S_comb_df$species = factor(def_16S_comb_df$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")

comb_16S_degree = ggplot(def_16S_comb_df, aes(x = species, y = degrees, fill = species))+
  geom_violin()+
  theme_pubr()+
  stat_summary(fun.y=mean, geom="point", size=3, colour = "black")+
  scale_fill_manual(values = cols)+
  xlab("Species")+
  ylab("Degree")+
  ggtitle("16S comb")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/comb_16S_degree.pdf", comb_16S_degree)

  ## Combined plot ##

def_16S_gut_df$location = "gut"
def_16S_comb_df$location = "comb"

def_16S_df_plot = rbind(def_16S_gut_df,def_16S_comb_df)

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")
alphas = c("comb" = 0.5, "gut" = 0.8)


combined_16S_degree = ggplot(def_16S_df_plot, aes(x = degrees, y = species, fill = species, alpha = location))+
  geom_density_ridges()+
  theme_ridges()+
  theme(legend.position = "none")+
  scale_fill_manual(values = cols)+
  scale_alpha_manual(values = alphas)+
  xlab("Degree")+
  ylab("Species")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/combined_16S_degree.pdf", combined_16S_degree)


degree_16S_res = aov(def_16S_df_plot$degrees ~ location * species, data = def_16S_df_plot)
summary(degree_16S_res)
cohens_f(degree_16S_res)
TukeyHSD(degree_16S_res, which = "location")


## AD ##

  # Gut #

# Mactorermes
refit_AD_belli_gut = getRefit(mb_res_AD_belli_gut)
deg_dist_AD_belli_gut = igraph::degree(adj2igraph(refit_AD_belli_gut))
deg_dist_AD_belli_gut_df = data.frame(degrees = deg_dist_AD_belli_gut, species = "bellicosus")
mean(igraph::degree(adj2igraph(refit_AD_belli_gut)))
sd(igraph::degree(adj2igraph(refit_AD_belli_gut)))

refit_AD_sub_gut = getRefit(mb_res_AD_sub_gut)
deg_dist_AD_sub_gut = igraph::degree(adj2igraph(refit_AD_sub_gut))
deg_dist_AD_sub_gut_df = data.frame(degrees = deg_dist_AD_sub_gut, species = "subhyalinus")
mean(igraph::degree(adj2igraph(refit_AD_sub_gut)))
sd(igraph::degree(adj2igraph(refit_AD_sub_gut)))

# Ancistrotermes
refit_AD_cavi_gut = getRefit(mb_res_AD_cavi_gut)
deg_dist_AD_cavi_gut = igraph::degree(adj2igraph(refit_AD_cavi_gut))
deg_dist_AD_cavi_gut_df = data.frame(degrees = deg_dist_AD_cavi_gut, species = "cavithorax")
mean(igraph::degree(adj2igraph(refit_AD_cavi_gut)))
sd(igraph::degree(adj2igraph(refit_AD_cavi_gut)))

refit_AD_guin_gut = getRefit(mb_res_AD_guin_gut)
deg_dist_AD_guin_gut = igraph::degree(adj2igraph(refit_AD_guin_gut))
deg_dist_AD_guin_gut_df = data.frame(degrees = deg_dist_AD_guin_gut, species = "guinensis")
mean(igraph::degree(adj2igraph(refit_AD_guin_gut)))
sd(igraph::degree(adj2igraph(refit_AD_guin_gut)))

# Odontotermes
refit_AD_odoA_gut = getRefit(mb_res_AD_odoA_gut)
deg_dist_AD_odoA_gut = igraph::degree(adj2igraph(refit_AD_odoA_gut))
deg_dist_AD_odoA_gut_df = data.frame(degrees = deg_dist_AD_odoA_gut, species = "Species A")
mean(igraph::degree(adj2igraph(refit_AD_odoA_gut)))
sd(igraph::degree(adj2igraph(refit_AD_odoA_gut)))

refit_AD_odoB_gut = getRefit(mb_res_AD_odoB_gut)
deg_dist_AD_odoB_gut = igraph::degree(adj2igraph(refit_AD_odoB_gut))
deg_dist_AD_odoB_gut_df = data.frame(degrees = deg_dist_AD_odoB_gut, species = "Species B")
mean(igraph::degree(adj2igraph(refit_AD_odoB_gut)))
sd(igraph::degree(adj2igraph(refit_AD_odoB_gut)))

# Pseudocanthotermes
refit_AD_mili_gut = getRefit(mb_res_AD_mili_gut)
deg_dist_AD_mili_gut = igraph::degree(adj2igraph(refit_AD_mili_gut))
deg_dist_AD_mili_gut_df = data.frame(degrees = deg_dist_AD_mili_gut, species = "militaris-spiniger")
mean(igraph::degree(adj2igraph(refit_AD_mili_gut)))
sd(igraph::degree(adj2igraph(refit_AD_mili_gut)))

def_AD_gut_df = rbind(deg_dist_AD_belli_gut_df, deg_dist_AD_sub_gut_df, 
                       deg_dist_AD_cavi_gut_df, deg_dist_AD_guin_gut_df, 
                       deg_dist_AD_odoA_gut_df, deg_dist_AD_odoB_gut_df,
                       deg_dist_AD_mili_gut_df)
def_AD_gut_df$species = factor(def_AD_gut_df$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")

gut_AD_degree = ggplot(def_AD_gut_df, aes(x = species, y = degrees, fill = species))+
  geom_violin()+
  theme_pubr()+
  stat_summary(fun.y=mean, geom="point", size=3, colour = "black")+
  scale_fill_manual(values = cols)
xlab("Species")+
  ylab("Degree")+
  ggtitle("AD gut")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/gut_AD_degree.pdf", gut_AD_degree)

# Comb #

# Mactorermes
refit_AD_belli_comb = getRefit(mb_res_AD_belli_comb)
deg_dist_AD_belli_comb = igraph::degree(adj2igraph(refit_AD_belli_comb))
deg_dist_AD_belli_comb_df = data.frame(degrees = deg_dist_AD_belli_comb, species = "bellicosus")
mean(igraph::degree(adj2igraph(refit_AD_belli_comb)))
sd(igraph::degree(adj2igraph(refit_AD_belli_comb)))

refit_AD_sub_comb = getRefit(mb_res_AD_sub_comb)
deg_dist_AD_sub_comb = igraph::degree(adj2igraph(refit_AD_sub_comb))
deg_dist_AD_sub_comb_df = data.frame(degrees = deg_dist_AD_sub_comb, species = "subhyalinus")
mean(igraph::degree(adj2igraph(refit_AD_sub_comb)))
sd(igraph::degree(adj2igraph(refit_AD_sub_comb)))

# Ancistrotermes
refit_AD_cavi_comb = getRefit(mb_res_AD_cavi_comb)
deg_dist_AD_cavi_comb = igraph::degree(adj2igraph(refit_AD_cavi_comb))
deg_dist_AD_cavi_comb_df = data.frame(degrees = deg_dist_AD_cavi_comb, species = "cavithorax")
mean(igraph::degree(adj2igraph(refit_AD_cavi_comb)))
sd(igraph::degree(adj2igraph(refit_AD_cavi_comb)))

refit_AD_guin_comb = getRefit(mb_res_AD_guin_comb)
deg_dist_AD_guin_comb = igraph::degree(adj2igraph(refit_AD_guin_comb))
deg_dist_AD_guin_comb_df = data.frame(degrees = deg_dist_AD_guin_comb, species = "guinensis")
mean(igraph::degree(adj2igraph(refit_AD_guin_comb)))
sd(igraph::degree(adj2igraph(refit_AD_guin_comb)))

# Odontotermes
refit_AD_odoA_comb = getRefit(mb_res_AD_odoA_comb)
deg_dist_AD_odoA_comb = igraph::degree(adj2igraph(refit_AD_odoA_comb))
deg_dist_AD_odoA_comb_df = data.frame(degrees = deg_dist_AD_odoA_comb, species = "Species A")
mean(igraph::degree(adj2igraph(refit_AD_odoA_comb)))
sd(igraph::degree(adj2igraph(refit_AD_odoA_comb)))

refit_AD_odoB_comb = getRefit(mb_res_AD_odoB_comb)
deg_dist_AD_odoB_comb = igraph::degree(adj2igraph(refit_AD_odoB_comb))
deg_dist_AD_odoB_comb_df = data.frame(degrees = deg_dist_AD_odoB_comb, species = "Species B")
mean(igraph::degree(adj2igraph(refit_AD_odoB_comb)))
sd(igraph::degree(adj2igraph(refit_AD_odoB_comb)))

# Pseudocanthotermes
refit_AD_mili_comb = getRefit(mb_res_AD_mili_comb)
deg_dist_AD_mili_comb = igraph::degree(adj2igraph(refit_AD_mili_comb))
deg_dist_AD_mili_comb_df = data.frame(degrees = deg_dist_AD_mili_comb, species = "militaris-spiniger")
mean(igraph::degree(adj2igraph(refit_AD_mili_comb)))
sd(igraph::degree(adj2igraph(refit_AD_mili_comb)))


def_AD_comb_df = rbind(deg_dist_AD_belli_comb_df, deg_dist_AD_sub_comb_df, 
                        deg_dist_AD_cavi_comb_df, deg_dist_AD_guin_comb_df, 
                        deg_dist_AD_odoA_comb_df, deg_dist_AD_odoB_comb_df,
                        deg_dist_AD_mili_comb_df)
def_AD_comb_df$species = factor(def_AD_comb_df$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23",
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")

comb_AD_degree = ggplot(def_AD_comb_df, aes(x = species, y = degrees, fill = species))+
  geom_violin()+
  theme_pubr()+
  stat_summary(fun.y=mean, geom="point", size=3, colour = "black")+
  scale_fill_manual(values = cols)+
  xlab("Species")+
  ylab("Degree")+
  ggtitle("AD comb")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/comb_AD_degree.pdf", comb_AD_degree)

## Combined plot ##

def_AD_gut_df$location = "gut"
def_AD_comb_df$location = "comb"

def_AD_df_plot = rbind(def_AD_gut_df,def_AD_comb_df)

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")
alphas = c("comb" = 0.5, "gut" = 0.8)


combined_AD_degree = ggplot(def_AD_df_plot, aes(x = degrees, y = species, fill = species, alpha = location))+
  geom_density_ridges()+
  theme_ridges()+
  theme(legend.position = "none")+
  scale_fill_manual(values = cols)+
  scale_alpha_manual(values = alphas)+
  xlab("Species")+
  ylab("Degree")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/combined_AD_degree.pdf", combined_AD_degree)

degree_AD_res = aov(def_AD_df_plot$degrees ~ location * species, data = def_AD_df_plot)
summary(degree_AD_res)
cohens_f(degree_AD_res)
TukeyHSD(degree_AD_res, which = "location")


#### degree end ####

#### Degree density ####

## 16S ##

  # Gut #

# Macrotermes
refit_16S_belli_gut = getRefit(mb_res_16S_belli_gut)
belli_gut_16S_density = c(species = "bellicosus", location = "gut", density = graph.density(adj2igraph(refit_16S_belli_gut)))

refit_16S_sub_gut = getRefit(mb_res_16S_sub_gut)
sub_gut_16S_density = c(species = "subhyalinus", location = "gut", density = graph.density(adj2igraph(refit_16S_sub_gut)))

# Ancistrotermes
refit_16S_cavi_gut = getRefit(mb_res_16S_cavi_gut)
cavi_gut_16S_density = c(species = "cavithorax", location = "gut", density = graph.density(adj2igraph(refit_16S_cavi_gut)))

refit_16S_guin_gut = getRefit(mb_res_16S_guin_gut)
guin_gut_16S_density = c(species = "guinensis", location = "gut", density = graph.density(adj2igraph(refit_16S_guin_gut)))
# Odontotermes
refit_16S_odoA_gut = getRefit(mb_res_16S_odoA_gut)
odoA_gut_16S_density = c(species = "Species A", location = "gut", density = graph.density(adj2igraph(refit_16S_odoA_gut)))

refit_16S_odoB_gut = getRefit(mb_res_16S_odoB_gut)
oboB_gut_16S_density = c(species = "Species B", location = "gut", density = graph.density(adj2igraph(refit_16S_odoB_gut)))

# Pseudocanthotermes
refit_16S_mili_gut = getRefit(mb_res_16S_mili_gut)
mili_gut_16S_density = c(species = "militaris-spiniger", location = "gut", density = graph.density(adj2igraph(refit_16S_mili_gut)))

deg_density_16S_gut = as.data.frame(rbind(belli_gut_16S_density, sub_gut_16S_density, 
                       cavi_gut_16S_density, guin_gut_16S_density, 
                       odoA_gut_16S_density, oboB_gut_16S_density,
                       mili_gut_16S_density))
deg_density_16S_gut$species = factor(deg_density_16S_gut$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))
deg_density_16S_gut$density = as.numeric(deg_density_16S_gut$density)
  # Comb #

# Macrotermes
refit_16S_belli_comb = getRefit(mb_res_16S_belli_comb)
belli_comb_16S_density = c(species = "bellicosus", location = "comb", density = graph.density(adj2igraph(refit_16S_belli_comb)))

refit_16S_sub_comb = getRefit(mb_res_16S_sub_comb)
sub_comb_16S_density = c(species = "subhyalinus", location = "comb", density = graph.density(adj2igraph(refit_16S_sub_comb)))

# Ancistrotermes
refit_16S_cavi_comb = getRefit(mb_res_16S_cavi_comb)
cavi_comb_16S_density = c(species = "cavithorax", location = "comb", density = graph.density(adj2igraph(refit_16S_cavi_comb)))

refit_16S_guin_comb = getRefit(mb_res_16S_guin_comb)
guin_comb_16S_density = c(species = "guinensis", location = "comb", density = graph.density(adj2igraph(refit_16S_guin_comb)))
# Odontotermes
refit_16S_odoA_comb = getRefit(mb_res_16S_odoA_comb)
odoA_comb_16S_density = c(species = "Species A", location = "comb", density = graph.density(adj2igraph(refit_16S_odoA_comb)))

refit_16S_odoB_comb = getRefit(mb_res_16S_odoB_comb)
oboB_comb_16S_density = c(species = "Species B", location = "comb", density = graph.density(adj2igraph(refit_16S_odoB_comb)))

# Pseudocanthotermes
refit_16S_mili_comb = getRefit(mb_res_16S_mili_comb)
mili_comb_16S_density = c(species = "militaris-spiniger", location = "comb", density = graph.density(adj2igraph(refit_16S_mili_comb)))

deg_density_16S_comb = as.data.frame(rbind(belli_comb_16S_density, sub_comb_16S_density, 
                                          cavi_comb_16S_density, guin_comb_16S_density, 
                                          odoA_comb_16S_density, oboB_comb_16S_density,
                                          mili_comb_16S_density))
deg_density_16S_comb$species = factor(deg_density_16S_comb$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))
deg_density_16S_comb$density = as.numeric(deg_density_16S_comb$density)

## Combined plot ##

deg_density_df = rbind(deg_density_16S_gut,deg_density_16S_comb)

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")
alphas = c("comb" = 0.4, "gut" = 0.9)


combined_16S_degree_density = ggplot(deg_density_df, aes(x = species, y = density, fill = species, alpha = location))+
  geom_bar(stat = "identity", position = "dodge")+
  coord_flip()+
  theme_pubr()+
  theme(legend.position = "none")+
  scale_fill_manual(values = cols)+
  scale_alpha_manual(values = alphas)+
  xlab("Species")+
  ylab("Degree density")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/degree_density_16S.pdf", combined_16S_degree_density, height = 9.74, width = 4)


degree_16S_res = aov(def_16S_df_plot$degrees ~ location * species, data = def_16S_df_plot)
summary(degree_16S_res)
cohens_f(degree_16S_res)
TukeyHSD(degree_16S_res, which = "location")

## AD ##

# Gut #

## AD ##

# Gut #

# Macrotermes
refit_AD_belli_gut = getRefit(mb_res_AD_belli_gut)
belli_gut_AD_density = c(species = "bellicosus", location = "gut", density = graph.density(adj2igraph(refit_AD_belli_gut)))

refit_AD_sub_gut = getRefit(mb_res_AD_sub_gut)
sub_gut_AD_density = c(species = "subhyalinus", location = "gut", density = graph.density(adj2igraph(refit_AD_sub_gut)))

# Ancistrotermes
refit_AD_cavi_gut = getRefit(mb_res_AD_cavi_gut)
cavi_gut_AD_density = c(species = "cavithorax", location = "gut", density = graph.density(adj2igraph(refit_AD_cavi_gut)))

refit_AD_guin_gut = getRefit(mb_res_AD_guin_gut)
guin_gut_AD_density = c(species = "guinensis", location = "gut", density = graph.density(adj2igraph(refit_AD_guin_gut)))
# Odontotermes
refit_AD_odoA_gut = getRefit(mb_res_AD_odoA_gut)
odoA_gut_AD_density = c(species = "Species A", location = "gut", density = graph.density(adj2igraph(refit_AD_odoA_gut)))

refit_AD_odoB_gut = getRefit(mb_res_AD_odoB_gut)
oboB_gut_AD_density = c(species = "Species B", location = "gut", density = graph.density(adj2igraph(refit_AD_odoB_gut)))

# Pseudocanthotermes
refit_AD_mili_gut = getRefit(mb_res_AD_mili_gut)
mili_gut_AD_density = c(species = "militaris-spiniger", location = "gut", density = graph.density(adj2igraph(refit_AD_mili_gut)))

deg_density_AD_gut = as.data.frame(rbind(belli_gut_AD_density, sub_gut_AD_density, 
                                          cavi_gut_AD_density, guin_gut_AD_density, 
                                          odoA_gut_AD_density, oboB_gut_AD_density,
                                          mili_gut_AD_density))
deg_density_AD_gut$species = factor(deg_density_AD_gut$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))
deg_density_AD_gut$density = as.numeric(deg_density_AD_gut$density)
# Comb #

# Macrotermes
refit_AD_belli_comb = getRefit(mb_res_AD_belli_comb)
belli_comb_AD_density = c(species = "bellicosus", location = "comb", density = graph.density(adj2igraph(refit_AD_belli_comb)))

refit_AD_sub_comb = getRefit(mb_res_AD_sub_comb)
sub_comb_AD_density = c(species = "subhyalinus", location = "comb", density = graph.density(adj2igraph(refit_AD_sub_comb)))

# Ancistrotermes
refit_AD_cavi_comb = getRefit(mb_res_AD_cavi_comb)
cavi_comb_AD_density = c(species = "cavithorax", location = "comb", density = graph.density(adj2igraph(refit_AD_cavi_comb)))

refit_AD_guin_comb = getRefit(mb_res_AD_guin_comb)
guin_comb_AD_density = c(species = "guinensis", location = "comb", density = graph.density(adj2igraph(refit_AD_guin_comb)))
# Odontotermes
refit_AD_odoA_comb = getRefit(mb_res_AD_odoA_comb)
odoA_comb_AD_density = c(species = "Species A", location = "comb", density = graph.density(adj2igraph(refit_AD_odoA_comb)))

refit_AD_odoB_comb = getRefit(mb_res_AD_odoB_comb)
oboB_comb_AD_density = c(species = "Species B", location = "comb", density = graph.density(adj2igraph(refit_AD_odoB_comb)))

# Pseudocanthotermes
refit_AD_mili_comb = getRefit(mb_res_AD_mili_comb)
mili_comb_AD_density = c(species = "militaris-spiniger", location = "comb", density = graph.density(adj2igraph(refit_AD_mili_comb)))

deg_density_AD_comb = as.data.frame(rbind(belli_comb_AD_density, sub_comb_AD_density, 
                                           cavi_comb_AD_density, guin_comb_AD_density, 
                                           odoA_comb_AD_density, oboB_comb_AD_density,
                                           mili_comb_AD_density))
deg_density_AD_comb$species = factor(deg_density_AD_comb$species, levels = c("militaris-spiniger", "Species B", "Species A", "guinensis", "cavithorax", "subhyalinus", "bellicosus"))
deg_density_AD_comb$density = as.numeric(deg_density_AD_comb$density)

## Combined plot ##

deg_density_df = rbind(deg_density_AD_gut,deg_density_AD_comb)

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")
alphas = c("comb" = 0.4, "gut" = 0.9)


combined_AD_degree_density = ggplot(deg_density_df, aes(x = species, y = density, fill = species, alpha = location))+
  geom_bar(stat = "identity", position = "dodge")+
  coord_flip()+
  theme_pubr()+
  theme(legend.position = "none")+
  scale_fill_manual(values = cols)+
  scale_alpha_manual(values = alphas)+
  xlab("Species")+
  ylab("Degree density")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/degree_density_AD.pdf", combined_AD_degree_density, height = 9.74, width = 4)


#### Degree density end ####

#### Edge weights ####

## 16S ##

  # Gut #

edges_df_gut_16S = data.frame(species = c("bellicosus_16S_gut", "subhyalinus_16S_gut", 
                                          "cavithorax_16S_gut", "guinensis_16S_gut", 
                                          "Species A_16S_gut", "Species B_16S_gut", 
                                          "militaris_16S_gut"),
                              percent = NA)
# Macrotermes
cor_16S_belli_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_belli_gut), mode='maxabs'))[,3], sample = "bellicosus_16S")
edges_df_gut_16S[edges_df_gut_16S$species == "bellicosus_16S_gut",]$percent = sum(cor_16S_belli_gut$edges > 0)/nrow(cor_16S_belli_gut)

cor_16S_sub_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_sub_gut), mode='maxabs'))[,3], sample = "subhyalinus_16S")
edges_df_gut_16S[edges_df_gut_16S$species == "subhyalinus_16S_gut",]$percent = sum(cor_16S_sub_gut$edges > 0)/nrow(cor_16S_sub_gut)

# Ancistro
cor_16S_cavi_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_cavi_gut), mode='maxabs'))[,3], sample = "cavithorax_16S")
edges_df_gut_16S[edges_df_gut_16S$species == "cavithorax_16S_gut",]$percent = sum(cor_16S_cavi_gut$edges > 0)/nrow(cor_16S_cavi_gut)

cor_16S_guin_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_guin_gut), mode='maxabs'))[,3], sample = "guinensis_16S")
edges_df_gut_16S[edges_df_gut_16S$species == "guinensis_16S_gut",]$percent = sum(cor_16S_guin_gut$edges > 0)/nrow(cor_16S_guin_gut)

# Odontotermes
cor_16S_odoA_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_odoA_gut), mode='maxabs'))[,3], sample = "Species A_16S")
edges_df_gut_16S[edges_df_gut_16S$species == "Species A_16S_gut",]$percent = sum(cor_16S_odoA_gut$edges > 0)/nrow(cor_16S_odoA_gut)

cor_16S_odoB_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_odoB_gut), mode='maxabs'))[,3], sample = "Species B_16S")
edges_df_gut_16S[edges_df_gut_16S$species == "Species B_16S_gut",]$percent = sum(cor_16S_odoB_gut$edges > 0)/nrow(cor_16S_odoB_gut)

# Pseudo
cor_16S_mili_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_mili_gut), mode='maxabs'))[,3], sample = "militaris_16S")
edges_df_gut_16S[edges_df_gut_16S$species == "militaris_16S_gut",]$percent = sum(cor_16S_mili_gut$edges > 0)/nrow(cor_16S_mili_gut)


  # Comb #

edges_df_comb_16S = data.frame(species = c("bellicosus_16S_comb", "subhyalinus_16S_comb", 
                                           "cavithorax_16S_comb", "guinensis_16S_comb", 
                                           "Species A_16S_comb", "Species B_16S_comb", 
                                           "militaris_16S_comb"),
                              percent = NA)
# Macrotermes
cor_16S_belli_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_belli_comb), mode='maxabs'))[,3], sample = "bellicosus_16S")
edges_df_comb_16S[edges_df_comb_16S$species == "bellicosus_16S_comb",]$percent = sum(cor_16S_belli_comb$edges > 0)/nrow(cor_16S_belli_comb)

cor_16S_sub_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_sub_comb), mode='maxabs'))[,3], sample = "subhyalinus_16S")
edges_df_comb_16S[edges_df_comb_16S$species == "subhyalinus_16S_comb",]$percent = sum(cor_16S_sub_comb$edges > 0)/nrow(cor_16S_sub_comb)

# Ancistro
cor_16S_cavi_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_cavi_comb), mode='maxabs'))[,3], sample = "cavithorax_16S")
edges_df_comb_16S[edges_df_comb_16S$species == "cavithorax_16S_comb",]$percent = sum(cor_16S_cavi_comb$edges > 0)/nrow(cor_16S_cavi_comb)

cor_16S_guin_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_guin_comb), mode='maxabs'))[,3], sample = "guinensis_16S")
edges_df_comb_16S[edges_df_comb_16S$species == "guinensis_16S_comb",]$percent = sum(cor_16S_guin_comb$edges > 0)/nrow(cor_16S_guin_comb)

# Odontotermes
cor_16S_odoA_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_odoA_comb), mode='maxabs'))[,3], sample = "Species A_16S")
edges_df_comb_16S[edges_df_comb_16S$species == "Species A_16S_comb",]$percent = sum(cor_16S_odoA_comb$edges > 0)/nrow(cor_16S_odoA_comb)

cor_16S_odoB_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_odoB_comb), mode='maxabs'))[,3], sample = "Species B_16S")
edges_df_comb_16S[edges_df_comb_16S$species == "Species B_16S_comb",]$percent = sum(cor_16S_odoB_comb$edges > 0)/nrow(cor_16S_odoB_comb)

# Pseudo
cor_16S_mili_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_16S_mili_comb), mode='maxabs'))[,3], sample = "militaris_16S")
edges_df_comb_16S[edges_df_comb_16S$species == "militaris_16S_comb",]$percent = sum(cor_16S_mili_comb$edges > 0)/nrow(cor_16S_mili_comb)




## AD ##

  # Gut #

edges_df_gut_AD = data.frame(species = c("bellicosus_AD_gut", "subhyalinus_AD_gut", 
                                          "cavithorax_AD_gut", "guinensis_AD_gut", 
                                          "Species A_AD_gut", "Species B_AD_gut", 
                                          "militaris_AD_gut"),
                              percent = NA)
# Macrotermes
cor_AD_belli_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_belli_gut), mode='maxabs'))[,3], sample = "bellicosus_AD")
edges_df_gut_AD[edges_df_gut_AD$species == "bellicosus_AD_gut",]$percent = sum(cor_AD_belli_gut$edges > 0)/nrow(cor_AD_belli_gut)

cor_AD_sub_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_sub_gut), mode='maxabs'))[,3], sample = "subhyalinus_AD")
edges_df_gut_AD[edges_df_gut_AD$species == "subhyalinus_AD_gut",]$percent = sum(cor_AD_sub_gut$edges > 0)/nrow(cor_AD_sub_gut)

# Ancistro
cor_AD_cavi_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_cavi_gut), mode='maxabs'))[,3], sample = "cavithorax_AD")
edges_df_gut_AD[edges_df_gut_AD$species == "cavithorax_AD_gut",]$percent = sum(cor_AD_cavi_gut$edges > 0)/nrow(cor_AD_cavi_gut)

cor_AD_guin_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_guin_gut), mode='maxabs'))[,3], sample = "guinensis_AD")
edges_df_gut_AD[edges_df_gut_AD$species == "guinensis_AD_gut",]$percent = sum(cor_AD_guin_gut$edges > 0)/nrow(cor_AD_guin_gut)

# Odontotermes
cor_AD_odoA_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_odoA_gut), mode='maxabs'))[,3], sample = "Species A_AD")
edges_df_gut_AD[edges_df_gut_AD$species == "Species A_AD_gut",]$percent = sum(cor_AD_odoA_gut$edges > 0)/nrow(cor_AD_odoA_gut)

cor_AD_odoB_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_odoB_gut), mode='maxabs'))[,3], sample = "Species B_AD")
edges_df_gut_AD[edges_df_gut_AD$species == "Species B_AD_gut",]$percent = sum(cor_AD_odoB_gut$edges > 0)/nrow(cor_AD_odoB_gut)

# Pseudo
cor_AD_mili_gut = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_mili_gut), mode='maxabs'))[,3], sample = "militaris_AD")
edges_df_gut_AD[edges_df_gut_AD$species == "militaris_AD_gut",]$percent = sum(cor_AD_mili_gut$edges > 0)/nrow(cor_AD_mili_gut)


  # Comb #

edges_df_comb_AD = data.frame(species = c("bellicosus_AD_comb", "subhyalinus_AD_comb", 
                                           "cavithorax_AD_comb", "guinensis_AD_comb", 
                                           "Species A_AD_comb", "Species B_AD_comb", 
                                           "militaris_AD_comb"),
                               percent = NA)
# Macrotermes
cor_AD_belli_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_belli_comb), mode='maxabs'))[,3], sample = "bellicosus_AD")
edges_df_comb_AD[edges_df_comb_AD$species == "bellicosus_AD_comb",]$percent = sum(cor_AD_belli_comb$edges > 0)/nrow(cor_AD_belli_comb)

cor_AD_sub_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_sub_comb), mode='maxabs'))[,3], sample = "subhyalinus_AD")
edges_df_comb_AD[edges_df_comb_AD$species == "subhyalinus_AD_comb",]$percent = sum(cor_AD_sub_comb$edges > 0)/nrow(cor_AD_sub_comb)

# Ancistro
cor_AD_cavi_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_cavi_comb), mode='maxabs'))[,3], sample = "cavithorax_AD")
edges_df_comb_AD[edges_df_comb_AD$species == "cavithorax_AD_comb",]$percent = sum(cor_AD_cavi_comb$edges > 0)/nrow(cor_AD_cavi_comb)

cor_AD_guin_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_guin_comb), mode='maxabs'))[,3], sample = "guinensis_AD")
edges_df_comb_AD[edges_df_comb_AD$species == "guinensis_AD_comb",]$percent = sum(cor_AD_guin_comb$edges > 0)/nrow(cor_AD_guin_comb)

# Odontotermes
cor_AD_odoA_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_odoA_comb), mode='maxabs'))[,3], sample = "Species A_AD")
edges_df_comb_AD[edges_df_comb_AD$species == "Species A_AD_comb",]$percent = sum(cor_AD_odoA_comb$edges > 0)/nrow(cor_AD_odoA_comb)

cor_AD_odoB_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_odoB_comb), mode='maxabs'))[,3], sample = "Species B_AD")
edges_df_comb_AD[edges_df_comb_AD$species == "Species B_AD_comb",]$percent = sum(cor_AD_odoB_comb$edges > 0)/nrow(cor_AD_odoB_comb)

# Pseudo
cor_AD_mili_comb = data.frame(edges = summary(symBeta(getOptBeta(mb_res_AD_mili_comb), mode='maxabs'))[,3], sample = "militaris_AD")
edges_df_comb_AD[edges_df_comb_AD$species == "militaris_AD_comb",]$percent = sum(cor_AD_mili_comb$edges > 0)/nrow(cor_AD_mili_comb)

## Plot edge weights ##

edge_cols = c("bellicosus_16S_comb" = "#E77675", "subhyalinus_16S_comb" = "#A31E23", "cavithorax_16S_comb" = "#E3A86C", "guinensis_16S_comb" = "#BA7328",
              "Species A_16S_comb" = "#91CB7B", "Species B_16S_comb" = "#3B7538", "militaris_16S_comb" = "#8499C7", 
              "bellicosus_16S_gut" = "#E77675", "subhyalinus_16S_gut" = "#A31E23", "cavithorax_16S_gut" = "#E3A86C", "guinensis_16S_gut" = "#BA7328",
              "Species A_16S_gut" = "#91CB7B","Species B_16S_gut" = "#3B7538", "militaris_16S_gut" = "#8499C7")
# Gut #

weight_df_16S= rbind(edges_df_gut_16S, edges_df_comb_16S)
weight_df_16S$species = factor(weight_df_16S$species, levels = c("bellicosus_16S_comb", "subhyalinus_16S_comb", "cavithorax_16S_comb", "guinensis_16S_comb",
                                                                 "Species A_16S_comb", "Species B_16S_comb", "militaris_16S_comb", 
                                                                 "militaris_16S_gut", "Species B_16S_gut", "Species A_16S_gut", "guinensis_16S_gut", 
                                                                 "cavithorax_16S_gut", "subhyalinus_16S_gut", "bellicosus_16S_gut"))

label_data = weight_df_16S %>%
  slice(match(c("bellicosus_16S_comb", "subhyalinus_16S_comb", "cavithorax_16S_comb", "guinensis_16S_comb", "Species A_16S_comb", 
                "Species B_16S_comb", "militaris_16S_comb", 
                "militaris_16S_gut", "Species B_16S_gut", "Species A_16S_gut", "guinensis_16S_gut", "cavithorax_16S_gut",
                "subhyalinus_16S_gut", "bellicosus_16S_gut"), species)) %>%
  rownames_to_column(var = "id")
number_of_bar = nrow(label_data)
angle =  90 - 360 * (as.numeric(label_data$id)-0.5) /number_of_bar
label_data$hjust = ifelse(angle < -90, 1, 0)
label_data$angle = ifelse(angle < -90, angle+180, angle)

edge_plot_16S = ggplot(weight_df_16S, aes(y = percent*100, x = species, fill = species))+
  geom_bar(stat="identity")+
  geom_text(data = label_data, aes(x = as.numeric(id), y = (percent*100)+5, label = round(100*percent, 1), 
                                   hjust = hjust), angle = label_data$angle, size = 4, inherit.aes = FALSE,fontface="bold")+
  coord_polar(start = 0)+
  ylim(-30, 100)+
  scale_fill_manual(values = edge_cols)+
  theme_pubclean()+
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank())

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/edge_plot_16S.pdf", edge_plot_16S)



# Comb #

edge_cols = c("bellicosus_AD_comb" = "#E77675", "subhyalinus_AD_comb" = "#A31E23", "cavithorax_AD_comb" = "#E3A86C", "guinensis_AD_comb" = "#BA7328",
              "Species A_AD_comb" = "#91CB7B", "Species B_AD_comb" = "#3B7538", "militaris_AD_comb" = "#8499C7", 
              "bellicosus_AD_gut" = "#E77675", "subhyalinus_AD_gut" = "#A31E23", "cavithorax_AD_gut" = "#E3A86C", "guinensis_AD_gut" = "#BA7328",
              "Species A_AD_gut" = "#91CB7B","Species B_AD_gut" = "#3B7538", "militaris_AD_gut" = "#8499C7")

weight_df_AD= rbind(edges_df_gut_AD, edges_df_comb_AD)
weight_df_AD$species = factor(weight_df_AD$species, levels = c("bellicosus_AD_comb", "subhyalinus_AD_comb", "cavithorax_AD_comb", "guinensis_AD_comb",
                                                                 "Species A_AD_comb", "Species B_AD_comb", "militaris_AD_comb", 
                                                                 "militaris_AD_gut", "Species B_AD_gut", "Species A_AD_gut", "guinensis_AD_gut", 
                                                                 "cavithorax_AD_gut", "subhyalinus_AD_gut", "bellicosus_AD_gut"))

label_data = weight_df_AD %>%
  slice(match(c("bellicosus_AD_comb", "subhyalinus_AD_comb", "cavithorax_AD_comb", "guinensis_AD_comb", "Species A_AD_comb", 
                "Species B_AD_comb", "militaris_AD_comb", 
                "militaris_AD_gut", "Species B_AD_gut", "Species A_AD_gut", "guinensis_AD_gut", "cavithorax_AD_gut",
                "subhyalinus_AD_gut", "bellicosus_AD_gut"), species)) %>%
  rownames_to_column(var = "id")
number_of_bar = nrow(label_data)
angle =  90 - 360 * (as.numeric(label_data$id)-0.5) /number_of_bar
label_data$hjust = ifelse(angle < -90, 1, 0)
label_data$angle = ifelse(angle < -90, angle+180, angle)

edge_plot_AD = ggplot(weight_df_AD, aes(y = percent*100, x = species, fill = species))+
  geom_bar(stat="identity")+
  geom_text(data = label_data, aes(x = as.numeric(id), y = (percent*100)+5, label = round(100*percent, 1), 
                                   hjust = hjust), angle = label_data$angle, size = 4, inherit.aes = FALSE,fontface="bold")+
  coord_polar(start = 0)+
  ylim(-30, 110)+
  scale_fill_manual(values = edge_cols)+
  theme_pubclean()+
  theme(legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank())


ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/edge_plot_AD.pdf", edge_plot_AD)

#### Edge weights end ####

#### Modularity ####

## 16S ##
#=======#

  # Gut #
modularity_df_gut_16S = data.frame(species = c("bellicosus_16S_gut", "subhyalinus_16S_gut", "cavithorax_16S_gut", "guinensis_16S_gut", "Species A_16S_gut",
                                               "Species B_16S_gut", "militaris_16S_gut"),
                                   modularity = NA)
# Macrotermes
refit_16S_belli_gut = getRefit(mb_res_16S_belli_gut)
graph_16S_belli_gut = graph.adjacency(refit_16S_belli_gut, mode = "undirected", weighted = TRUE)
modules_16S_belli_gut = cluster_fast_greedy(graph_16S_belli_gut)
modularity_df_gut_16S[modularity_df_gut_16S$species == "bellicosus_16S_gut",]$modularity = modularity(modules_16S_belli_gut)

refit_16S_sub_gut = getRefit(mb_res_16S_sub_gut)
graph_16S_sub_gut = graph.adjacency(refit_16S_sub_gut, mode = "undirected", weighted = TRUE)
modules_16S_sub_gut = cluster_fast_greedy(graph_16S_sub_gut)
modularity_df_gut_16S[modularity_df_gut_16S$species == "subhyalinus_16S_gut",]$modularity = modularity(modules_16S_sub_gut)

# Ancistrotermes
refit_16S_cavi_gut = getRefit(mb_res_16S_cavi_gut)
graph_16S_cavi_gut = graph.adjacency(refit_16S_cavi_gut, mode = "undirected", weighted = TRUE)
modules_16S_cavi_gut = cluster_fast_greedy(graph_16S_cavi_gut)
modularity_df_gut_16S[modularity_df_gut_16S$species == "cavithorax_16S_gut",]$modularity = modularity(modules_16S_cavi_gut)

refit_16S_guin_gut = getRefit(mb_res_16S_guin_gut)
graph_16S_guin_gut = graph.adjacency(refit_16S_guin_gut, mode = "undirected", weighted = TRUE)
modules_16S_guin_gut = cluster_fast_greedy(graph_16S_guin_gut)
modularity_df_gut_16S[modularity_df_gut_16S$species == "guinensis_16S_gut",]$modularity = modularity(modules_16S_guin_gut)

# Odontotermes
refit_16S_odoA_gut = getRefit(mb_res_16S_odoA_gut)
graph_16S_odoA_gut = graph.adjacency(refit_16S_odoA_gut, mode = "undirected", weighted = TRUE)
modules_16S_odoA_gut = cluster_fast_greedy(graph_16S_odoA_gut)
modularity_df_gut_16S[modularity_df_gut_16S$species == "Species A_16S_gut",]$modularity = modularity(modules_16S_odoA_gut)

refit_16S_odoB_gut = getRefit(mb_res_16S_odoB_gut)
graph_16S_odoB_gut = graph.adjacency(refit_16S_odoB_gut, mode = "undirected", weighted = TRUE)
modules_16S_odoB_gut = cluster_fast_greedy(graph_16S_odoB_gut)
modularity_df_gut_16S[modularity_df_gut_16S$species == "Species B_16S_gut",]$modularity = modularity(modules_16S_odoB_gut)

# Pseudocanthotermes
refit_16S_mili_gut = getRefit(mb_res_16S_mili_gut)
graph_16S_mili_gut = graph.adjacency(refit_16S_mili_gut, mode = "undirected", weighted = TRUE)
modules_16S_mili_gut = cluster_fast_greedy(graph_16S_mili_gut)
modularity_df_gut_16S[modularity_df_gut_16S$species == "militaris_16S_gut",]$modularity = modularity(modules_16S_mili_gut)

  # Comb #

modularity_df_comb_16S = data.frame(species = c("bellicosus_16S_comb", "subhyalinus_16S_comb", "cavithorax_16S_comb", "guinensis_16S_comb", "Species A_16S_comb",
                                                "Species B_16S_comb", "militaris_16S_comb"),
                                   modularity = NA)
# Macrotermes
refit_16S_belli_comb = getRefit(mb_res_16S_belli_comb)
graph_16S_belli_comb = graph.adjacency(refit_16S_belli_comb, mode = "undirected", weighted = TRUE)
modules_16S_belli_comb = cluster_fast_greedy(graph_16S_belli_comb)
modularity_df_comb_16S[modularity_df_comb_16S$species == "bellicosus_16S_comb",]$modularity = modularity(modules_16S_belli_comb)

refit_16S_sub_comb = getRefit(mb_res_16S_sub_comb)
graph_16S_sub_comb = graph.adjacency(refit_16S_sub_comb, mode = "undirected", weighted = TRUE)
modules_16S_sub_comb = cluster_fast_greedy(graph_16S_sub_comb)
modularity_df_comb_16S[modularity_df_comb_16S$species == "subhyalinus_16S_comb",]$modularity = modularity(modules_16S_sub_comb)

# Ancistrotermes
refit_16S_cavi_comb = getRefit(mb_res_16S_cavi_comb)
graph_16S_cavi_comb = graph.adjacency(refit_16S_cavi_comb, mode = "undirected", weighted = TRUE)
modules_16S_cavi_comb = cluster_fast_greedy(graph_16S_cavi_comb)
modularity_df_comb_16S[modularity_df_comb_16S$species == "cavithorax_16S_comb",]$modularity = modularity(modules_16S_cavi_comb)

refit_16S_guin_comb = getRefit(mb_res_16S_guin_comb)
graph_16S_guin_comb = graph.adjacency(refit_16S_guin_comb, mode = "undirected", weighted = TRUE)
modules_16S_guin_comb = cluster_fast_greedy(graph_16S_guin_comb)
modularity_df_comb_16S[modularity_df_comb_16S$species == "guinensis_16S_comb",]$modularity = modularity(modules_16S_guin_comb)

# Odontotermes
refit_16S_odoA_comb = getRefit(mb_res_16S_odoA_comb)
graph_16S_odoA_comb = graph.adjacency(refit_16S_odoA_comb, mode = "undirected", weighted = TRUE)
modules_16S_odoA_comb = cluster_fast_greedy(graph_16S_odoA_comb)
modularity_df_comb_16S[modularity_df_comb_16S$species == "Species A_16S_comb",]$modularity = modularity(modules_16S_odoA_comb)

refit_16S_odoB_comb = getRefit(mb_res_16S_odoB_comb)
graph_16S_odoB_comb = graph.adjacency(refit_16S_odoB_comb, mode = "undirected", weighted = TRUE)
modules_16S_odoB_comb = cluster_fast_greedy(graph_16S_odoB_comb)
modularity_df_comb_16S[modularity_df_comb_16S$species == "Species B_16S_comb",]$modularity = modularity(modules_16S_odoB_comb)

# Pseudocanthotermes
refit_16S_mili_comb = getRefit(mb_res_16S_mili_comb)
graph_16S_mili_comb = graph.adjacency(refit_16S_mili_comb, mode = "undirected", weighted = TRUE)
modules_16S_mili_comb = cluster_fast_greedy(graph_16S_mili_comb)
modularity_df_comb_16S[modularity_df_comb_16S$species == "militaris_16S_comb",]$modularity = modularity(modules_16S_mili_comb)


## AD ##
#======#

# Gut #
modularity_df_gut_AD = data.frame(species = c("bellicosus_AD_gut", "subhyalinus_AD_gut", "cavithorax_AD_gut", "guinensis_AD_gut", "Species A_AD_gut", 
                                              "Species B_AD_gut", "militaris_AD_gut"),
                                  modularity = NA)

# Macrotermes
refit_AD_belli_gut = getRefit(mb_res_AD_belli_gut)
graph_AD_belli_gut = graph.adjacency(refit_AD_belli_gut, mode = "undirected", weighted = TRUE)
modules_AD_belli_gut = cluster_fast_greedy(graph_AD_belli_gut)
modularity_df_gut_AD[modularity_df_gut_AD$species == "bellicosus_AD_gut",]$modularity = modularity(modules_AD_belli_gut)

refit_AD_sub_gut = getRefit(mb_res_AD_sub_gut)
graph_AD_sub_gut = graph.adjacency(refit_AD_sub_gut, mode = "undirected", weighted = TRUE)
modules_AD_sub_gut = cluster_fast_greedy(graph_AD_sub_gut)
modularity_df_gut_AD[modularity_df_gut_AD$species == "subhyalinus_AD_gut",]$modularity = modularity(modules_AD_sub_gut)

# Ancistrotermes
refit_AD_cavi_gut = getRefit(mb_res_AD_cavi_gut)
graph_AD_cavi_gut = graph.adjacency(refit_AD_cavi_gut, mode = "undirected", weighted = TRUE)
modules_AD_cavi_gut = cluster_fast_greedy(graph_AD_cavi_gut)
modularity_df_gut_AD[modularity_df_gut_AD$species == "cavithorax_AD_gut",]$modularity = modularity(modules_AD_cavi_gut)

refit_AD_guin_gut = getRefit(mb_res_AD_guin_gut)
graph_AD_guin_gut = graph.adjacency(refit_AD_guin_gut, mode = "undirected", weighted = TRUE)
modules_AD_guin_gut = cluster_fast_greedy(graph_AD_guin_gut)
modularity_df_gut_AD[modularity_df_gut_AD$species == "guinensis_AD_gut",]$modularity = modularity(modules_AD_guin_gut)

# Odontotermes
refit_AD_odoA_gut = getRefit(mb_res_AD_odoA_gut)
graph_AD_odoA_gut = graph.adjacency(refit_AD_odoA_gut, mode = "undirected", weighted = TRUE)
modules_AD_odoA_gut = cluster_fast_greedy(graph_AD_odoA_gut)
modularity_df_gut_AD[modularity_df_gut_AD$species == "Species A_AD_gut",]$modularity = modularity(modules_AD_odoA_gut)

refit_AD_odoB_gut = getRefit(mb_res_AD_odoB_gut)
graph_AD_odoB_gut = graph.adjacency(refit_AD_odoB_gut, mode = "undirected", weighted = TRUE)
modules_AD_odoB_gut = cluster_fast_greedy(graph_AD_odoB_gut)
modularity_df_gut_AD[modularity_df_gut_AD$species == "Species B_AD_gut",]$modularity = modularity(modules_AD_odoB_gut)

# Pseudocanthotermes
refit_AD_mili_gut = getRefit(mb_res_AD_mili_gut)
graph_AD_mili_gut = graph.adjacency(refit_AD_mili_gut, mode = "undirected", weighted = TRUE)
modules_AD_mili_gut = cluster_fast_greedy(graph_AD_mili_gut)
modularity_df_gut_AD[modularity_df_gut_AD$species == "militaris_AD_gut",]$modularity = modularity(modules_AD_mili_gut)

# Comb #

modularity_df_comb_AD = data.frame(species = c("bellicosus_AD_comb", "subhyalinus_AD_comb", "cavithorax_AD_comb", "guinensis_AD_comb", "Species A_AD_comb", 
                                               "Species B_AD_comb", "militaris_AD_comb"),
                                    modularity = NA)
# Macrotermes
refit_AD_belli_comb = getRefit(mb_res_AD_belli_comb)
graph_AD_belli_comb = graph.adjacency(refit_AD_belli_comb, mode = "undirected", weighted = TRUE)
modules_AD_belli_comb = cluster_fast_greedy(graph_AD_belli_comb)
modularity_df_comb_AD[modularity_df_comb_AD$species == "bellicosus_AD_comb",]$modularity = modularity(modules_AD_belli_comb)

refit_AD_sub_comb = getRefit(mb_res_AD_sub_comb)
graph_AD_sub_comb = graph.adjacency(refit_AD_sub_comb, mode = "undirected", weighted = TRUE)
modules_AD_sub_comb = cluster_fast_greedy(graph_AD_sub_comb)
modularity_df_comb_AD[modularity_df_comb_AD$species == "subhyalinus_AD_comb",]$modularity = modularity(modules_AD_sub_comb)

# Ancistrotermes
refit_AD_cavi_comb = getRefit(mb_res_AD_cavi_comb)
graph_AD_cavi_comb = graph.adjacency(refit_AD_cavi_comb, mode = "undirected", weighted = TRUE)
modules_AD_cavi_comb = cluster_fast_greedy(graph_AD_cavi_comb)
modularity_df_comb_AD[modularity_df_comb_AD$species == "cavithorax_AD_comb",]$modularity = modularity(modules_AD_cavi_comb)

refit_AD_guin_comb = getRefit(mb_res_AD_guin_comb)
graph_AD_guin_comb = graph.adjacency(refit_AD_guin_comb, mode = "undirected", weighted = TRUE)
modules_AD_guin_comb = cluster_fast_greedy(graph_AD_guin_comb)
modularity_df_comb_AD[modularity_df_comb_AD$species == "guinensis_AD_comb",]$modularity = modularity(modules_AD_guin_comb)

# Odontotermes
refit_AD_odoA_comb = getRefit(mb_res_AD_odoA_comb)
graph_AD_odoA_comb = graph.adjacency(refit_AD_odoA_comb, mode = "undirected", weighted = TRUE)
modules_AD_odoA_comb = cluster_fast_greedy(graph_AD_odoA_comb)
modularity_df_comb_AD[modularity_df_comb_AD$species == "Species A_AD_comb",]$modularity = modularity(modules_AD_odoA_comb)

refit_AD_odoB_comb = getRefit(mb_res_AD_odoB_comb)
graph_AD_odoB_comb = graph.adjacency(refit_AD_odoB_comb, mode = "undirected", weighted = TRUE)
modules_AD_odoB_comb = cluster_fast_greedy(graph_AD_odoB_comb)
modularity_df_comb_AD[modularity_df_comb_AD$species == "Species B_AD_comb",]$modularity = modularity(modules_AD_odoB_comb)

# Pseudocanthotermes
refit_AD_mili_comb = getRefit(mb_res_AD_mili_comb)
graph_AD_mili_comb = graph.adjacency(refit_AD_mili_comb, mode = "undirected", weighted = TRUE)
modules_AD_mili_comb = cluster_fast_greedy(graph_AD_mili_comb)
modularity_df_comb_AD[modularity_df_comb_AD$species == "militaris_AD_comb",]$modularity = modularity(modules_AD_mili_comb)

### Modularity stats ###

t.test(modularity_df_comb_16S$modularity, modularity_df_gut_16S$modularity)
sd(modularity_df_comb_16S$modularity)
sd(modularity_df_gut_16S$modularity)

t.test(modularity_df_comb_AD$modularity, modularity_df_gut_AD$modularity)
sd(modularity_df_comb_AD$modularity)
sd(modularity_df_gut_AD$modularity)

t.test(modularity_df_comb_16S$modularity, modularity_df_comb_AD$modularity)
t.test(modularity_df_gut_16S$modularity, modularity_df_gut_AD$modularity)


### Plot modularity ###

# Gut #
mod_cols = c("bellicosus_16S_comb" = "#E77675", "subhyalinus_16S_comb" = "#A31E23", "cavithorax_16S_comb" = "#E3A86C", "guinensis_16S_comb" = "#BA7328",
             "Species A_16S_comb" = "#91CB7B", "Species B_16S_comb" = "#3B7538", "militaris_16S_comb" = "#8499C7", 
             "bellicosus_16S_gut" = "#E77675", "subhyalinus_16S_gut" = "#A31E23", "cavithorax_16S_gut" = "#E3A86C", "guinensis_16S_gut" = "#BA7328",
             "Species A_16S_gut" = "#91CB7B","Species B_16S_gut" = "#3B7538", "militaris_16S_gut" = "#8499C7")

modularity_df_16S = rbind(modularity_df_gut_16S, modularity_df_comb_16S)
modularity_df_16S$species = factor(modularity_df_16S$species, levels = c("bellicosus_16S_comb", "subhyalinus_16S_comb", "cavithorax_16S_comb", "guinensis_16S_comb",
                                                                         "Species A_16S_comb", "Species B_16S_comb", "militaris_16S_comb", 
                                                                         "militaris_16S_gut", "Species B_16S_gut", "Species A_16S_gut", "guinensis_16S_gut", 
                                                                         "cavithorax_16S_gut", "subhyalinus_16S_gut", "bellicosus_16S_gut"))
label_data = modularity_df_16S %>%
  slice(match(c("bellicosus_16S_comb", "subhyalinus_16S_comb", "cavithorax_16S_comb", "guinensis_16S_comb",
                "Species A_16S_comb", "Species B_16S_comb", "militaris_16S_comb", 
                "militaris_16S_gut", "Species B_16S_gut", "Species A_16S_gut", "guinensis_16S_gut", 
                "cavithorax_16S_gut", "subhyalinus_16S_gut", "bellicosus_16S_gut"), species)) %>%
  rownames_to_column(var = "id")
# calculate the ANGLE of the labels
number_of_bar = nrow(label_data)
angle =  90 - 360 * (as.numeric(label_data$id)-0.5) /number_of_bar
label_data$hjust = ifelse(angle < -90, 1, 0)
label_data$angle = ifelse(angle < -90, angle+180, angle)
                         
mod_plot_16S = ggplot(modularity_df_16S, aes(y = modularity, x = species, fill = species))+
  geom_bar(stat="identity")+
  geom_text(data = label_data, aes(x = as.numeric(id), y = modularity+1, label = round(modularity, 3), hjust = hjust), 
            angle = label_data$angle, fontface="bold", size = 4, inherit.aes = FALSE)+
  coord_polar(start = 0)+
  ylim(-10,2)+
  scale_fill_manual(values = mod_cols)+
  theme_pubclean()+
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank())

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/mod_plot_16S.pdf", mod_plot_16S)

# Comb #

mod_cols = c("bellicosus_AD_comb" = "#E77675", "subhyalinus_AD_comb" = "#A31E23", "cavithorax_AD_comb" = "#E3A86C", "guinensis_AD_comb" = "#BA7328",
             "Species A_AD_comb" = "#91CB7B", "Species B_AD_comb" = "#3B7538", "militaris_AD_comb" = "#8499C7", 
             "bellicosus_AD_gut" = "#E77675", "subhyalinus_AD_gut" = "#A31E23", "cavithorax_AD_gut" = "#E3A86C", "guinensis_AD_gut" = "#BA7328",
             "Species A_AD_gut" = "#91CB7B","Species B_AD_gut" = "#3B7538", "militaris_AD_gut" = "#8499C7")

modularity_df_AD = rbind(modularity_df_gut_AD, modularity_df_comb_AD)
modularity_df_AD$species = factor(modularity_df_AD$species, levels = c("bellicosus_AD_comb", "subhyalinus_AD_comb", "cavithorax_AD_comb", "guinensis_AD_comb",
                                                                         "Species A_AD_comb", "Species B_AD_comb", "militaris_AD_comb", 
                                                                         "militaris_AD_gut", "Species B_AD_gut", "Species A_AD_gut", "guinensis_AD_gut", 
                                                                         "cavithorax_AD_gut", "subhyalinus_AD_gut", "bellicosus_AD_gut"))
label_data = modularity_df_AD %>%
  slice(match(c("bellicosus_AD_comb", "subhyalinus_AD_comb", "cavithorax_AD_comb", "guinensis_AD_comb",
                "Species A_AD_comb", "Species B_AD_comb", "militaris_AD_comb", 
                "militaris_AD_gut", "Species B_AD_gut", "Species A_AD_gut", "guinensis_AD_gut", 
                "cavithorax_AD_gut", "subhyalinus_AD_gut", "bellicosus_AD_gut"), species)) %>%
  rownames_to_column(var = "id")
# calculate the ANGLE of the labels
number_of_bar = nrow(label_data)
angle =  90 - 360 * (as.numeric(label_data$id)-0.5) /number_of_bar
label_data$hjust = ifelse(angle < -90, 1, 0)
label_data$angle = ifelse(angle < -90, angle+180, angle)

mod_plot_AD = ggplot(modularity_df_AD, aes(y = modularity, x = species, fill = species))+
  geom_bar(stat="identity")+
  geom_text(data = label_data, aes(x = as.numeric(id), y = modularity+1, label = round(modularity, 3), hjust = hjust), 
            angle = label_data$angle, fontface="bold", size = 4, inherit.aes = FALSE)+
  coord_polar(start = 0)+
  ylim(-10,2)+
  scale_fill_manual(values = mod_cols)+
  theme_pubclean()+
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank())


ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 3/mod_plot_AD.pdf", mod_plot_AD)


#### Modularity end ####

#### Metrics end ####


#### Plot networks ####
#=====================#

filt_vert = function(g){
  degree_vec = igraph::degree(g, mode = "total")
  vertices = sort(degree_vec,decreasing = T)[1:200]
  return(min(vertices))
} 


degree_cols = c("0" = "grey", "1" = "#D73027", "2" = "#D73027", "3" = "#D73027", "4" = "#D73027", "5" = "#D73027", "6" = "#D73027", "7" = "#D73027", 
                "8" = "#D73027","9" = "#D73027", "10" = "#D73027",
                "11" = "#FC8D59","12" = "#FC8D59", "13" = "#FC8D59", "14" = "#FC8D59", "15" = "#FC8D59", "16" = "#FC8D59","17" = "#FC8D59", 
                "18" = "#FC8D59", "19" = "#FC8D59", "20" = "#FC8D59",
                "21" = "#FEE090", "22" = "#FEE090", "23" = "#FEE090", "24" = "#FEE090","25" = "#FEE090", "26" = "#FEE090", "27" = "#FEE090", 
                "28" = "#FEE090", "29" = "#FEE090", "30" = "#FEE090",
                "31" = "#E0F3F8", "32" = "#E0F3F8", "33" = "#E0F3F8", "34" = "#E0F3F8", "35" = "#E0F3F8", "36" = "#E0F3F8", "37" = "#E0F3F8" ,
                "38" = "#E0F3F8", "39" = "#E0F3F8", "40" = "#E0F3F8",
                "41" = "#abc4ed", "42" = "#abc4ed", "43" = "#abc4ed", "44" = "#abc4ed", "45" = "#abc4ed" ,"46" = "#abc4ed" ,"47" = "#abc4ed", 
                "48" = "#abc4ed", "49" = "#abc4ed", "50" = "#abc4ed", 
                "51" = "#4575B4", "52" = "#4575B4", "53" = "#4575B4" ,"54" = "#4575B4" ,"55" = "#4575B4", "56" = "#4575B4","57" = "#4575B4", 
                "58" = "#4575B4", "59" = "#4575B4", "60" = "#4575B4", 
                "61" = "#024ec9" ,"62" = "#024ec9", "63" = "#024ec9", "64" = "#024ec9", "65" = "#024ec9", "66" = "#024ec9", "67" = "#024ec9",
                "68" = "#024ec9", "69" = "#024ec9", "70" = "#024ec9")
  ## 16S ##

# Gut #

# M. bellicosus

n_16S_belli_gut = symBeta(getOptBeta(mb_res_16S_belli_gut))
graph_16S_belli_gut = graph.adjacency(n_16S_belli_gut, mode = "undirected", weighted = TRUE)

graph_16S_belli_gut_sub = delete_edges(graph_16S_belli_gut, E(graph_16S_belli_gut)[weight >= max(weight)*.50])

degree_vec = igraph::degree(graph_16S_belli_gut, mode = "total")

min_degree = filt_vert(graph_16S_belli_gut_sub)

E(graph_16S_belli_gut)$color = ifelse(E(graph_16S_belli_gut)$weight  > 0, "#85bc38", "#e2523c")

plot_16S_belli_gut = ggnet2(graph_16S_belli_gut_sub, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

gsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_belli_gut.pdf", plot_16S_belli_gut)


# M. subhyalinus

n_16S_sub_gut = symBeta(getOptBeta(mb_res_16S_sub_gut))
graph_16S_sub_gut = graph.adjacency(n_16S_sub_gut, mode = "undirected", weighted = TRUE)

degree_vec = igraph::degree(graph_16S_sub_gut, mode = "total")

min_degree = filt_vert(graph_16S_sub_gut)

E(graph_16S_sub_gut)$color = ifelse(E(graph_16S_sub_gut)$weight > 0, "#85bc38", "#e2523c")

plot_16S_sub_gut = ggnet2(graph_16S_sub_gut, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_sub_gut.pdf", plot_16S_sub_gut)


# A. cavi
n_16S_cavi_gut = symBeta(getOptBeta(mb_res_16S_cavi_gut))
graph_16S_cavi_gut = graph.adjacency(n_16S_cavi_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_cavi_gut, mode = "total")

min_degree = filt_vert(graph_16S_cavi_gut)

E(graph_16S_cavi_gut)$color = ifelse(E(graph_16S_cavi_gut)$weight > 0, "#85bc38", "#e2523c")

plot_16S_cavi_gut = ggnet2(graph_16S_cavi_gut, size = "degree", color = degree_vec, palette = degree_cols,
                          legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                          edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_cavi_gut.pdf", plot_16S_cavi_gut)


# A. guin

n_16S_guin_gut = symBeta(getOptBeta(mb_res_16S_guin_gut))
graph_16S_guin_gut = graph.adjacency(n_16S_guin_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_guin_gut, mode = "total")

min_degree = filt_vert(graph_16S_guin_gut)

E(graph_16S_guin_gut)$color = ifelse(E(graph_16S_guin_gut)$weight > 0, "#85bc38", "#e2523c")


plot_16S_guin_gut = ggnet2(graph_16S_guin_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_guin_gut.pdf", plot_16S_guin_gut)

# O. species A

n_16S_odoA_gut = symBeta(getOptBeta(mb_res_16S_odoA_gut))
graph_16S_odoA_gut = graph.adjacency(n_16S_odoA_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_odoA_gut, mode = "total")

min_degree = filt_vert(graph_16S_odoA_gut)

E(graph_16S_odoA_gut)$color = ifelse(E(graph_16S_odoA_gut)$weight > 0, "#85bc38", "#e2523c")


plot_16S_odoA_gut = ggnet2(graph_16S_odoA_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_odoA_gut.pdf", plot_16S_odoA_gut)

# O. species B

n_16S_odoB_gut = symBeta(getOptBeta(mb_res_16S_odoB_gut))
graph_16S_odoB_gut = graph.adjacency(n_16S_odoB_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_odoB_gut, mode = "total")

min_degree = filt_vert(graph_16S_odoB_gut)

E(graph_16S_odoB_gut)$color = ifelse(E(graph_16S_odoB_gut)$weight > 0, "#85bc38", "#e2523c")


plot_16S_odoB_gut = ggnet2(graph_16S_odoB_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_odoB_gut.pdf", plot_16S_odoB_gut)

# P. mili

n_16S_mili_gut = symBeta(getOptBeta(mb_res_16S_mili_gut))
graph_16S_mili_gut = graph.adjacency(n_16S_mili_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_mili_gut, mode = "total")

min_degree = filt_vert(graph_16S_mili_gut)

E(graph_16S_mili_gut)$color = ifelse(E(graph_16S_mili_gut)$weight > 0, "#85bc38", "#e2523c")


plot_16S_mili_gut = ggnet2(graph_16S_mili_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_mili_gut.pdf", plot_16S_mili_gut)



  # Comb #

# M. bellicosus

n_16S_belli_comb = symBeta(getOptBeta(mb_res_16S_belli_comb))
graph_16S_belli_comb = graph.adjacency(n_16S_belli_comb, mode = "undirected", weighted = TRUE)


degree_vec = igraph::degree(graph_16S_belli_comb, mode = "total")

min_degree = filt_vert(graph_16S_belli_comb)

E(graph_16S_belli_comb)$color = ifelse(E(graph_16S_belli_comb)$weight  > 0, "#85bc38", "#e2523c")

plot_16S_belli_comb = ggnet2(graph_16S_belli_comb, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_belli_comb.pdf", plot_16S_belli_comb)


# M. subhyalinus

n_16S_sub_comb = symBeta(getOptBeta(mb_res_16S_sub_comb))
graph_16S_sub_comb = graph.adjacency(n_16S_sub_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_sub_comb, mode = "total")

min_degree = filt_vert(graph_16S_sub_comb)

E(graph_16S_sub_comb)$color = ifelse(E(graph_16S_sub_comb)$weight > 0, "#85bc38", "#e2523c")

plot_16S_sub_comb = ggnet2(graph_16S_sub_comb, size = "degree", color = degree_vec, palette = degree_cols,
                          legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                          edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_sub_comb.pdf", plot_16S_sub_comb)


# A. cavi
n_16S_cavi_comb = symBeta(getOptBeta(mb_res_16S_cavi_comb))
graph_16S_cavi_comb = graph.adjacency(n_16S_cavi_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_cavi_comb, mode = "total")

min_degree = filt_vert(graph_16S_cavi_comb)

E(graph_16S_cavi_comb)$color = ifelse(E(graph_16S_cavi_comb)$weight > 0, "#85bc38", "#e2523c")

plot_16S_cavi_comb = ggnet2(graph_16S_cavi_comb, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_cavi_comb.pdf", plot_16S_cavi_comb)


# A. guin

n_16S_guin_comb = symBeta(getOptBeta(mb_res_16S_guin_comb))
graph_16S_guin_comb = graph.adjacency(n_16S_guin_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_guin_comb, mode = "total")

min_degree = filt_vert(graph_16S_guin_comb)

E(graph_16S_guin_comb)$color = ifelse(E(graph_16S_guin_comb)$weight > 0, "#85bc38", "#e2523c")


plot_16S_guin_comb = ggnet2(graph_16S_guin_comb, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_guin_comb.pdf", plot_16S_guin_comb)

# O. species A

n_16S_odoA_comb = symBeta(getOptBeta(mb_res_16S_odoA_comb))
graph_16S_odoA_comb = graph.adjacency(n_16S_odoA_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_odoA_comb, mode = "total")

min_degree = filt_vert(graph_16S_odoA_comb)

E(graph_16S_odoA_comb)$color = ifelse(E(graph_16S_odoA_comb)$weight > 0, "#85bc38", "#e2523c")


plot_16S_odoA_comb = ggnet2(graph_16S_odoA_comb, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_odoA_comb.pdf", plot_16S_odoA_comb)

# O. species B

n_16S_odoB_comb = symBeta(getOptBeta(mb_res_16S_odoB_comb))
graph_16S_odoB_comb = graph.adjacency(n_16S_odoB_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_odoB_comb, mode = "total")

min_degree = filt_vert(graph_16S_odoB_comb)

E(graph_16S_odoB_comb)$color = ifelse(E(graph_16S_odoB_comb)$weight > 0, "#85bc38", "#e2523c")


plot_16S_odoB_comb = ggnet2(graph_16S_odoB_comb, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_odoB_comb.pdf", plot_16S_odoB_comb)

# P. mili

n_16S_mili_comb = symBeta(getOptBeta(mb_res_16S_mili_comb))
graph_16S_mili_comb = graph.adjacency(n_16S_mili_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_16S_mili_comb, mode = "total")

min_degree = filt_vert(graph_16S_mili_comb)

E(graph_16S_mili_comb)$color = ifelse(E(graph_16S_mili_comb)$weight > 0, "#85bc38", "#e2523c")


plot_16S_mili_comb = ggnet2(graph_16S_mili_comb, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_16S_mili_comb.pdf", plot_16S_mili_comb)




  ## AD ##
  #======#

# Gut #

# M. bellicosus

n_AD_belli_gut = symBeta(getOptBeta(mb_res_AD_belli_gut))
graph_AD_belli_gut = graph.adjacency(n_AD_belli_gut, mode = "undirected", weighted = TRUE)


degree_vec = igraph::degree(graph_AD_belli_gut, mode = "total")

min_degree = filt_vert(graph_AD_belli_gut)

E(graph_AD_belli_gut)$color = ifelse(E(graph_AD_belli_gut)$weight  > 0, "#85bc38", "#e2523c")

plot_AD_belli_gut = ggnet2(graph_AD_belli_gut, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_belli_gut_alt.pdf", plot_AD_belli_gut)


# M. subhyalinus

n_AD_sub_gut = symBeta(getOptBeta(mb_res_AD_sub_gut))
graph_AD_sub_gut = graph.adjacency(n_AD_sub_gut, mode = "undirected", weighted = TRUE)

degree_vec = igraph::degree(graph_AD_sub_gut, mode = "total")

min_degree = filt_vert(graph_AD_sub_gut)

E(graph_AD_sub_gut)$color = ifelse(E(graph_AD_sub_gut)$weight > 0, "#85bc38", "#e2523c")

plot_AD_sub_gut = ggnet2(graph_AD_sub_gut, size = "degree", color = degree_vec, palette = degree_cols,
                          legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                          edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_sub_gut.pdf", plot_AD_sub_gut)


# A. cavi
n_AD_cavi_gut = symBeta(getOptBeta(mb_res_AD_cavi_gut))
graph_AD_cavi_gut = graph.adjacency(n_AD_cavi_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_cavi_gut, mode = "total")

min_degree = filt_vert(graph_AD_cavi_gut)

E(graph_AD_cavi_gut)$color = ifelse(E(graph_AD_cavi_gut)$weight > 0, "#85bc38", "#e2523c")

plot_AD_cavi_gut = ggnet2(graph_AD_cavi_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_cavi_gut.pdf", plot_AD_cavi_gut)


# A. guin

n_AD_guin_gut = symBeta(getOptBeta(mb_res_AD_guin_gut))
graph_AD_guin_gut = graph.adjacency(n_AD_guin_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_guin_gut, mode = "total")

min_degree = filt_vert(graph_AD_guin_gut)

E(graph_AD_guin_gut)$color = ifelse(E(graph_AD_guin_gut)$weight > 0, "#85bc38", "#e2523c")


plot_AD_guin_gut = ggnet2(graph_AD_guin_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_guin_gut.pdf", plot_AD_guin_gut)

# O. species A

n_AD_odoA_gut = symBeta(getOptBeta(mb_res_AD_odoA_gut))
graph_AD_odoA_gut = graph.adjacency(n_AD_odoA_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_odoA_gut, mode = "total")

min_degree = filt_vert(graph_AD_odoA_gut)

E(graph_AD_odoA_gut)$color = ifelse(E(graph_AD_odoA_gut)$weight > 0, "#85bc38", "#e2523c")


plot_AD_odoA_gut = ggnet2(graph_AD_odoA_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_odoA_gut.pdf", plot_AD_odoA_gut)

# O. species B

n_AD_odoB_gut = symBeta(getOptBeta(mb_res_AD_odoB_gut))
graph_AD_odoB_gut = graph.adjacency(n_AD_odoB_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_odoB_gut, mode = "total")

min_degree = filt_vert(graph_AD_odoB_gut)

E(graph_AD_odoB_gut)$color = ifelse(E(graph_AD_odoB_gut)$weight > 0, "#85bc38", "#e2523c")


plot_AD_odoB_gut = ggnet2(graph_AD_odoB_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_odoB_gut.pdf", plot_AD_odoB_gut)

# P. mili

n_AD_mili_gut = symBeta(getOptBeta(mb_res_AD_mili_gut))
graph_AD_mili_gut = graph.adjacency(n_AD_mili_gut, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_mili_gut, mode = "total")

min_degree = filt_vert(graph_AD_mili_gut)

E(graph_AD_mili_gut)$color = ifelse(E(graph_AD_mili_gut)$weight > 0, "#85bc38", "#e2523c")


plot_AD_mili_gut = ggnet2(graph_AD_mili_gut, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_mili_gut.pdf", plot_AD_mili_gut)



# Comb #

# M. bellicosus

n_AD_belli_comb = symBeta(getOptBeta(mb_res_AD_belli_comb))
graph_AD_belli_comb = graph.adjacency(n_AD_belli_comb, mode = "undirected", weighted = TRUE)


degree_vec = igraph::degree(graph_AD_belli_comb, mode = "total")

min_degree = filt_vert(graph_AD_belli_comb)

E(graph_AD_belli_comb)$color = ifelse(E(graph_AD_belli_comb)$weight  > 0, "#85bc38", "#e2523c")

plot_AD_belli_comb = ggnet2(graph_AD_belli_comb, size = "degree", color = degree_vec, palette = degree_cols,
                             legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                             edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_belli_comb.pdf", plot_AD_belli_comb)


# M. subhyalinus

n_AD_sub_comb = symBeta(getOptBeta(mb_res_AD_sub_comb))
graph_AD_sub_comb = graph.adjacency(n_AD_sub_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_sub_comb, mode = "total")

min_degree = filt_vert(graph_AD_sub_comb)

E(graph_AD_sub_comb)$color = ifelse(E(graph_AD_sub_comb)$weight > 0, "#85bc38", "#e2523c")

plot_AD_sub_comb = ggnet2(graph_AD_sub_comb, size = "degree", color = degree_vec, palette = degree_cols,
                           legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                           edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_sub_comb.pdf", plot_AD_sub_comb)


# A. cavi
n_AD_cavi_comb = symBeta(getOptBeta(mb_res_AD_cavi_comb))
graph_AD_cavi_comb = graph.adjacency(n_AD_cavi_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_cavi_comb, mode = "total")

min_degree = filt_vert(graph_AD_cavi_comb)

E(graph_AD_cavi_comb)$color = ifelse(E(graph_AD_cavi_comb)$weight > 0, "#85bc38", "#e2523c")

plot_AD_cavi_comb = ggnet2(graph_AD_cavi_comb, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_cavi_comb.pdf", plot_AD_cavi_comb)


# A. guin

n_AD_guin_comb = symBeta(getOptBeta(mb_res_AD_guin_comb))
graph_AD_guin_comb = graph.adjacency(n_AD_guin_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_guin_comb, mode = "total")

min_degree = filt_vert(graph_AD_guin_comb)

E(graph_AD_guin_comb)$color = ifelse(E(graph_AD_guin_comb)$weight > 0, "#85bc38", "#e2523c")


plot_AD_guin_comb = ggnet2(graph_AD_guin_comb, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_guin_comb.pdf", plot_AD_guin_comb)

# O. species A

n_AD_odoA_comb = symBeta(getOptBeta(mb_res_AD_odoA_comb))
graph_AD_odoA_comb = graph.adjacency(n_AD_odoA_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_odoA_comb, mode = "total")

min_degree = filt_vert(graph_AD_odoA_comb)

E(graph_AD_odoA_comb)$color = ifelse(E(graph_AD_odoA_comb)$weight > 0, "#85bc38", "#e2523c")


plot_AD_odoA_comb = ggnet2(graph_AD_odoA_comb, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_odoA_comb.pdf", plot_AD_odoA_comb)

# O. species B

n_AD_odoB_comb = symBeta(getOptBeta(mb_res_AD_odoB_comb))
graph_AD_odoB_comb = graph.adjacency(n_AD_odoB_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_odoB_comb, mode = "total")

min_degree = filt_vert(graph_AD_odoB_comb)

E(graph_AD_odoB_comb)$color = ifelse(E(graph_AD_odoB_comb)$weight > 0, "#85bc38", "#e2523c")


plot_AD_odoB_comb = ggnet2(graph_AD_odoB_comb, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_odoB_comb.pdf", plot_AD_odoB_comb)

# P. mili

n_AD_mili_comb = symBeta(getOptBeta(mb_res_AD_mili_comb))
graph_AD_mili_comb = graph.adjacency(n_AD_mili_comb, mode = "undirected", weighted = TRUE)

degree_vec = degree(graph_AD_mili_comb, mode = "total")

min_degree = filt_vert(graph_AD_mili_comb)

E(graph_AD_mili_comb)$color = ifelse(E(graph_AD_mili_comb)$weight > 0, "#85bc38", "#e2523c")


plot_AD_mili_comb = ggnet2(graph_AD_mili_comb, size = "degree", color = degree_vec, palette = degree_cols,
                            legend.position = "none", mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75),
                            edge.alpha = 0.8, edge.color = "color", size.min = min_degree)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2_v2/plot_minimal_AD_mili_comb.pdf", plot_AD_mili_comb)

#### network plots end ####