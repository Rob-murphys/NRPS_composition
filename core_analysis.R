library(tibble)
library(dplyr)
library(vegan)
### Import Data ###
#=================#

# AD
feature_table_AD.norm = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_normalised_OBU_table-filtered_2.rds") 

metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  mutate(sample_type = case_when(caste_phase %in% c("YC", "MC", "OC", "C") ~ "Comb",
                                 caste_phase %in% c("SW", "LW", "SS", "LS", "S") ~ "Gut"))
## Removing Microtermes ##
metadata_noMicro = metadata %>%
  filter(genus != "Microtermes")

feature_table_AD.norm_noMicro = feature_table_AD.norm %>%
  filter(rownames(.) %in% metadata_noMicro$sample.id)

#### Get core gut AD ####
#=======================#
total_AD_core_gut = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW")) %>%
  select(!c("caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(total_AD_core_gut)

## Macrotermes ##
core_gut_AD_belli = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW"), species == "bellicosus") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_gut_AD_belli)

core_gut_AD_sub = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW"), species == "subhyalinus") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_gut_AD_sub)

## Ancistrotermes ##
core_gut_AD_cavi = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW"), species == "cavithorax") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_gut_AD_cavi)

core_gut_AD_guin = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW"), species == "guinensis") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_gut_AD_guin)

## Odontotermes ##

core_gut_AD_odoA = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW"), species == "Species A") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_gut_AD_odoA)

core_gut_AD_odoB = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW"), species == "Species B") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_gut_AD_odoB)

## Pseudocanthotermes ##
core_gut_AD_mili = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW"), species == "militaris-spiniger") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_gut_AD_mili)

### combined dataframes ###
## Macrotermes
gut_AD_belli = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("SW", "LW"), species == "bellicosus") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "bellicosus") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_belli), "core"))

gut_AD_sub = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("SW", "LW"), species == "subhyalinus") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "subhyalinus") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_sub), "core"))
## Pseudocanthotermes
gut_AD_mili = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("SW", "LW"), species == "militaris-spiniger") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "militaris-spiniger") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_mili), "core"))

## Odontotermers
gut_AD_odoA = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("SW", "LW"), species == "Species A") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "Species A") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_odoA), "core"))

gut_AD_odoB = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("SW", "LW"), species == "Species B") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "Species B") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_odoB), "core"))

## Ancistrotermes
gut_AD_cavi = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("SW", "LW"), species == "cavithorax") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "cavithorax") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_cavi), "core"))

gut_AD_guin = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("SW", "LW"), species == "guinensis") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "guinensis") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_guin), "core"))
##### core gut end ####

#### Get core comb AD ####
total_AD_core_comb = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C")) %>%
  select(!c("caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(total_AD_core_comb)

## Macrotermes ##
core_comb_AD_belli = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC"), species == "bellicosus") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_belli)

core_comb_AD_sub = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC"), species == "subhyalinus") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_sub)

## Ancistrotermes ##
core_comb_AD_cavi = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "cavithorax") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_cavi)

core_comb_AD_guin = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "guinensis") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_guin)

## Odontotermes ##

core_comb_AD_odoA = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "Species A") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_odoA)

core_comb_AD_odoB = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "Species B") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_odoB)

## Pseudocanthotermes ##
core_comb_AD_mili = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "militaris-spiniger") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_mili)

### combined dataframes ###

## Total
comb_AD_total = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C")) %>%
  select(!c("caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this species
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per species for OBUs present in that species
  mutate(core = "no") %>%
  mutate(core = replace(core, variable %in% colnames(total_AD_core_comb), "yes"))

## Macrotermes
comb_AD_belli = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "bellicosus") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "bellicosus") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_belli), "core"))

comb_AD_sub = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "subhyalinus") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "subhyalinus") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_sub), "core"))
## Pseudocanthotermes
comb_AD_mili = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "militaris-spiniger") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "militaris-spiniger") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_mili), "core"))

## Odontotermers
comb_AD_odoA = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "Species A") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "Species A") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_odoA), "core"))

comb_AD_odoB = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "Species B") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "Species B") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_odoB), "core"))

## Ancistrotermes
comb_AD_cavi = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "cavithorax") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "cavithorax") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_cavi), "core"))

comb_AD_guin = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "guinensis") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per species for OBUs present in that species
  mutate(species = "guinensis") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_guin), "core"))

df_AD_comb_species = rbind(comb_AD_belli, comb_AD_sub, comb_AD_mili, comb_AD_odoA, comb_AD_odoB, comb_AD_cavi, comb_AD_guin) %>% # rbindng all the above dataframes together
  mutate(colour = replace(colour, colour != "core", "not_core")) %>%
  select(variable, colour) %>%
  distinct(.)
##### core comb end ####

#### total relative abundance of core ####
#========================================#

## Gut ##
#=======#

# Macrotermes
(sum(gut_AD_belli[gut_AD_belli$colour == "core",]$mn)/sum(gut_AD_belli$mn))*100
(nrow(gut_AD_belli[gut_AD_belli$colour == "core",])/nrow(gut_AD_belli))*100

(sum(gut_AD_sub[gut_AD_sub$colour == "core",]$mn)/sum(gut_AD_sub$mn))*100
(nrow(gut_AD_sub[gut_AD_sub$colour == "core",])/nrow(gut_AD_sub))*100

# Pseudocanthotermes
(sum(gut_AD_mili[gut_AD_mili$colour == "core",]$mn)/sum(gut_AD_mili$mn))*100
(nrow(gut_AD_mili[gut_AD_mili$colour == "core",])/nrow(gut_AD_mili))*100

# Odontotermes
(sum(gut_AD_odoA[gut_AD_odoA$colour == "core",]$mn)/sum(gut_AD_odoA$mn))*100
(nrow(gut_AD_odoA[gut_AD_odoA$colour == "core",])/nrow(gut_AD_odoA))*100

(sum(gut_AD_odoB[gut_AD_odoB$colour == "core",]$mn)/sum(gut_AD_odoB$mn))*100
(nrow(gut_AD_odoB[gut_AD_odoB$colour == "core",])/nrow(gut_AD_odoB))*100

# Ancistrotermes
(sum(gut_AD_cavi[gut_AD_cavi$colour == "core",]$mn)/sum(gut_AD_cavi$mn))*100
(nrow(gut_AD_cavi[gut_AD_cavi$colour == "core",])/nrow(gut_AD_cavi))*100

(sum(gut_AD_guin[gut_AD_guin$colour == "core",]$mn)/sum(gut_AD_guin$mn))*100
(nrow(gut_AD_guin[gut_AD_guin$colour == "core",])/nrow(gut_AD_guin))*100

## Comb ##
#========#

# Total
(sum(comb_AD_total[comb_AD_total$core == "yes",]$mn)/sum(comb_AD_total$mn))*100
(nrow(comb_AD_total[comb_AD_total$core == "yes",])/nrow(comb_AD_total))*100

# Macrotermes
(sum(comb_AD_belli[comb_AD_belli$colour == "core",]$mn)/sum(comb_AD_belli$mn))*100
(nrow(comb_AD_belli[comb_AD_belli$colour == "core",])/nrow(comb_AD_belli))*100

(sum(comb_AD_sub[comb_AD_sub$colour == "core",]$mn)/sum(comb_AD_sub$mn))*100
(nrow(comb_AD_sub[comb_AD_sub$colour == "core",])/nrow(comb_AD_sub))*100

# Pseudocanthotermes
(sum(comb_AD_mili[comb_AD_mili$colour == "core",]$mn)/sum(comb_AD_mili$mn))*100
(nrow(comb_AD_mili[comb_AD_mili$colour == "core",])/nrow(comb_AD_mili))*100

# Odontotermes
(sum(comb_AD_odoA[comb_AD_odoA$colour == "core",]$mn)/sum(comb_AD_odoA$mn))*100
(nrow(comb_AD_odoA[comb_AD_odoA$colour == "core",])/nrow(comb_AD_odoA))*100

(sum(comb_AD_odoB[comb_AD_odoB$colour == "core",]$mn)/sum(comb_AD_odoB$mn))*100
(nrow(comb_AD_odoB[comb_AD_odoB$colour == "core",])/nrow(comb_AD_odoB))*100

# Ancistrotermes
(sum(comb_AD_cavi[comb_AD_cavi$colour == "core",]$mn)/sum(comb_AD_cavi$mn))*100
(nrow(comb_AD_cavi[comb_AD_cavi$colour == "core",])/nrow(comb_AD_cavi))*100

(sum(comb_AD_guin[comb_AD_guin$colour == "core",]$mn)/sum(comb_AD_guin$mn))*100
(nrow(comb_AD_guin[comb_AD_guin$colour == "core",])/nrow(comb_AD_guin))*100


#### total relative abundance of core end ####

#### Core comb from gut ####

gut2comb_df = data.frame(species = c("bellicosus", "subhyalinus", "cavithorax", "guinensis", 
                                     "Species A", "Species B", "militaris-spiniger"),
                         per_gut2comb = NA)
gut2comb_df$species = factor(gut2comb_df$species, levels = c("bellicosus", "subhyalinus", "cavithorax", "guinensis", 
                                                   "Species A", "Species B", "militaris-spiniger"))

# Macrotermes
gut2comb_df[gut2comb_df$species == "bellicosus",]$per_gut2comb = sum(colnames(core_comb_AD_belli) %in% colnames(core_gut_AD_belli))/ncol(core_comb_AD_belli)

gut2comb_df[gut2comb_df$species == "subhyalinus",]$per_gut2comb = sum(colnames(core_comb_AD_sub) %in% colnames(core_gut_AD_sub))/ncol(core_comb_AD_sub)

# Pseudocanthotermes
gut2comb_df[gut2comb_df$species == "militaris-spiniger",]$per_gut2comb = sum(colnames(core_comb_AD_mili) %in% colnames(core_gut_AD_mili))/ncol(core_comb_AD_mili)

# Odontotermes
gut2comb_df[gut2comb_df$species == "Species A",]$per_gut2comb = sum(colnames(core_comb_AD_odoA) %in% colnames(core_gut_AD_odoA))/ncol(core_comb_AD_odoA)

gut2comb_df[gut2comb_df$species == "Species B",]$per_gut2comb = sum(colnames(core_comb_AD_odoB) %in% colnames(core_gut_AD_odoB))/ncol(core_comb_AD_odoB)

# Ancistrotermes
gut2comb_df[gut2comb_df$species == "cavithorax",]$per_gut2comb = sum(colnames(core_comb_AD_cavi) %in% colnames(core_gut_AD_cavi))/ncol(core_comb_AD_cavi)

gut2comb_df[gut2comb_df$species == "guinensis",]$per_gut2comb = sum(colnames(core_comb_AD_guin) %in% colnames(core_gut_AD_guin))/ncol(core_comb_AD_guin)


#### Core comb from gut end ####


#### Phylogenetic comparison ####
#===============================#

dismat = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/dna-sequences.dismat", header = FALSE) %>%
  column_to_rownames(var = "V1")
colnames(dismat) = rownames(dismat)


# Bellicosus
dismat_filt_belli = dismat %>%
  filter(row.names(.) %in% gut_AD_belli$variable) %>%
  select(matches(as.character(gut_AD_belli$variable)))

phylo_res_belli = adonis2(dismat_filt_belli ~ colour, data = gut_AD_belli, by = "terms")

# Subhyalinus
dismat_filt_sub = dismat %>%
  filter(row.names(.) %in% gut_AD_sub$variable) %>%
  select(matches(as.character(gut_AD_belli$variable)))

phylo_res_sub = adonis2(dismat_filt_sub ~ colour, data = gut_AD_sub, by = "terms")


# Cavithorax
dismat_filt_cavi = dismat %>%
  filter(row.names(.) %in% gut_AD_cavi$variable) %>%
  select(matches(as.character(gut_AD_cavi$variable)))

phylo_res_cavi = adonis2(dismat_filt_cavi ~ colour, data = gut_AD_cavi, by = "terms")


# Guinensis
dismat_filt_guin = dismat %>%
  filter(row.names(.) %in% gut_AD_guin$variable) %>%
  select(matches(as.character(gut_AD_guin$variable)))

phylo_res_guin = adonis2(dismat_filt_guin ~ colour, data = gut_AD_guin, by = "terms")


# Species A
dismat_filt_odoA = dismat %>%
  filter(row.names(.) %in% gut_AD_odoA$variable) %>%
  select(matches(as.character(gut_AD_odoA$variable)))

phylo_res_odoA = adonis2(dismat_filt_odoA ~ colour, data = gut_AD_odoA, by = "terms")


# Species B
dismat_filt_odoB = dismat %>%
  filter(row.names(.) %in% gut_AD_odoB$variable) %>%
  select(matches(as.character(gut_AD_odoB$variable)))

phylo_res_odoB = adonis2(dismat_filt_odoB ~ colour, data = gut_AD_odoB, by = "terms")


# Militaris
dismat_filt_mili = dismat %>%
  filter(row.names(.) %in% gut_AD_mili$variable) %>%
  select(matches(as.character(gut_AD_mili$variable)))

phylo_res_mili = adonis2(dismat_filt_mili ~ colour, data = gut_AD_mili, by = "terms")


#### Phylogenetic comparison end ####