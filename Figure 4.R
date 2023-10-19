library(tibble)
library(dplyr)
library(ggpubr)
library(reshape2)
library(ggbeeswarm)

### Import Data ###
#=================#

# 16S
feature_table_16S.norm = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_normalised_feature_table.rds")

# AD
feature_table_AD.norm = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_normalised_OBU_table-filtered_2.rds") 

metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  mutate(sample_type = case_when(caste_phase %in% c("YC", "MC", "OC", "C") ~ "Comb",
                                 caste_phase %in% c("SW", "LW", "SS", "LS", "S") ~ "Gut"))

#### Get core gut AD ####
#=======================#

core_gut_AD = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("SW", "LW")) %>%
  select(!c("caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))

  ### By species ###

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

### Plotting gut cores ###
#========================#
# Here are are wanting to plot the core gut microbiomes in relation to every other OBU to show they are in general, highly abundant

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
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_belli), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_sub), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD), "total_core"))
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
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_mili), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_odoA), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_odoB), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_cavi), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD_guin), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_gut_AD), "total_core"))

plot_df_AD_gut_species = rbind(gut_AD_belli, gut_AD_sub, gut_AD_mili, gut_AD_odoA, gut_AD_odoB, gut_AD_cavi, gut_AD_guin) # rbindng all the above dataframes together

plot_df_AD_gut_species$species = factor(plot_df_AD_gut_species$species, levels = c("bellicosus", "subhyalinus", "cavithorax", "guinensis", 
                                                                                   "Species A", "Species B", "militaris-spiniger", "core"))
# Making the plot

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7",
         "core" = "black",  "total_core" = "black")

point_size = c("bellicosus" = 3 , "subhyalinus" = 3, 
         "cavithorax" = 3, "guinensis" = 3, 
         "Species A" = 3, "Species B" = 3, 
         "militaris-spiniger" = 3,
         "core" = 4,  "total_core" = 4)

s = c("bellicosus" = 16, "subhyalinus" = 16, 
           "cavithorax" = 16, "guinensis" = 16, 
           "Species A" = 16, "Species B" = 16, 
           "militaris-spiniger" = 16,
           "core" = 16,  "total_core" = 17)

alpha_values = c("bellicosus" = 0.7, "subhyalinus" = 0.7, 
                 "cavithorax" = 0.7, "guinensis" = 0.7, 
                 "Species A" = 0.7, "Species B" = 0.7, 
                 "militaris-spiniger" = 0.7,
                 "core" = 1,  "total_core" = 1)

core_plot_gut_species = ggplot(plot_df_AD_gut_species, aes(y = log10(mn), x = species, colour = colour, shape = colour))+
  geom_beeswarm(cex = 0.5, aes(size = colour, alpha = colour))+
  scale_color_manual(values = cols)+
  scale_shape_manual(values = s)+
  scale_size_manual(values = point_size)+
  scale_alpha_manual(values = alpha_values)+
  theme_pubclean()+
  labs(x = "Termite host species",
       y = "Log of mean normalised abundance")+
  scale_x_discrete(labels = c(expression(italic("M. bellicosus")),
                              expression(italic("M. subhyalinus")),
                              expression(italic("A. cavithorax")),
                              expression(italic("A. guinensis")),
                              expression(italic("O. species A")),
                              expression(italic("O. species B")),
                              expression(italic("P. militaris"))))+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        legend.position = "None")
ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 4/pannel A.pdf", core_plot_gut_species, width = 20, height = 12)


#### Get core gut AD end ####


#### Get core comb AD ####
#========================#

core_comb_AD = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C")) %>%
  select(!c("caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
  
### By species ###

## Macrotermes ##
core_comb_AD_belli = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "bellicosus") %>%
  select(!c("sample_type", "species", "caste_phase")) %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.90))
dim(core_comb_AD_belli)

core_comb_AD_sub = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  column_to_rownames(var = "rn") %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "subhyalinus") %>%
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

### Plotting comb cores ###
#========================#
# Here are are wanting to plot the core comb microbiomes in relation to every other OBU to show they are in general, highly abundant

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
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_belli), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_sub), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD), "total_core"))
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
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_mili), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_odoA), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_odoB), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD), "total_core"))

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
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_cavi), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD), "total_core"))

comb_AD_guin = feature_table_AD.norm %>%
  rownames_to_column(var = "rn") %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "sample_type", "species", "caste_phase")], by = c("rn" = "sample.id")) %>%
  filter(caste_phase %in% c("YC", "MC", "OC", "C"), species == "guinensis") %>% # filtering only works and correct genus
  select(!c("sample_type", "species", "caste_phase")) %>% # remove now unsated columns
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)) %>% # removing OBUs not present in this genus
  melt(.) %>% # melting the datframe
  select(!rn) %>% # removing sample name
  group_by(variable) %>% # group by OBU
  summarise(mn = mean(value)) %>% # getting the mean of each OBU per genus for OBUs present in that genus
  mutate(species = "guinensis") %>%
  mutate(colour = species) %>%
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD_guin), "core")) %>% # making a new column for colouring in the ggplot
  mutate(colour = replace(colour, variable %in% colnames(core_comb_AD), "total_core"))

plot_df_AD_comb_species = rbind(comb_AD_belli, comb_AD_sub, comb_AD_mili, comb_AD_odoA, comb_AD_odoB, comb_AD_cavi, comb_AD_guin) # rbindng all the above dataframes together

plot_df_AD_comb_species$species = factor(plot_df_AD_comb_species$species, levels = c("bellicosus", "subhyalinus", "cavithorax", "guinensis", 
                                                                                   "Species A", "Species B", "militaris-spiniger", "core"))
# Making the plot
cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7",
         "core" = "black",  "total_core" = "black")

point_size = c("bellicosus" = 3 , "subhyalinus" = 3, 
               "cavithorax" = 3, "guinensis" = 3, 
               "Species A" = 3, "Species B" = 3, 
               "militaris-spiniger" = 3,
               "core" = 4,  "total_core" = 4)

s = c("bellicosus" = 16, "subhyalinus" = 16, 
      "cavithorax" = 16, "guinensis" = 16, 
      "Species A" = 16, "Species B" = 16, 
      "militaris-spiniger" = 16,
      "core" = 16,  "total_core" = 17)

alpha_values = c("bellicosus" = 0.7, "subhyalinus" = 0.7, 
                 "cavithorax" = 0.7, "guinensis" = 0.7, 
                 "Species A" = 0.7, "Species B" = 0.7, 
                 "militaris-spiniger" = 0.7,
                 "core" = 1,  "total_core" = 1)

core_plot_comb_species = ggplot(plot_df_AD_comb_species, aes(y = log10(mn), x = species, colour = colour, shape = colour))+
  geom_beeswarm(cex = 0.5, aes(size = colour, alpha = colour))+
  scale_color_manual(values = cols)+
  scale_shape_manual(values = s)+
  scale_size_manual(values = point_size)+
  scale_alpha_manual(values = alpha_values)+
  theme_pubclean()+
  labs(x = "Termite host species",
       y = "Log of mean normalised abundance")+
  scale_x_discrete(labels = c(expression(italic("M. bellicosus")),
                              expression(italic("M. subhyalinus")),
                              expression(italic("A. cavithorax")),
                              expression(italic("A. guinensis")),
                              expression(italic("O. species A")),
                              expression(italic("O. species B")),
                              expression(italic("P. militaris"))))+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        legend.position = "None")
ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 4/pannel B.pdf", core_plot_comb_species, width = 20, height = 12)


#### Get core comb AD end ####


#### Core comb from gut ####
#==========================#
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

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")

gut2comb_df_plot = ggplot(gut2comb_df, aes(x = species, y = per_gut2comb, fill = species))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = cols)+
  theme_pubclean()+
  scale_x_discrete(labels = c(expression(italic("M. bellicosus")),
                              expression(italic("M. subhyalinus")),
                              expression(italic("A. cavithorax")),
                              expression(italic("A. guinensis")),
                              expression(italic("O. species A")),
                              expression(italic("O. species B")),
                              expression(italic("P. militaris"))))+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 13, angle = 30),
        legend.position = "None")+
  xlab("Species")+
  ylab("Percent of core comb OBUs also in termite guts")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 4/pannel c.pdf", gut2comb_df_plot)


#### Core comb from gut end ####