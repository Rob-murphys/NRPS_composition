library(dplyr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(reshape2)

#### Import data ####
#===================#

## Blast results
blast_res = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/dna-sequences_blastx_table.out", 
                       header = FALSE, sep = "\t", quote="")
colnames(blast_res) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

## AD
feature_table_AD.norm = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_normalised_OBU_table-filtered.rds")

## Meatdata ##

metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

metadata_oneS = metadata %>%
  mutate(caste_phase2 = caste_phase) %>%
  mutate(caste_phase2 = replace(caste_phase2, caste_phase2 %in% c("S", "SS", "LS"), "S"))

## Removing Microtermes ##
metadata_oneS_noMicro = metadata_oneS %>%
  filter(genus != "Microtermes")

feature_table_AD.norm_noMicro = feature_table_AD.norm %>%
  filter(rownames(.) %in% metadata_oneS_noMicro$sample.id)
#### Import data end ####

#### Recover best hit per sequence based on evalue ####
#=====================================================#
blast_res_singles = NULL
for(id in unique(blast_res$qseqid)){
  df = blast_res %>%
    filter(qseqid == id) %>%
    arrange(evalue)
  blast_res_singles = rbind(blast_res_singles, df[1,])
}

# Filter out evalues below e-10 to ensure we are conservative in our estimation of homology
blast_res_singles_filt = blast_res_singles %>%
  filter(bitscore >= 50)
#### Recover best hit per sequence based on evalue end ####


#### Create sample type specific OBU lists ####
SW_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 == "SW") %>% # filtering sample type
  select(!c("caste_phase2")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

LW_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 == "LW") %>% # filtering sample type
  select(!c("caste_phase2")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

S_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 == "S") %>% # filtering sample type
  select(!c("caste_phase2")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

YC_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 == "YC") %>% # filtering sample type
  select(!c("caste_phase2")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

MC_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 == "MC") %>% # filtering sample type
  select(!c("caste_phase2")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

OC_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 == "OC") %>% # filtering sample type
  select(!c("caste_phase2")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

C_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "caste_phase2")], by = c("rowname" = "sample.id")) %>%
  filter(caste_phase2 == "C") %>% # filtering sample type
  select(!c("caste_phase2")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus
##### ####

#### Create species type specific OBU lists ####
belli_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "bellicosus") %>% # filtering sample type
  select(!c("species")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

sub_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "subhyalinus") %>% # filtering sample type
  select(!c("species")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

cavi_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "cavithorax") %>% # filtering sample type
  select(!c("species")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

guin_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "guinensis") %>% # filtering sample type
  select(!c("species")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

odoA_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species A") %>% # filtering sample type
  select(!c("species")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

odoB_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "Species B") %>% # filtering sample type
  select(!c("species")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

mili_OBUs = colnames(feature_table_AD.norm_noMicro %>%
  rownames_to_column() %>%
  left_join(metadata_oneS_noMicro[colnames(metadata_oneS_noMicro) %in% c("sample.id", "species")], by = c("rowname" = "sample.id")) %>%
  filter(species == "militaris-spiniger") %>% # filtering sample type
  select(!c("species")) %>% # remove columns
  column_to_rownames() %>%
  select(where(~ sum(. != 0, na.rm = TRUE) / length(.) > 0.0)))  # removing OBUs not present in this genus

##### ####


#### Group blast results by bit score ####
#========================================#

# Grouping by bitscore starting at 50 as this almost always infers homology
df_id_filt = blast_res_singles_filt %>%
  filter(pident < 50)

df_id_filt_50 = blast_res_singles_filt %>%
  filter(pident >= 50, pident < 70)

df_id_filt_70 = blast_res_singles_filt %>%
  filter(pident >= 70, pident < 85 )

df_id_filt_85 = blast_res_singles_filt %>%
  filter(pident >= 85, pident < 95)

df_id_filt_95 = blast_res_singles_filt %>%
  filter(pident >= 95, pident <= 100)

#### Create sample type specific grouped blast result dataframe ####
#==================================================================#
OBU_lists = list(SW_OBUs, LW_OBUs, S_OBUs, YC_OBUs, MC_OBUs, OC_OBUs, C_OBUs)
names(OBU_lists) = c("SW", "LW", "S", "YC", "MC", "OC", "C")

grouped_id_df = data.frame(id = c(0, 50,70,85,95), SW=NA, LW=NA, S=NA, YC=NA, MC=NA, OC=NA, C=NA)
for(type in 1:length(unique(metadata_oneS_noMicro$caste_phase2))){
  grouped_id_df[grouped_id_df$id == 0, type+1] = round((nrow(df_id_filt[df_id_filt$qseqid %in% OBU_lists[[type]],]) / length(OBU_lists[[type]]))*100, 2)
  grouped_id_df[grouped_id_df$id == 50, type+1] = round((nrow(df_id_filt_50[df_id_filt_50$qseqid %in% OBU_lists[[type]],]) / length(OBU_lists[[type]]))*100, 2)
  grouped_id_df[grouped_id_df$id == 70, type+1] = round((nrow(df_id_filt_70[df_id_filt_70$qseqid %in% OBU_lists[[type]],]) / length(OBU_lists[[type]]))*100, 2)
  grouped_id_df[grouped_id_df$id == 85, type+1] = round((nrow(df_id_filt_85[df_id_filt_85$qseqid %in% OBU_lists[[type]],]) / length(OBU_lists[[type]]))*100, 2)
  grouped_id_df[grouped_id_df$id == 95, type+1] = round((nrow(df_id_filt_95[df_id_filt_95$qseqid %in% OBU_lists[[type]],]) / length(OBU_lists[[type]]))*100, 2)
}

plot_df = melt(grouped_id_df, id="id")

ggplot(plot_df, aes(x=variable, y = value, fill = as.factor(id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_pubclean()
fill_cols = c("0" = "#c85b82", "50" = "#83a344", "70" = "#8d72c9", "85" = "#c86f3c", "95" = "#4aac8d")

sample_type_novelty_plot = ggplot(plot_df, aes(x=variable, y = value, fill = as.factor(id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_pubclean()+
  scale_fill_manual(values = fill_cols)+
  labs(x = "Sample type",
       y = "Percent of total OBUS in sample type")+
  scale_x_discrete(labels = c("Small worker", "Large worker", "Soldier", "Young comb", "Mature comb", "Old comb", "Comb"))+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        legend.position = "None",
        axis.text.x = element_text(angle = 30))
ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 4/pannel C.pdf", sample_type_novelty_plot)

#### ####

#### Do core OBUs have hits? ####

## Species cores

# list of core OBUs. Requiers running part of core analysis script


core_OBUs = c(colnames(core_gut_AD_belli), colnames(core_gut_AD_sub), 
              colnames(core_gut_AD_cavi), colnames(core_gut_AD_guin), 
              colnames(core_gut_AD_odoA), colnames(core_gut_AD_odoB), 
              colnames(core_gut_AD_mili),
              colnames(core_comb_AD_belli), colnames(core_comb_AD_sub), 
              colnames(core_comb_AD_cavi), colnames(core_comb_AD_guin), 
              colnames(core_comb_AD_odoA), colnames(core_comb_AD_odoB), 
              colnames(core_comb_AD_mili))

core_blast_hits = blast_res_singles_filt[blast_res_singles_filt$qseqid %in% core_OBUs,]

core_blast_hits_filt = core_blast_hits %>%
  filter(pident >= 50)

write.csv(core_blast_hits_filt, "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/core_blast_hits_filt.csv", row.names = F)

nrow(core_blast_hits_filt) / length(core_OBUs)

## Resampling random OBUs from list of all OBUs to see how frequently they have hits
all_OBUs = colnames(feature_table_AD.norm_noMicro)
novelty_random = NULL
for(x in 1:1000){
  sample_OBUs = sample(all_OBUs, length(core_OBUs), replace = FALSE)
  df = blast_res_singles_filt %>%
    filter(pident >= 50)
  y = sum(sample_OBUs %in% df$qseqid) / length(sample_OBUs)
  novelty_random = append(novelty_random, y)
}
mean(novelty_random)


## Totlal core

core_blast_hits = blast_res_singles_filt[blast_res_singles_filt$qseqid %in% colnames(total_AD_core_comb),]


#### Do core OBUs have hits end ####

#### Checking for already identified NRPS or hybrids ####
#=======================================================#

actinomycin = blast_res_singles_filt[grep("BGC0000296", blast_res_singles_filt$sseqid),]

bacillaene = blast_res_singles_filt[grep("BGC0001089", blast_res_singles_filt$sseqid),]

dentigerumycin = blast_res_singles_filt[grep("BGC0002735", blast_res_singles_filt$sseqid),]

microtermolide = blast_res_singles_filt[grep("BGC0001332", blast_res_singles_filt$sseqid),]

macrotermycins = blast_res_singles_filt[grep("BGC0001658", blast_res_singles_filt$sseqid),] # nothing
