library(tibble)
library(ggplot2)
library(dplyr)
library(vegan)
library(cowplot)
library(lsmeans)
library(rstatix)
library(ggbeeswarm)
library(effectsize)
library(ggpubr)

#### Import Data ####
#===================#

## Un-normalized data ##

feature_table_AD = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_OBU_table-filtered.rds") 

feature_table_16S = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_feature_table.rds") 

## Meatdata ##

metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  mutate(sample_type = case_when(caste_phase %in% c("YC", "MC", "OC", "C") ~ "Comb",
                                 caste_phase %in% c("SW", "LW", "SS", "LS", "S") ~ "Gut"))

levels(metadata$caste_phase) = c(levels(metadata$caste_phase), "W")

metadata.sorted_oneS_W = metadata %>%
  mutate(caste_phase2 = caste_phase) %>%
  mutate(caste_phase2 = replace(caste_phase2, caste_phase2 %in% c("S", "SS", "LS"), "S"),
         caste_phase2 = replace(caste_phase2, caste_phase2 %in% c("SW", "LW"), "W"))

## Removing Microtermes ##

feature_table_16S_noMicro = feature_table_16S %>%
  rownames_to_column() %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "genus")], by = c("rowname" = "sample.id")) %>%
  column_to_rownames() %>%
  filter(genus != "Microtermes") %>%
  dplyr::select(!genus)

feature_table_AD_noMicro = feature_table_AD %>%
  rownames_to_column() %>%
  left_join(metadata[colnames(metadata) %in% c("sample.id", "genus")], by = c("rowname" = "sample.id")) %>%
  column_to_rownames() %>%
  filter(genus != "Microtermes") %>%
  dplyr::select(!genus)

#### Import Data end ####


#### Calculate Shannon diversity ####
#===================================#

## 16S ##
shannon_df_16S = as.data.frame(vegan::diversity(feature_table_16S_noMicro, index = "shannon")) %>%
  rename(shannon = `vegan::diversity(feature_table_16S_noMicro, index = "shannon")`) %>%
  rownames_to_column() %>%
  left_join(metadata.sorted_oneS_W, by = c("rowname" = "sample.id")) %>%
  column_to_rownames()
shannon_df_16S$caste_phase = factor(shannon_df_16S$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))
shannon_df_16S$caste_phase2 = factor(shannon_df_16S$caste_phase2, levels = c("C", "YC", "MC", "OC", "W", "S"))

## AD ##
shannon_df_AD = as.data.frame(vegan::diversity(feature_table_AD_noMicro, index = "shannon")) %>%
  rename(shannon = `vegan::diversity(feature_table_AD_noMicro, index = "shannon")`) %>%
  rownames_to_column() %>%
  left_join(metadata.sorted_oneS_W, by = c("rowname" = "sample.id")) %>%
  column_to_rownames(var = "rowname")
shannon_df_AD$caste_phase = factor(shannon_df_AD$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))
shannon_df_AD$caste_phase2 = factor(shannon_df_AD$caste_phase2, levels = c("C", "YC", "MC", "OC", "W", "S"))


#### Chao 1 richness ####
#=======================#

## 16S ##
chao_df_16S = feature_table_16S_noMicro %>%
  round() %>%
  estimateR(.) %>%
  as.data.frame() %>%
  dplyr::slice(.data=., 2) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample.id") %>%
  left_join(metadata, by = c("sample.id" = "sample.id")) %>%
  column_to_rownames(var = "sample.id") %>%
  rename(chao1 = S.chao1) %>%
  filter(!is.na(species))
chao_df_16S$caste_phase = factor(chao_df_16S$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))

## AD ##
chao_df_AD = feature_table_AD_noMicro %>%
  round() %>%
  estimateR(.) %>%
  as.data.frame() %>%
  dplyr::slice(.data=., 2) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample.id") %>%
  left_join(metadata, by = c("sample.id" = "sample.id")) %>%
  column_to_rownames(var = "sample.id") %>%
  rename(chao1 = S.chao1) %>%
  filter(!is.na(species))
chao_df_AD$caste_phase = factor(chao_df_AD$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))

#### Shannon plot ####
cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")
shannon_df_16S$data_type = "16S"
shannon_df_AD$data_type = "AD"

shannon_combind = rbind(shannon_df_16S, shannon_df_AD)

shan_plot = ggplot(shannon_combind, aes(x = caste_phase, y = shannon, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  facet_wrap(~data_type, nrow = 1, scales = "fixed")+
  xlab("Sample type")+
  ylab("Shannon diversity")+
  scale_colour_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none", axis.title.y = element_text(size = 18))

#### Shannon plot ####
cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")
chao_df_16S$data_type = "16S"
chao_df_AD$data_type = "AD"

chao_combind = rbind(chao_df_16S, chao_df_AD)

chao_plot = ggplot(chao_combind, aes(x = caste_phase, y = chao1, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  facet_wrap(~data_type, nrow = 1, scales = "fixed")+
  xlab("Sample type")+
  ylab("Chao1 richness")+
  scale_colour_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

#### Cobmine the plots ####

combined_plot = ggarrange(shan_plot, chao_plot, nrow = 2, align = "v")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure S2/combined_plot.pdf", combined_plot, width = 15, height = 8)
