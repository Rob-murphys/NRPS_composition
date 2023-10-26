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

### Import Data ###
#=================#

## Un-normalized data ##

feature_table_AD = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/post_filt_feature_tables/AD_profiles_paired_redo_mean_OBU_table-filtered.rds") 

feature_table_16S = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/post_filt_feature_tables/16S_profiles_redo_mean_feature_table.rds") 

## Meatdata ##

metadata = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_noRepeats.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
  mutate(sample_type = case_when(caste_phase %in% c("YC", "MC", "OC", "C") ~ "Comb",
                                 caste_phase %in% c("SW", "LW", "SS", "LS", "S") ~ "Gut"))

# Sort metadata rows to match feature table
metadata.sorted = metadata[match(rownames(feature_table_AD), metadata$sample.id),]

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
  rename(chao1_16S = S.chao1) %>%
  filter(!is.na(species))
chao_df_16S$caste_phase = factor(chao_df_16S$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))

saveRDS(chao_df_16S, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/diversity_dataframes/16S_chao_richness.rds")

# Plot it

chao_plot_16S = ggplot(chao_df_16S, aes(x = caste_phase, y = chao1_16S, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  xlab("Sample type")+
  ylab("Chao 1 richness")+
  theme_pubclean()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave2("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2/Pannel A1.pdf", chao_plot_16S)


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
  rename(chao1_AD = S.chao1) %>%
  filter(!is.na(species))
chao_df_AD$caste_phase = factor(chao_df_AD$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))

saveRDS(chao_df_AD, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/diversity_dataframes/AD_chao_richness.rds")

# Plot it
chao_plot_AD = ggplot(chao_df_AD, aes(x = caste_phase, y = chao1_AD, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  xlab("Sample type")+
  ylab("Chao1 richness")+
  theme_pubclean()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave2("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2/Pannel A2.pdf", chao_plot_AD)


### Any differences between 16S and AD and then gut and comb within 16S and AD ###
t.test(chao_df_16S$chao1_16S, chao_df_AD$chao1_AD)

t.test(chao_df_16S[chao_df_16S$sample_type == "Gut",]$chao1_16S, chao_df_16S[chao_df_16S$sample_type == "Comb",]$chao1_16S)
t.test(chao_df_AD[chao_df_AD$sample_type == "Gut",]$chao1_AD, chao_df_AD[chao_df_AD$sample_type == "Comb",]$chao1_AD)

#### ANOVA of shannon diversity ####
#==================================#

## 16S ##

chao_16S.aov = aov(chao1_16S ~ caste_phase * species, data = chao_df_16S)
summary(chao_16S.aov)
chao_HSD_16S = TukeyHSD(chao_16S.aov, which = "caste_phase")
cohens_f(chao_16S.aov)

## AD ##
chao_AD.aov = aov(chao1_AD ~ caste_phase * species, data = chao_df_AD)
summary(chao_AD.aov)
chao_HSD_AD = TukeyHSD(chao_AD.aov, which = "caste_phase")
cohens_f(chao_AD.aov)


#### Linear modeling of the relationship between 16S and AD ####
#==============================================================#

## Merge the data together for plotting ##
chao_df_16S_values = chao_df_16S %>% 
  select(chao1_16S) %>%
  rownames_to_column(var = "rn")

chao_df = chao_df_AD %>% 
  rownames_to_column(var = "rn") %>%
  left_join(chao_df_16S_values, by = c("rn" = "rn"))

## The modeling ##
fit_gut = lm(chao1_AD ~ chao1_16S, data = chao_df[chao_df$sample_type == "Gut",])
summary(fit_gut)
confint(fit_gut, "chao1_16S", level=0.95)

fit_comb = lm(chao1_AD ~ chao1_16S, data = chao_df[chao_df$sample_type == "Comb",])
summary(fit_comb)
confint(fit_comb, "chao1_16S", level=0.95)

## Testing for a difference in the slopes ##
fit_all = lm(chao1_AD ~ chao1_16S * sample_type, data = chao_df)
anova(fit_all)
fit_all.trends = lstrends(fit_all, "sample_type", var = "chao1_16S")
pairs(fit_all.trends)

## Plotting these models ##
chao_AD_16S_lm_plot = ggplot(chao_df, aes(x = chao1_16S, y = chao1_AD, colour = sample_type))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Chao 1 richness for 16S ASVs")+
  ylab("Chao 1 richness for AD OBUs")+
  theme_pubclean()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 2/Pannel C1.pdf", chao_AD_16S_lm_plot)

#### Correlation between 16S and AD ####
#======================================#

metadata_noMicro = metadata %>%
  filter(genus != "Microtermes")

# Fungus Comb
cor.test(chao_df_AD[chao_df_AD$sample_type == "Comb",]$chao1_AD, chao_df_16S[chao_df_16S$sample_type == "Comb",]$chao1_16S, method = "pearson")

# Gut
cor.test(chao_df_AD[chao_df_AD$sample_type == "Gut",]$chao1_AD, chao_df_16S[chao_df_16S$sample_type == "Gut",]$chao1_16S, method = "pearson")

