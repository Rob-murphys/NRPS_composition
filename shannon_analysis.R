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
#=====================================#

### 16S ###

# Calculate the diversity index
shannon_df_16S = as.data.frame(vegan::diversity(feature_table_16S_noMicro, index = "shannon")) %>%
  rename(shannon_16S = `vegan::diversity(feature_table_16S_noMicro, index = "shannon")`) %>%
  rownames_to_column() %>%
  left_join(metadata.sorted_oneS_W, by = c("rowname" = "sample.id")) %>%
  column_to_rownames()
shannon_df_16S$caste_phase = factor(shannon_df_16S$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))
shannon_df_16S$caste_phase2 = factor(shannon_df_16S$caste_phase2, levels = c("C", "YC", "MC", "OC", "W", "S"))

saveRDS(shannon_df_16S,  file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/diversity_dataframes/16S_shannon_richness.rds")

## Figure 1 plot ##
cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7")

shan_plot_16S = ggplot(shannon_df_16S, aes(x = caste_phase2, y = shannon_16S, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  xlab("Sample type")+
  ylab("Shannon diversity")+
  scale_colour_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
      axis.text.x = element_text(colour = "black", size = 18), 
      legend.position = "none", axis.title.y = element_text(size = 18), 
      axis.title.x = element_text(size = 18, colour = "black"))

ggsave2("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 1/pannel D1.pdf", shan_plot_16S)

## Supplamentary plot ##
shan_plot_16S = ggplot(shannon_df_16S, aes(x = caste_phase, y = shannon_16S, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  xlab("Sample type")+
  ylab("Shannon diversity")+
  scale_colour_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave2("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure S1 A.pdf", shan_plot_16S)

### AD ###

# Calculate the diversity index
shannon_df_AD = as.data.frame(vegan::diversity(feature_table_AD_noMicro, index = "shannon")) %>%
  rename(shannon_AD = `vegan::diversity(feature_table_AD_noMicro, index = "shannon")`) %>%
  rownames_to_column() %>%
  left_join(metadata.sorted_oneS_W, by = c("rowname" = "sample.id")) %>%
  column_to_rownames(var = "rowname")
shannon_df_AD$caste_phase = factor(shannon_df_AD$caste_phase, levels = c("C", "YC", "MC", "OC", "SW", "LW", "SS", "LS", "S"))
shannon_df_AD$caste_phase2 = factor(shannon_df_AD$caste_phase2, levels = c("C", "YC", "MC", "OC", "W", "S"))

saveRDS(shannon_df_AD,  file = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/AD/redo/diversity_dataframes/AD_shannon_richness.rds")

## Figure 1 plot ##
shan_plot_AD = ggplot(shannon_df_AD, aes(x = caste_phase2, y = shannon_AD, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  xlab("Sample type")+
  ylab("Shannon diversity")+
  scale_colour_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))
ggsave2("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 1/pannel D2.pdf", shan_plot_AD)

## Supplamentary plot ##
# Plot it

shan_plot_AD = ggplot(shannon_df_AD, aes(x = caste_phase, y = shannon_AD, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  xlab("Sample type")+
  ylab("Shannon diversity")+
  scale_colour_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))
ggsave2("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure S1/pannel B2.pdf", shan_plot_AD)

#### Combined 16S and AD plot ####
comb_16S_df = shannon_df_16S %>%
  rename(shannon = shannon_16S) %>%
  mutate(data = "16S")

comb_AD_df = shannon_df_AD %>%
  rename(shannon = shannon_AD) %>%
  mutate(data = "AD")

comb_16S_AD_df = rbind(comb_16S_df, comb_AD_df)

comb_16S_AD = ggplot(comb_16S_AD_df, aes(x = caste_phase2, y = shannon, shape = sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = species), size = 3)+
  facet_wrap(~data, scales = "fixed", nrow = 1)+
  xlab("Sample type")+
  ylab("Shannon diversity")+
  scale_colour_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave2("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 1/pannel D1-2.pdf", comb_16S_AD, width = 18.42, height = 9.83)


#### ANOVA of shannon diversity ####
#==================================#

## 16S ##
shannon_16S.aov = aov(shannon_16S ~ caste_phase * species, data = shannon_df_16S)
summary(shannon_16S.aov)
shannon_HSD_16S = TukeyHSD(shannon_16S.aov, which = "caste_phase")
cohens_f(shannon_16S.aov)

## AD ##
shannon_AD.aov = aov(shannon_AD ~ caste_phase * species, data = shannon_df_AD)
summary(shannon_AD.aov)
shannon_HSD_AD = TukeyHSD(shannon_AD.aov, which = "caste_phase")
cohens_f(shannon_AD.aov)


#### Linear modeling of the relationship between 16S and AD ####
#==============================================================#

## Merge the data together for plotting ##
shannon_df_16S_values = shannon_df_16S %>% 
  dplyr::select(shannon_16S) %>%
  rownames_to_column(var = "rn")

shannon_df = shannon_df_AD %>% 
  rownames_to_column(var = "rn") %>%
  left_join(shannon_df_16S_values, by = c("rn" = "rn")) %>%
  column_to_rownames(var = "rn")

## The modeling ##

fit_shannon_gut = lm(shannon_AD ~ shannon_16S, data = shannon_df[shannon_df$sample_type == "Gut",])
summary(fit_shannon_gut)
coef(fit_shannon_gut)
confint(fit_shannon_gut, "shannon_16S", level=0.95)

fit_shannon_comb = lm(shannon_AD ~ shannon_16S, data = shannon_df[shannon_df$sample_type == "Comb",])
summary(fit_shannon_comb)
coef(fit_shannon_comb)
confint(fit_shannon_comb, "shannon_16S", level=0.95)

fit_all = lm(shannon_AD ~ shannon_16S * sample_type, data = shannon_df)
anova(fit_all)
fit_all.trends = lstrends(fit_all, "sample_type", var = "shannon_16S")
pairs(fit_all.trends)

## Plotting the data ##

cols = c("bellicosus" = "#E77675", "subhyalinus" = "#A31E23", 
         "cavithorax" = "#E3A86C", "guinensis" = "#BA7328", 
         "Species A" = "#91CB7B", "Species B" = "#3B7538", 
         "militaris-spiniger" = "#8499C7",
         "Gut" = "#D7DF23", "Comb" = "#8B5E3C")

shannon_AD_16S_lm_plot = ggplot(shannon_df, aes(x = shannon_16S, y = shannon_AD, colour = species, shape = sample_type))+
  geom_point(size = 2.5)+
  scale_colour_manual(values = cols) +
  geom_smooth(method = "lm", formula = y ~ x, aes(x = shannon_16S, y = shannon_AD, colour = sample_type), inherit.aes = F, linewidth = 2, se = T)+
  xlab("Shannon diversity for 16S ASVs")+
  ylab("Shannon diversity for AD OBUs")+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Figures/figure 1/pannel D3.pdf", shannon_AD_16S_lm_plot)

#### Correlation between 16S and AD ####
#======================================#
metadata_noMicro = metadata %>%
  filter(genus != "Microtermes")

# Fungus Comb
cor.test(shannon_df_AD[shannon_df_AD$sample_type == "Comb",]$shannon_AD, shannon_df_16S[shannon_df_16S$sample_type == "Comb",]$shannon_16S, method = "pearson")

# Gut
cor.test(shannon_df_AD[shannon_df_AD$sample_type == "Gut",]$shannon_AD, shannon_df_16S[shannon_df_16S$sample_type == "Gut",]$shannon_16S, method = "pearson")

# Total
cor.test(shannon_df_AD$shannon_AD, shannon_df_16S$shannon_16S, method = "pearson")
