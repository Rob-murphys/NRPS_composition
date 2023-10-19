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

#### Degree plots ####

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


#### degree end plots ####