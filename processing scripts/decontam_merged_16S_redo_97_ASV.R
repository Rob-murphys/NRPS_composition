setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project")

source("D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/R scripts/decontam_functions.R")

feature_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/pre_filt_feature_tables/16S_profiles_redo_97_ASV_merged_table-dada2.qza"
meta_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata.txt"
new_meta_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_metadata_decontam_16S.txt"
taxa_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/qiime_outputs/16S/redo/16S_profiles_redo_merged_seqs-taxonamy.qza"
# quant_path = "D:/OneDrive - University of Copenhagen/PhD/Projects/BGC_domain_amplicon_project/Metadata/BGC_sample_DNA_conc_16S.csv"

# new_meta_path = metadata_add(meta_path = meta_path, quant_path = quant_path)

physeq = physeq_import(feature_path = feature_path, new_meta_path = new_meta_path, taxa_path = taxa_path)

lib_size(physeq = physeq)

# Removing Mitochondria and chrolroplast
physeq.removed = remove_unwanted(physeq) # 86 removed

# Removing singletons and those with less that 100 reads
physeq.removed.low = remove_low(physeq.removed) # removed 78117 taxa

# Identifying and removing contaminats based on prevalence method
physeq.noncontam = id_contam(physeq.removed.low = physeq.removed.low) # identified 44 contaminants

generate_biom(physeq.noncontam = physeq.noncontam)




# pos_controls = c("batch-1-ve-pos", "batch-2-ve-pos", "batch-3-ve-pos", "batch-4-ve-pos", "batch-5-ve-pos", "batch-6-ve-pos")
# physeq.pos = subset_samples(physeq, (sample_names(physeq) %in% pos_controls))
# physeq.pos.low = filter_taxa(physeq.pos, function(x) sum(x > 0) > 1, TRUE)
# plot_bar(physeq.pos.low, fill = "Genus")
# 
# physeq.pos.low.Genus = tax_glom(physeq.pos.low, taxrank = "Genus") %>%
#   transform_sample_counts(., function(x) x / sum(x) ) %>%
#   psmelt()


