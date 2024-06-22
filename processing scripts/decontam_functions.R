library(qiime2R)
library(phyloseq)
library(decontam)
library(ggplot2)
library(dplyr)
library(biomformat)

metadata_add = function(meta_path, quant_path){
  # add DNA quant and if a control of true sample to metadata
  metadata = read.csv(meta_path, sep = "\t") # read in metadata
  quantdata = read.csv(quant_path) # dna quant data
  lookup = c(quant_reading = "X16S.total.DNA.in.ng." , quant_reading = "AD..total.DNA.in.ng.", "sample-id" = "sample.id", "barcode-number" = "barcode.number", "barcode-sequence" = "barcode.sequence", "caste_phase" = "caste.phase") # dictionary to change col names
  # Left join metadata and quant data then add true sample or control column and rename columns
  metadata = metadata %>% left_join(quantdata, by = c("sample.id" = "sample.id")) %>% # left join
    mutate(Sample_or_Control = case_when(nest == "negative-control" ~ "Control Sample", # Create new column with control or sample
                                         nest == "positive-control" ~ "Control Sample",
                                         nest == "environmental"  ~ "Control Sample",
                                         TRUE ~ "True Sample")) %>%
    rename(any_of(lookup)) %>% # rename columns
    mutate(`barcode-number` = as.character(`barcode-number`), quant_reading = as.character(quant_reading)) %>% # Need columns to be characters to add the q2:types
    add_row(.before = 1 ,`sample-id` = "#q2:types", `barcode-number` = "Categorical", `barcode-sequence` = "Categorical", # Add q2:types
            genus = "Categorical", species = "Categorical", nest = "Categorical", `caste_phase` = "Categorical", 
            quant_reading = "Numeric", Sample_or_Control = "Categorical")
  
  new_meta_path = paste(strsplit(meta_path, split = ".txt")[[1]][1], "decontam.txt", sep = "_") # creating a path for the decontam metadata
  write.table(file = new_meta_path, x = metadata, sep = "\t", row.names = FALSE, quote = FALSE) # writing to file as phyloseq needs to read from file
  return(new_meta_path)
}

physeq_import = function(feature_path, new_meta_path, taxa_path){
  if(taxa_path != ""){
  # Importing qiime2 object as phyloseq object
  physeq = qza_to_phyloseq(
    features = feature_path,
    metadata = new_meta_path,
    taxonomy = taxa_path)
  }
  else{
    physeq = qza_to_phyloseq(
      features = feature_path,
      metadata = new_meta_path)
  }
  physeq = subset_samples(physeq, !(sample_names(physeq) %in% c("batch-1-ve-pos", "batch-2-ve-pos", "batch-3-ve-pos", "batch-4-ve-pos", "batch-5-ve-pos", "batch-6-ve-pos"))) # removing positive controls as they can't be used in any analysis
  return(physeq)
}

lib_size = function(physeq){
  # Inspect library sizes
  df = as.data.frame(sample_data(physeq))# Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize = sample_sums(physeq)
  df = df[order(df$LibrarySize),]
  df$Index = seq(nrow(df))
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + 
    geom_point()
}

id_contam = function(physeq.removed.low){
  # Identify contaminants with both frequence and prevelancy methods
  sample_data(physeq.removed.low)$is.neg = sample_data(physeq.removed.low)$Sample_or_Control == "Control Sample"
  contamdf.either = isContaminant(physeq.removed.low, method = "prevalence", neg="is.neg")
  print("Contaminant overview:")
  print(table(contamdf.either$contaminant))
  physeq.noncontam = prune_taxa(!contamdf.either$contaminant, physeq.removed.low) # pruning contaninant ASVs
  return(physeq.noncontam)
}

generate_biom = function(physeq.noncontam){
  # Write phyloseq object to biom so we can read it back into qiime2
  otu = as(otu_table(physeq.noncontam),"matrix")
  otu_biom = biomformat::make_biom(data = otu)
  biom_file = paste(strsplit(feature_path, ".qza")[[1]][1],"decontam.biom", sep = "_")
  biomformat::write_biom(otu_biom, biom_file)
}

remove_unwanted = function(physeq){
  physeq.removed = subset_taxa(physeq, (Order != "Chloroplast"|is.na(Order))) %>%
    subset_taxa((Family != "Mitochondria"|is.na(Family)))
  cat(paste0("Original: ", ntaxa(physeq),
               "\nCleaned: ", ntaxa(physeq.removed),
               "\nRemoved: ", ntaxa(physeq) - ntaxa(physeq.removed)))
  return(physeq.removed)
}

remove_low = function(physeq.removed){
  physeq.removed.low = physeq.removed %>% subset_samples(., !Sample_or_Control == "Control Sample") %>%  
    # Removing singleton taxa. It is asking how many times is the read count greater than 0 and if that is less that one remove it
    filter_taxa(., function(x) sum(x > 0) > 1, prune=TRUE) %>% # Sum the number of times x (x being a taxa) has a count greater than 0. if that is not greater than 1 then remove it
    # Removing features with less that 100 reads by asking how many times this feature has more 100 reads, and if that occurs at least once, keep it
    filter_taxa(., function (x) {sum(x > 100) > 0}, prune=TRUE) %>% # Remove features with less that 100 reads
    # Pruning the original object of these taxa. We do this as the object output from filter_taxa is without controls because we subset them out for this step
    taxa_names() %>%
    prune_taxa(physeq.removed)
  cat(paste0("Original: ", ntaxa(physeq.removed),
             "\nCleaned: ", ntaxa(physeq.removed.low),
             "\nRemoved: ", ntaxa(physeq) - ntaxa(physeq.removed.low)))
  return(physeq.removed.low)
}