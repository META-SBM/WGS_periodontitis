library(phyloseq)
library(tidyverse)
library(tibble)
# metadata table
meta <- read_excel('../porph_all_data_tables/all_metadata.xlsx')
meta$ID <- gsub('MPAR','MPAR_',meta$ID)
meta <- as.data.frame(sapply(meta, function(x) replace(x, x == "", "no_data")))
colnames(meta) <- gsub(" ", "_", colnames(meta))
row.names(meta) <-meta$ID
row.names(meta) <-meta$ID
names(meta)[duplicated(names(meta))]
meta <- meta %>%
  mutate(Age = as.numeric(Age),  # Convert age to numeric
         Age_group = case_when(
           between(Age, 19, 35) ~ '19_35',
           between(Age, 36, 45) ~ '36_45',
           between(Age, 46, 55) ~ '46_55',
           between(Age, 56, 65) ~ '56_65',
           is.na(Age) ~ 'no_data',  
           TRUE ~ 'above_65'
         ))

numeric_cols <- sapply(meta, is.numeric)
meta[numeric_cols] <- lapply(meta[numeric_cols], as.factor)
calculate_bmi <- function(weight, height) {
  if (is.na(weight) || is.na(height)) {
    return(NA)  # Skip calculation if either value is NA
  } else {
    bmi <- weight / ((height/100)^2)  # Calculate BMI only if both values are provided
    return(bmi)
  }
}
meta$BMI <- mapply(calculate_bmi, as.numeric(meta$Weight), as.numeric(meta$Height))
meta <- meta %>%
  mutate(BMI = as.numeric(BMI),  # Convert age to numeric
         BMI_group = case_when(
           between(BMI, 18.5, 24.9999) ~ 'normal',
           between(BMI, 25, 29.99999) ~ 'overweight',
           between(BMI, 30, 99) ~ 'obesity',
           between(BMI, 0, 18.49999) ~ 'underweight',
           is.na(BMI) ~ 'no_data',  
           TRUE ~ 'no_data'
         ))
# WGS otu and tax table
otu <- read.csv('./porph_all_data_tables/merged_data_raw (4).tsv', sep='\t' )
row.names(otu) <- otu$X
colnames(otu) <- gsub("\\.", "_", colnames(otu))
otu$X <- NULL
tax <- data.frame(matrix(ncol = 6 , nrow = length(colnames(otu))))
colnames(tax) <- c("Kingdom","Phylum", "Family", "Order", "Genus", "Species")
tax$Species = colnames(otu)
row.names(tax) = tax$Species
tax <- tax %>%
  separate(col = Species, into = c("Genus", "Species"), sep = "_")
tax$Species = colnames(otu)

# create WGS phyloseq object
ps <- phyloseq(otu_table(as.data.frame(as.matrix(otu)),taxa_are_rows = F),
               sample_data(as.data.frame(as.matrix(meta))),
               tax_table(as.matrix(tax)))

# filtering WGS phyloseq
ps <- subset_taxa(ps, !(Kingdom %in% c("Archaea", "Eukaryota")))
sample_data(ps) <- meta
samples_2_del = c('MPAR_80','MPAR_87') # it's easier to remove them here rather then clearing previous tables.
ps <- subset_samples(ps, !(ID %in% samples_2_del))
saveRDS(ps, './porph_all_data_tables/ps_wgs.rds')

#create HUMANn phyloseq
humann_data <- read.csv("../porph_all_data_tables/MPAR_all_abundancies.tsv", sep='\t')
humann_data <- humann_data %>%
  column_to_rownames(var = "Pathway")
humann_data <- humann_data[!grepl("UNINTEGRATED|UNMAPPED", rownames(humann_data)), ]
humann_data <- humann_data[!grepl("archaea|fungi|yeast|mammals|vertebrates|plants|mitochondria|invertebrates|mammalian", rownames(humann_data)), ]

# making also HUMANn joint phyloseq, without splitting pathways to specific bacteria
humann_data_joint <- humann_data[!grepl("\\|", row.names(humann_data)), ]  # Rows without '|'
hm_phyloseq <- phyloseq(otu_table(humann_data,taxa_are_rows = TRUE), sample_data(meta))
hm_phyloseq_joint <- phyloseq(otu_table(humann_data_joint,taxa_are_rows = TRUE), sample_data(meta))
hm_phyloseq <- subset_samples(hm_phyloseq, !(ID %in% samples_2_del))
hm_phyloseq_joint <- subset_samples(hm_phyloseq_joint, !(ID %in% samples_2_del))
saveRDS(hm_phyloseq_joint, '/home/ignatsonets/ps_humann_joint.rds')
saveRDS(hm_phyloseq, '/home/ignatsonets/ps_humann.rds')

# MAGs phyloseq and its filtering
mags_count_data <- read.csv("../porph_all_data_tables/instrain_coverage_res.csv", sep=',',row.names = 1)
mags_tax_data <- read.csv("../porph_all_data_tables/gtdb_v2.tsv", sep=',')
mags_tax_data$Species <- gsub(" ", "_", mags_tax_data$Species)
mags_tax_data<- mags_tax_data %>% mutate_all(na_if,"")
tax_matrix <- as.matrix(mags_tax_data[, -1])
rownames(tax_matrix) <- mags_tax_data[, 1]    # 1st col is rownames now
colnames(tax_matrix) <- colnames(mags_tax_data)[-1]
mags_tax_table <- tax_table(tax_matrix)
mags_ps <- phyloseq(otu_table(mags_count_data,taxa_are_rows = TRUE), sample_data(meta),tax_table(mags_tax_table))
mags_ps <- subset_samples(mags_ps, !(ID %in% samples_2_del))
mags_ps@otu_table[is.na(mags_ps@otu_table)] <- 0
# filtering out of bad quality MAGs
keep_genomes <- read.csv("/home/ignatsonets/Downloads/porph_all_data_tables/mpar_filtered_mags.txt")
keep_genomes <- as.character(keep_genomes)
keep_genomes_clean <- gsub('[\\"c()]', '', keep_genomes)  # Remove \", "c(", ")", etc.
keep_genomes_clean <- trimws(keep_genomes_clean)
keep_genomes_clean <- unlist(strsplit(keep_genomes_clean, ", "))
if (!taxa_are_rows(mags_ps)) {
  otu_table(mags_ps) <- t(otu_table(mags_ps))
}
# Subset taxa (genomes)
mags_ps_filtered <- prune_taxa(taxa_names(mags_ps) %in% keep_genomes_clean, mags_ps)
saveRDS(mags_ps_filtered, '/home/ignatsonets/ps_mags_filt.rds')
