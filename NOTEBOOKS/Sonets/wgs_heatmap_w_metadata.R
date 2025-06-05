library(phyloseq)
library(microbiome)
library(pheatmap)

#Read phyloseq
ps <-  readRDS('/home/ignatsonets/ps_wgs.rds')
ps_obj <- ps

#Function for making count table

prepare_count_table <- function(ps_obj, taxa_name_func = function(x) {paste0(x['ASV'], '__', 'p__', x['Phylum'], '.f__' , x['Family'], '.g__' , gsub('/', '', gsub("-", "", x['Genus'])))}) {
  otu_matrix <- as(otu_table(ps_obj), 'matrix')
  if (taxa_are_rows(ps_obj)) {
    otu_matrix <- t(otu_matrix)
  }
  taxa_matrix <- as(tax_table(ps_obj), 'matrix')
  taxa_matrix <- cbind(ASV = rownames(taxa_matrix), taxa_matrix)  
  taxa_matrix_good_names <- apply(taxa_matrix, MARGIN = 1, taxa_name_func)
  colnames(otu_matrix) <- taxa_matrix_good_names
  return(otu_matrix)
}

# ==============================================================================
# Filter specimens
print('After specimen selection:')
ps_obj
# ==============================================================================
# Filter taxa
prevalence <- 10/100
detection <- 1
taxas <- core_members(ps_obj, detection = detection, prevalence = prevalence, include.lowest = F)
ps_obj <- prune_taxa(taxas, ps_obj)

print('After taxa filtering:')
ps_obj


# ==============================================================================
# Agglomeration
tax_glom_level <- 'Species'
if (tax_glom_level != 'ASV'){
  ps_obj <- speedyseq::tax_glom(ps_obj, tax_glom_level, NArm =T)
}

print('After agglomeration:')
ps_obj


# ==============================================================================
# Transform counts
transformation <- 'clr'
ps_obj <- microbiome::transform(ps_obj, transformation)

taxa_name_func <- function(x) {paste0(x['ASV'],  '_' , gsub('/', '', gsub("-", "", x['Species'])))}
# ==============================================================================
# Prepare metadata table
ps_meta <- as(sample_data(ps_obj), 'data.frame')
cols_of_interest =c('Antibiotics_treatment_in_last_6_months','Treatment_of_periodontitis', 'Sex', 'Age_group', 'BMI_group', 'Stage_of_severity_of_periodontitis', 'Smoking_status','Arterial hypertension')
meta_for_heatmap <- ps_meta[, (names(ps_meta) %in% cols_of_interest)]
count_mtrx <- prepare_count_table(ps_obj, taxa_name_func)
ordered_counts <-count_mtrx[rownames(meta_for_heatmap),]

# ==============================================================================
# Plot Heatmap
clustering_distance_rows <- 'euclidean'
# clustering_distance_rows <- 'binary'
# clustering_distance_rows <- 'correlation'

clustering_distance_cols <- 'euclidean'
# clustering_distance_cols <- 'binary'
# clustering_distance_cols <- 'correlation'

save_path <- '/home/ignatsonets/TEST_wgs_hm_md_v3_'
filename <- paste0(save_path, 'WGS_tax_glom_level_no_MPAR_8087', tax_glom_level, '_transform_', transformation, '_det', detection, '_prev',prevalence,  '_clust_dist_cols_', clustering_distance_cols,  '_clust_dist_rows_', clustering_distance_rows, '.pdf')
# Heatmap
p <-pheatmap(t(count_mtrx),
             cluster_rows = T, treeheight_row  = 70, treeheight_col = 70,
             cluster_cols = T,
             fontsize=14,
             annotation_col = meta_for_heatmap,
             silent=F,
             clustering_distance_rows=clustering_distance_rows,
             clustering_distance_cols=clustering_distance_cols,
             clustering_method='ward.D2',
             gaps_col = breaks,
             fontsize_row=10,
             fontsize_col=10,
             show_rownames = T,
             show_colnames = T,
             border_color=NA,
             width=20, height = 30,
             filename = filename)
p
