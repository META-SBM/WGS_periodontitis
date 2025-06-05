library(vegan)
library(phyloseq)
# Read phyloseq
ps <-  readRDS('/home/ignatsonets/ps_wgs.rds')
ps_obj <- ps
ps_meta <- as(sample_data(ps_obj), 'data.frame')
create_plot_rel <- function(cols, palette,ps_meta) {
  dataex <- ps_meta %>%
    dplyr::group_by_at(cols) %>%
    dplyr::summarise(N = n()) %>%
    dplyr::mutate(rel_count = N / sum(N)) 
  
  # Perform chi-squared test
  chisq_result <- chisq.test(ps_meta[[cols[1]]], ps_meta[[cols[2]]])
  chisq_pvalue <- format(chisq_result$p.value, digits = 3)
  
  ggplot(dataex, aes(x = !!sym(cols[[1]]), y = rel_count, fill = !!sym(cols[[2]]))) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text(aes(label = scales::percent(rel_count, accuracy = 1)), position = position_stack(vjust = 0.5), size = 6)  +
    labs(title = paste("Chi-squared p-value:", chisq_pvalue ), fill = cols[[2]]) +
    scale_fill_manual(values = palette) +
    theme_linedraw()
}
# Create a vector of colors
colors <- c("#34b9ed", "#1d8dcf", "#024cc1", "#a100bf", "#cc7daa", "#bf4380", "#ff594c", "#f7b126", "#f0e74c", 
            "#E05A3F", "#9467BD", "#FFC266", "#4E79A7", "#F28E2B", "#76B7B2", "#D62728", "#7F7F7F")
variables <-  c('Age_group', 'BMI_group','Sex','Arterial_hypertension','Cardiovascular_diseases','Rheumatoid_arthritis','Diabetes_mellitus','IBD','Antibiotics_treatment_in_last_6_months','Smoking_status','Treatment_of_periodontitis','Patient_category','Stage_of_severity_of_periodontitis')
# Loop over variables and plot each with "phenotype" as the x-axis variable
plots <- lapply(variables, function(var) {
  create_plot_rel(c('Age_group', var), colors, ps_meta)
})

# Display one of the plots, for example the first
print(plots)
