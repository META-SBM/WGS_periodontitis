library(phyloseq)
library(ggplot2)
library(microbiome)
library(ggpubr)
library(grid)
#install.packages("PNWColors")
library(PNWColors)
library(gridExtra)

# function
plot_alpha_div <- function(
    ps_obj,
    group = 'batch',
    color = 'batch',
    measure = 'Shannon'
) {
  # Generate boxplot of alpha diversity without transformation
  p <- plot_richness(ps_obj, x = group, color = color, measures = measure) +
    geom_boxplot() 
  
  return(p)
} 
# Read phyloseq and define the parameters

ps <-  readRDS('/home/ignatsonets/ps_wgs.rds')
#ps <-  readRDS('/home/ignatsonets/ps_mags_filt.rds')
ps@sam_data$BMI_group <- factor(ps@sam_data$BMI_group, levels=c('normal', 'overweight', 'obesity'))
ps_obj <- ps
level <- 'Species'
det <- 1
prev<- 10/100
measure = 'Shannon'
method <- 'wilcox.test'
# Prepare the phyloseq object by pruning taxa
taxas <- core_members(ps_obj, detection = det, prevalence = prev, include.lowest = F)
ps_obj <- prune_taxa(taxas, ps_obj)
taxa_num <- length(colnames(ps_obj@otu_table))
ps_obj <- speedyseq::tax_glom(ps_obj,taxrank = level , NArm=T )
# Define the columns to keep for comparisons
#cols_to_keep <- c('Age_group', 'BMI_group','Sex','Arterial_hypertension','Cardiovascular_diseases','Rheumatoid_arthritis','Diabetes_mellitus','Antibiotics_treatment_in_last_6_months','Smoking_status','Stage_of_severity_of_periodontitis','Treatment_of_periodontitis','Patient_category')
cols_to_keep <- c('Age_group', 'BMI_group','Sex','Antibiotics_treatment_in_last_6_months',
                  'Stage_of_severity_of_periodontitis','Treatment_of_periodontitis')
# Generate alpha diversity plots for each column in cols_to_keep
plot_list <- list()
size=16
# Define the color palette
palette <- sample(pnw_palette("Sunset2", 8, type = "continuous"))

for (col in cols_to_keep) {
  # Dynamically generate comparison groups
  unique_vals <- unique(ps_obj@sam_data[[col]])
  if (length(unique_vals) > 1) {
    my_comparisons <- lapply(combn(unique_vals, 2, simplify = FALSE), as.vector)
    
    # Generate the alpha diversity plot for the current column
    p <- plot_alpha_div(ps_obj, group = col, color = col, measure = measure) +
      geom_violin(trim=F, alpha=0.1) +
      geom_boxplot(width=0.5, alpha=0.75, position=position_dodge(0.9)) +
      geom_jitter(size=1.5, alpha=0.5, position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) +
      scale_color_manual(values=palette)+
      stat_compare_means(comparisons = my_comparisons, method = method, label = "p.signif", size=5) +
      theme_bw(base_size=16)+
      theme(
        plot.title = element_text(color = "black", size=size),
        axis.text.y = element_text(color = "black", size = size),
        axis.text.x = element_text(angle=45, hjust=1,size=size,color = 'black'),
        legend.position = "none",
        axis.title.y  = element_text(color = "black", size = size-3,angle=90),
        axis.title.x  = element_text(color = "black", size = size-3),
        #legend.key.size = unit(0.5, 'cm'),
        text = element_text(size = size,colour ='black'),
        #plot.margin=unit(c(1,1,1,1), "cm")
        )
    # Add the plot to the list
    plot_list[[col]] <- p
  }
}

# Convert all plots to grobs
grob_list <- lapply(plot_list, ggplotGrob)

# Find maximum height among grobs
heights_list <- lapply(grob_list, function(g) g$heights)
max_height <- Reduce(unit.pmax, heights_list)

# Set all grobs to same height
grob_list_equal <- lapply(grob_list, function(g) {
  g$heights <- max_height
  g
})


# Display all plots (adjust ncol as needed)
#p111 <- do.call(grid.arrange, c(plot_list, ncol = 3))
p111 <- do.call(grid.arrange, c(grob_list_equal, ncol=3))
ggsave(
  filename = paste0("/home/ignatsonets/shannons_new_WGS_combined_plot_v4_fixed_order.png"),
  plot = p111,
  width = 18,
  height = 14,
  units = "in"
)
p111
