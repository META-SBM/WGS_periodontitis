library(microViz)
library(phyloseq)
library(ggplot2)
library(patchwork)

# Read phyloseq
ps <- readRDS('/home/ignatsonets/ps_wgs.rds')
ps_obj <- ps

# Define severity groups and taxonomic ranks
group_col <- "Stage_of_severity_of_periodontitis"
ranks <- c("Family", "Genus", "Species")

# Generate and save plots for each taxonomic rank
for (rank in ranks) {
  # Create plot for this rank
  p_list <- comp_barplot(
    ps_obj,
    tax_level = rank,
    n_taxa = 15,
    group_by = group_col,
    merge_other = TRUE,
    sample_order = "aitchison",
    bar_outline_colour = NA,
    bar_width = 0.9
  )
  
  # Create a unified theme for all elements
  unified_theme <- theme(
    # Grid lines
    panel.grid.major = element_line(color = "gray90", size = 0.1),
    panel.grid.minor = element_blank(),
    
    # Titles and text
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    
    # Legend
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  )
  
  # Apply unified theme to each plot
  p_list <- lapply(p_list, function(p) {
    p + 
      unified_theme
  })
  
  # Add grid to each plot in the list
  p_list <- lapply(p_list, function(p) {
    p + theme(
      panel.grid.major = element_line(color = "gray90", size = 0.1),
      panel.grid.minor = element_blank()
    )
  })
  
  # Combine the group plots horizontally with shared legend
  combined_plot <- wrap_plots(p_list, ncol = length(p_list), guides = "collect") + 
    plot_annotation(title = paste("Top taxons at", rank, "level")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  combined_plot <-  combined_plot & coord_flip()
  # Save the plot
  ggsave(paste0("top_15_", tolower(rank), "_by_severity.png"), 
         combined_plot, 
         width = 18, height = 10, dpi = 300)
  
  # Print the plot
  print(combined_plot)
}
