# upd 03/06/25
# libs
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(vegan)
library(dplyr)
library(gtable)
library(grid)
library(lemon)
library(PNWColors)  # For color palettes
library(microbiome)
library(patchwork)
# Function for beta-diversity plots
create_ordination_plots <- function(ps_obj, method = "PCoA", distance_method = "wunifrac", group, size = 10, palette) {
  # Calculate distance matrix
  dist <- phyloseq::distance(ps_obj, method = distance_method)
  # Perform ordination
  ordination <- ordinate(ps_obj, method = method, distance = dist)
  P <- cbind(as.data.frame(ordination$vectors), as.data.frame(as.matrix(sample_data(ps_obj))))
  
  # Ensure axes are numeric
  P$Axis.1 <- as.numeric(as.character(P$Axis.1))
  P$Axis.2 <- as.numeric(as.character(P$Axis.2))
  
  # Variance explained
  variance_explained <- 100 * ordination$values$Relative_eig[1:2]
  
  # Means per group
  means <- P %>%
    group_by(!!rlang::sym(group)) %>%
    summarise(
      mean_Axis1 = mean(Axis.1, na.rm = TRUE),
      mean_Axis2 = mean(Axis.2, na.rm = TRUE)
    )
  
  # Merge means
  P <- P %>%
    dplyr::left_join(means, by = c(group))
  
  # ANOSIM
  anosim_res <- anosim(dist, P[[group]])
  anosim_text <- paste("ANOSIM R:", round(anosim_res$statistic, 3), 
                       "p-value:", round(anosim_res$signif, 3))
  
  # Plot ordination
  pl <- ggplot(P, aes_string(x = "Axis.1", y = "Axis.2", color = group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse() +
    geom_segment(aes(x = Axis.1, y = Axis.2, xend = mean_Axis1, yend = mean_Axis2), alpha = 0.3) +
    geom_label(data = means, aes(x = mean_Axis1, y = mean_Axis2, label = !!rlang::sym(group)),
               fill = "white", color = "black", fontface = "bold", size = 7) +
    geom_text(x = Inf, y = Inf, label = anosim_text, hjust = 1.1, vjust = 1.1, size = size-10, color = "black") +
    scale_color_manual(values = palette) +
    xlab(paste0('PCoA component 1 [', round(variance_explained[1], 1), '%]')) +
    ylab(paste0('PCoA component 2 [', round(variance_explained[2], 1), '%]')) +
    theme(
      legend.position = "right",
      text = element_text(size = size, color = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),
      axis.text.x = element_text(size = size+10, color='black'),
      axis.text.y = element_text(size = size+10, color='black'),
      axis.title.x = element_text(size=25),
      axis.title.y = element_text(size=25),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      legend.key.size = unit(2.0, "cm"),
      legend.text=element_text(size=20),
      legend.title = element_text(size=20)
    )
  legend <- g_legend(pl)
  pl <- pl + theme(legend.position="none")
  
  # Pairwise comparisons for Wilcoxon tests
  unique_vals <- unique(P[[group]])
  my_comparisons <- combn(unique_vals, 2, simplify = FALSE) %>% lapply(as.vector)
  
  # Density plots for Axis.1
  x_dens <- ggplot(P, aes_string(x = "Axis.1", y = group, color = group)) + # Add fill=group if necessary
    geom_violin(trim = FALSE, alpha = 0.1) +
    #geom_boxplot(width = 0.5, alpha = 0.75, position = position_dodge(0.9))
    geom_boxplot(width = 0.5, alpha = 0.75) + # Add boxplot inside violin if you want
    #geom_jitter(size = 2.5, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9))
    geom_jitter(size = 2.5, alpha = 0.5) +
    scale_color_manual(values = palette) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black", size = size+10),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size+10, color = 'black'),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    ) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                       label = "p.signif", y.position = max(P$Axis.1) * 1.1, size = size-8)
  
  # Density plots for Axis.2
  y_dens <- ggplot(P, aes_string(x = "Axis.2", y = group, color = group)) + coord_flip() + # Add fill=group if necessary
    geom_violin(trim = FALSE, alpha = 0.1) +
    #geom_boxplot(width = 0.5, alpha = 0.75, position = position_dodge(0.9))
    geom_boxplot(width = 0.5, alpha = 0.75) + # Add boxplot inside violin if you want
    #geom_jitter(size = 2.5, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9))
    geom_jitter(size = 2.5, alpha = 0.5) +
    scale_color_manual(values = palette) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black", size = size+10),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size+10, color = 'black'),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    ) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                       label = "p.signif", y.position = max(P$Axis.2) * 1.1, size = size-8)
  
  # theme
  theme_new <- theme_classic() + 
    theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "cm"), text=element_text(size=20))
  
  theme_set(theme_new)
  
  # Convert to grobs
  gA <- ggplotGrob(x_dens)
  gB <- ggplotGrob(pl)
  gD <- ggplotGrob(y_dens)
  gL <- legend
  
  # Set widths
  xWidth = unit.pmax(gB$widths)
  yHeight = unit.pmax(gB$heights)
  gA$widths<- xWidth
  gA$heights[2:5] <- yHeight
  gB$widths <- xWidth
  gB$heights <- yHeight
  gD$widths[2:5]<- xWidth
  gD$heights<- yHeight
  
  # Arrange plots
  grid.arrange(gD, gB, gL, gA, ncol=2, nrow=2, widths=c(7, 8), heights=c(8, 7))
}

# Function for alpha-diversity plot
plot_alpha_div <- function(ps_obj, group = 'batch', measure = 'Shannon', palette) {
  # Extract data frame for richness
  richness_df <- data.frame(
    sample_data(ps_obj)[[group]],
    measure = estimate_richness(ps_obj, measures = measure)[[measure]],
    sample_names = sample_names(ps_obj)
  )
  colnames(richness_df)[1] <- "group"
  
  # Convert group to factor with levels matching palette names
  richness_df$group <- factor(richness_df$group, levels=names(palette))
  
  # Define pairwise comparisons
  unique_groups <- levels(richness_df$group)
  comparisons <- combn(unique_groups, 2, simplify=FALSE)
  
  p <- ggplot(richness_df, aes(x = group, y = measure, color = group, fill=group)) + # Add fill for boxplot if necessary
    geom_boxplot(alpha=0.5,) +
    geom_jitter(width=0.2, size=2,) +
    scale_color_manual(values=palette) +
    scale_fill_manual(values=palette) +
    theme_minimal() + 
    ylab('Shannon diversity index')+
    theme(
      axis.text.x = element_text(angle=45, hjust=1,size=20),
      axis.title.y = element_text(size=20),
      legend.position = "none"
    ) +
   
    # Add significance comparisons
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.signif", size=size-5, tip.length=0.01)
  return(p)
}

# Main script
# read phyloseq and specify parameters
ps <- readRDS('/home/ignatsonets/ps_wgs.rds')
#ps <- readRDS('/home/ignatsonets/ps_humann_joint.rds')
#ps <- readRDS('/home/ignatsonets/ps_mags_filt.rds')

ps_obj <- ps
det <- 1
prev <- 10/100
taxas <- core_members(ps_obj, detection=det, prevalence=prev, include.lowest=F)
ps_obj <- prune_taxa(taxas, ps_obj)

# Define parameters
measure <- "Shannon"
trans <- 'compositional'
distance_method <- 'bray'
size <- 16
output_dir <- "/home/ignatsonets/AB_plots_v7_big_font_newgrid/"
#cols_to_keep <- c('Age_group', 'BMI_group', 'Sex') 
cols_to_keep <- c('Age_group', 'BMI_group','Sex','Arterial_hypertension',
                  'Cardiovascular_diseases',
                  'Rheumatoid_arthritis','Diabetes_mellitus','IBD',
                  'Antibiotics_treatment_in_last_6_months',
                 'Smoking_status','Treatment_of_periodontitis',
                 'Patient_category',
                 'Stage_of_severity_of_periodontitis')# your variables
# Transform data
ps_obj <- microbiome::transform(ps_obj, trans)

# Define consistent alpha plot theme
alpha_theme <- function(size) {
  theme(
    axis.text.y = element_text(color = "black", size = size+5),
    axis.text.x = element_text(angle=90, hjust=1, size = size+5, color='black'),
    legend.position = "none",
    axis.title.y = element_text(color='black', size=25, angle=90),
    axis.title.x = element_text(color='black', size=size+5),
    panel.border = element_rect(fill=NA, color='black'),
    panel.background = element_rect(fill='white', color='black'),
    panel.grid.major = element_line(size=0.1, linetype='solid', color='gray'),
    panel.grid.minor = element_line(size=0.1, linetype='solid', color='gray'),
    text = element_text(size=size, color='black')
  )
}

# Loop through variables
for (group_var in cols_to_keep) {
  tryCatch({
    # Generate palette with named levels matching groups
    group_levels <- unique(as.character(sample_data(ps_obj)[[group_var]]))
    group_levels <- sort(group_levels)
    palette_colors <- PNWColors::pnw_palette('Bay', length(group_levels))
    names(palette_colors) <- group_levels
    
    # Create beta plot
    p_beta <- create_ordination_plots(
      ps_obj,
      method = "PCoA",
      distance_method = distance_method,
      group = group_var,
      size = size,
      palette = palette_colors
    )
    
    # Create alpha plot
    p_alpha <- plot_alpha_div(
      ps_obj,
      group = group_var,
      measure = measure,
      palette = palette_colors
    ) + alpha_theme(size)
    
    # Adjust plot margins
    p_alpha <- p_alpha + theme(plot.margin = unit(c(0, 1.5, 0, 0), "cm"))
    
    # Combine plots side by side with tags
    P <- (p_alpha | p_beta) + 
      plot_layout(widths = c(0.5, 1)) + 
      plot_annotation(tag_levels='A', title = group_var, 
                      theme = theme(plot.title = element_text(size = 28, face = "bold"))) &
      theme(
        plot.tag = element_text(size=28, face='bold'),
        plot.tag.position = c(0.015, 0.97)
      )
    
    # Save plot
    ggsave(
      filename = paste0(output_dir, group_var, "_v7_MAGs_combined_plot_big_font_newgrid.png"),
      plot = P,
      width = 22,
      height = 18,
      units = "in"
    )
    
    message("Successfully created plot for: ", group_var)
    
  }, error=function(e){
    warning("Failed for ", group_var, ": ", e$message)
  })
}
