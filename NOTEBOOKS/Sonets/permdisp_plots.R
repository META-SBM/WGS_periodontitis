library(vegan)
library(phyloseq)

# Reqd phyloseq
ps <-  readRDS('/home/ignatsonets/ps_mags_filt.rds')
#ps <-  readRDS('/home/ignatsonets/ps_humann_joint.rds')
#ps <-  readRDS('/home/ignatsonets/ps_wgs.rds')
ps <- prune_samples(sample_sums(ps) > 0, ps)
otu_table <- as.data.frame(otu_table(ps))  # Samples as rows
metadata <- as.data.frame(sample_data(ps))

cols_to_keep <- c('Age_group','BMI_group','Sex','Cardiovascular_diseases','Rheumatoid_arthritis',
                  'Antibiotics_treatment_in_last_6_months','Smoking_status',
                  'Stage_of_severity_of_periodontitis','Treatment_of_periodontitis')
for (var in cols_to_keep) {
  if (!var %in% colnames(metadata)) next
  non_na <- !is.na(metadata[[var]])
  meta_sub <- metadata[non_na, , drop = FALSE]
  otu_sub <- otu_table[non_na, ]
  
  valid_samples <- rowSums(otu_sub) > 0
  otu_sub <- otu_sub[valid_samples, ]
  meta_sub <- meta_sub[valid_samples, , drop = FALSE]
  
  group <- droplevels(factor(meta_sub[[var]]))
  if (nlevels(group) < 2 || min(table(group)) < 2) next
  
  dist_matrix <- vegdist(otu_sub, method = "bray")
  dispersion <- betadisper(dist_matrix, group = group)
  perm_test <- permutest(dispersion, permutations = 999)
  cat("\nPERMDISP (method='bray') results for", var, ":\n")
  print(perm_test)

  # Plot with adjusted margins and labels
  png(paste0("/home/ignatsonets/PERMDISP_mags_", var, ".png"), width = 1200, height = 600)
  par(mfrow = c(1, 2), mar = c(8, 5, 4, 2))  # Larger bottom margin (8)
  
  # Boxplot: Rotated (90°) and smaller group labels
  boxplot(dispersion, 
          main = paste("Dispersion by", var),
          ylab = "Distance to Centroid",
          col = rainbow(nlevels(group)),
          las = 2,          # Rotate group labels 90°
          cex.axis = 0.8,   # Smaller group labels
          xlab = "",        # Remove redundant x-axis title
          xaxt = "n")       # Suppress default x-axis
  
  # Add rotated group labels with adjusted spacing
  text(
    x = 1:nlevels(group),
    y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),  # Position below plot
    labels = levels(group),
    srt = 90,           # 90° rotation
    adj = 1,            # Right-align labels
    xpd = TRUE,         # Allow drawing outside plot area
    cex = 0.8           # Smaller label size
  )
  
  # Ordination plot (unchanged)
  plot(dispersion, hull = FALSE, ellipse = TRUE,
       main = paste("Multivariate Dispersion:", var),
       col = 1:nlevels(group), pch = 16)
  legend("topright", legend = levels(group),
         col = 1:nlevels(group), pch = 16, cex = 0.8)
  
  dev.off()
  
  cat("\nResults for", var, ":\n")
  print(perm_test)
}
