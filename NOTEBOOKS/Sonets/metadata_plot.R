library(tidyverse)
library(phyloseq)

# Keep all severity groups
ps <-  readRDS('/home/ignatsonets/ps_wgs.rds')
ps_obj <- ps
df <- as.data.frame(as.matrix(ps_obj@sam_data))

# Include severity in the variables
cols_to_keep <- c('Stage_of_severity_of_periodontitis',
                 'Age_group', 'BMI_group','Sex','Arterial_hypertension',
                 'Cardiovascular_diseases','Rheumatoid_arthritis','Diabetes_mellitus',
                 'IBD','Antibiotics_treatment_in_last_6_months','Smoking_status',
                 'Treatment_of_periodontitis','Patient_category')

df_f <- select(df, all_of(cols_to_keep))

# Create long format while preserving severity
df_long <- df_f %>%
  pivot_longer(cols = -Stage_of_severity_of_periodontitis,
               names_to = "Column_Name",
               values_to = "Value")

# Calculate counts and percentages WITHIN each category value
df_summary <- df_long %>%
  group_by(Column_Name, Value, Stage_of_severity_of_periodontitis) %>%
  tally() %>%
  group_by(Column_Name, Value) %>%  # Group by category value for percentages
  mutate(
    total = sum(n),
    percent = n/total * 100
  ) %>%
  ungroup()

size <- 20  # Doubled base text size
pdf(paste0('/home/ignatsonets/metadata_distribution.pdf'), width = 35, height =60) 
p <- ggplot(df_summary, aes(x = Value, y = n, fill = Stage_of_severity_of_periodontitis)) +
  geom_col() +
  geom_text(aes(label = paste0(n, "\n(", round(percent, 1), "%)")),
            position = position_stack(vjust = 0.5),
            size = size-4, color = "black") +  
  facet_wrap(~ Column_Name, scales = "free", ncol = 4) +
  labs(title = "Case Composition by Category Value",
       x = "Category Values",
       y = "Total Count",
       fill = "Severity") +
  theme(
    axis.text.y = element_text(size = size, color = "black"),
    axis.text.x = element_text(size = size, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = size + 4, color = "black"),  # Increased title size
    legend.text = element_text(size = size),
    legend.title = element_text(size = size + 2),  # Larger legend title
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    strip.text = element_text(size = size + 5, color = "black"),  # Larger facet titles
    strip.background = element_rect(fill = "gray80", color = "black"),
    plot.title = element_text(size = size + 10, face = "bold")  # Main title size
  ) +
  scale_fill_manual(values = c("mild" = "#2c7bb6", "severe" = "#d7191c")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
p
dev.off()
ggsave(p, path='/home/ignatsonets/', filename='metadata_distribution.jpeg',device = 'jpeg',units="in", width=35, height=60, dpi=300 ,limitsize = FALSE)
