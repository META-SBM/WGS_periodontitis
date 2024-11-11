# phyloseq tutorial
# 27-Oct-2024
# load libs
library("phyloseq")
library("ggplot2")      
library("readxl")      
library("dplyr")        
library("tibble") 
library("vegan")
library("scales")
library("microbiome")

# tutorial see here: https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

# load data
otu_count_df<- read.csv("/home/ignatsonets/metaphlan_count_df.csv", sep='\t')
#otu_count_df <- slice(otu_count_df, -1)
taxa_df<- read.csv("/home/ignatsonets/metaphlan_taxa_df.csv", sep='\t')
samples_df <- read.csv("/home/ignatsonets/metaphlan_full_metadata_orig.csv")
samples_df$Стадия.тяжести.пародонтита <- as.factor(samples_df$Стадия.тяжести.пародонтита)
#taxa_df <- slice(taxa_df, -1)
# remove underscore to reach consistency in sample IDs
colnames(otu_count_df)<-gsub("_","",colnames(otu_count_df))

# NA filling for empty values in samples metadata df
samples_df <- samples_df %>% mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

# defining row names
otu_count_df <- otu_count_df %>%
  tibble::column_to_rownames("OTU") 
samples_df <- samples_df %>% 
  tibble::column_to_rownames("ID") 
taxa_df <- taxa_df %>% 
  tibble::column_to_rownames("OTU")
otu_count_df <- as.matrix(otu_count_df)
taxa_df <- as.matrix(taxa_df)

# creating phyloseq
OTU = otu_table(otu_count_df, taxa_are_rows = TRUE)
TAX = tax_table(taxa_df)
SAMPLES = sample_data(samples_df)
ps_phyloseq <- phyloseq(OTU, TAX, SAMPLES)

# playing with phyloseq
sample_names(ps_phyloseq)
rank_names(ps_phyloseq)
sample_variables(ps_phyloseq)

# filtering options, for example
ps_phyloseq_sample_subset <- subset_samples(ps_phyloseq, Стадия.тяжести.пародонтита  =="средняя")
ps_phyloseq_taxa_subset <- subset_taxa(ps_phyloseq, Phylum %in% c("Proteobacteria", "Firmicutes"))

# reads normalization
total = median(sample_sums(ps_phyloseq))
standf = function(x, t=total) round(t * (x / sum(x)))
ps_phyloseq_norm = transform_sample_counts(ps_phyloseq, standf)

# barplots just 4 fun
plot1 <- plot_bar(ps_phyloseq, fill = "Class")
plot1
plot_norm <- plot_bar(ps_phyloseq_norm, fill = "Class")
plot_norm

# heatmaps with filtering by abundance
ps_phyloseq_abund <- filter_taxa(ps_phyloseq_norm, function(x) sum(x > total*0.05) > 0, TRUE)
ps_phyloseq_abund
heatmap_1 <- plot_heatmap(ps_phyloseq, method = "NMDS", distance = "bray", taxa.label = "Class", taxa.order = "Class")
heatmap_1
heatmap_2 <- plot_heatmap(ps_phyloseq_abund, method = "NMDS", distance = "bray",)
heatmap_2
heatmap_3 <- plot_heatmap(ps_phyloseq_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Phylum", taxa.order = "Phylum", 
             trans=NULL, low="beige", high="red", na.value="beige")
heatmap_3

# alpha-diversity 
alpha1 <- plot_richness(ps_phyloseq_abund, measures=c("Chao1", "Shannon"))
alpha1

#ordination
ps_phyloseq.ord <- ordinate(ps_phyloseq_abund, "NMDS", "bray")
ord1 <- plot_ordination(ps_phyloseq_abund, ps_phyloseq.ord, type="taxa", color="Family", shape= "Genus", 
                title="OTUs")
ord1
ord2 <-   plot_ordination(ps_phyloseq_abund, ps_phyloseq.ord, type="taxa", color="Genus", 
                          title="OTUs", label="Genus") + 
  facet_wrap(~Family, 3)
ord2

# network
net1 <- plot_net(ps_phyloseq_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
                 maxdist = 0.8, color="Family", point_label="Genus") 
net1

# trying to create with taxa_are_rows=F
OTU2 = otu_table(otu_count_df, taxa_are_rows = FALSE)
TAX = tax_table(taxa_df)
samples = sample_data(samples_df)
ps_phyloseq2 <- phyloseq(OTU2, TAX, SAMPLES)
# fails as intended, I suppose
# Error in validObject(.Object) : invalid class “phyloseq” object: 
# Component taxa/OTU names do not match.
# Taxa indices are critical to analysis.

#taxa agglomeration
ps_genus <- tax_glom(ps_phyloseq, taxrank="Genus")
ps_species <- tax_glom(ps_phyloseq, taxrank="Species")


# to replace OTU table: assign-otu_table
# to replace tax table: assign-tax_table
# to replace sample table: assign-sample_data

# 7-11 Nov 2024: adonis, compositional analysis & clr-transformation
# compositional analysis => Bray-Curtis distance
# clr-transformation => Euclidean
# do we need to agglomerate?? try anyway

# trying from this tutorial: https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# and this: https://microbiome.github.io/tutorials/PERMANOVA.html
# may be useful: https://chrischizinski.github.io/rstats/adonis/
# also in different corners of Internet there are transposing of matrix, dunno why
# genus
ps_genus <- tax_glom(ps_phyloseq, taxrank="Genus")
ps_phyloseq_g_bray.ord <- ordinate(ps_genus, method = "bray")
ps_phyloseq_g_bray.matrix <- phyloseq::distance(ps_genus, method = "bray")
# omitting NAs for now, add more sample data later(same approach for species data and all data)
adonis_genus <- adonis2(ps_phyloseq_g_bray.matrix ~ Стадия.тяжести.пародонтита, data = samples_df, na.action = 'na.omit')

# species
ps_species <- tax_glom(ps_phyloseq, taxrank="Species")
ps_phyloseq_s_bray.ord <- ordinate(ps_species, method = "bray")
ps_phyloseq_s_bray.matrix <- phyloseq::distance(ps_species, method = "bray")
adonis_species <- adonis2(ps_phyloseq_s_bray.matrix ~ Стадия.тяжести.пародонтита, data = samples_df, na.action = 'na.omit')

# all data, no agglomeration
ps_phyloseq_bray.ord <- ordinate(ps_phyloseq, method = "bray")
ps_phyloseq_bray.matrix <- phyloseq::distance(ps_phyloseq, method = "bray")
adonis_all <- adonis2(ps_phyloseq_bray.matrix ~ Стадия.тяжести.пародонтита, data = samples_df, na.action = 'na.omit')
# it works, I suppose, but not significant

# barplots for all data
coef <- coef(adonis_all)["Total",]
# coeff NULL, wtf
# not significant?
# code below should work, but due to NULL, stop here.
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

# CLR to see if it even works
# tutorial here: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
# also may be useful: https://microbiome.github.io/OMA/docs/devel/pages/21_microbiome_community.html
# probably we need microbiome package, install if necessary
ps_porph_clr <- microbiome::transform(ps_phyloseq, "clr", method = "euclidean")
ps_clr <- adonis2(ps_porph_clr ~ phyloseq::sample_data(ps_phyloseq)$Стадия.тяжести.пародонтита, na.action='na.omit')
# here na.omit doesn't work, strange
ps_porph_g_clr <- microbiome::transform(ps_genus, "clr", method = "euclidean")
ps_porph_s_clr <- microbiome::transform(ps_species, "clr", method = "euclidean")

# todo: get advices, get codereview