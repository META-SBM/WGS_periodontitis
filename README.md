# WGS_porphyromonas

Repo for WGS metagenomics analysis of periodontal pockets of patients with mild and severe periodontitis.

**UPD 26/08/2025**: README is ready, code is more polished and usable now.

Available code are inside of 2 subfolders:

- **Galeeva/**
    + deseq.Rmd: code for DESEq2 analysis of WGS and MAGs datasets;
    + final_permanova_v2.Rmd: code for PERMANOVA analysis of all 3 datasets;
- **Sonets/**
    + MPAR_readcount_barplot.ipynb: code for Figure S1;
    + alpha_beta_combined_plots.R: code for creating plots of alpha- and beta-diversity analysis;
    + alpha_div_grid_plots.R: code for alpha-diversity plots;
    + chi_square_plots.R: code for chi-square analysis(veryfying PERMANOVA results);
    + create.ps.R: create phyloseq objects of WGS/MAGs/HUMAnN datasets from raw_data;
    + merge_instrain_cov.py: creating 1 big table from InStrain results for MAGs phyloseq object;
    + metadata_plot.R: code for Figure 1;
    + metaphlan2phyloseq: create table from MetaPhLAn results for WGS phyloseq object;
    + multiple_plots.py: code for Figs S6/7/8 (combined alpha-/beta-diversity plots for different covariates for 3 datasets. Modify this code to suit your needs);
    + permdips_plots.R: code for PERMDISP analysis (and plots if desired);
    + top_taxa.R: code for Fig 2 and Figs S2,3 (top-15 taxons on different levels);
    + wgs_heatmap_w_metadata: code for heatmap of top prevalent species and MAGs with clusterization and metadata(cheage ps for MAGs heatmap).

If you have any questions and/or troubles, please contact ignatsonets@gmail.com for advice, troubleshhoting and help. We hope that this repo will help with your future research!
       
