# Bradley June 2025

# Packages
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)

# Options
scale=F
indir="results/"
modulename="schogginsISG_fixed"
add_meta_file="input/eQTL_interactions.tsv" # Copied from what Tobi used here - /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/interaction_files/2024_12_27-multi_tissue/eQTL_interactions.tsv
umap.palette.df = read.csv("/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/IBDVerse-sc-eQTL-code/data/palette.csv")
annot.mapping <- read_csv(paste0('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/IBDVerse-sc-eQTL-code/data/all_IBDverse_annotation_mastersheet.csv'), n_max = 103) %>% 
    dplyr::rename(label_machine = leiden,
                label_new = JAMBOREE_ANNOTATION,
                category_new = Category) %>% 
    dplyr::select(label_machine, label_new, category_new) %>%
    mutate(
        label_new = str_replace_all(label_new, '_', ' ')
    ) %>% 
    rename(
        predicted_labels = label_machine
    )


# Define other variables
out=paste0("plotout/", modulename, "/")
if(scale){
    indir = paste0(indir, "scaled-")
    out = paste0(out, "scaled-")
}
allres=paste0(indir, modulename, "-module_scores.txt.gz")
tires=paste0(indir, modulename, "-TI-module_scores.txt.gz")
umap.category.palette = deframe(dplyr::select(umap.palette.df, category, category_color))


#######################
# What is the overall expression of this signature like?
#######################
# Load
ms = fread(allres) %>% 
    mutate( # Couple adjustment
        predicted_labels = gsub("_ct", "", predicted_labels),
        category = gsub("_ct", "", predicted_category),
        disease_status = factor(disease_status, levels = c("Healthy", "CD"))
    ) %>% # Add the proper names
    left_join(select(annot.mapping, c('predicted_labels', 'label_new')))

# Merge with other meta
add_meta = read.delim(add_meta_file) %>% 
    select(Genotyping_ID, ses_cd_binned, ses_inflamed) %>% 
    mutate(
        ses_cd_binned = case_when( # ses_cd_binned is the ses scored binned into 0-3, 4-6, 6-9.
            ses_cd_binned == 0 ~ "mild",
            ses_cd_binned == 1 ~ "moderate",
            TRUE ~ "severe"
        ), 
        ses_inflamed = case_when( # ses_inflamed is a binary score of >=3 is 1 (inflamed), <3 is uninflamed
            ses_inflamed == 1 ~ "inflamed",
            TRUE ~ "uninflamed"
        ) 
    )

ms = ms %>% 
    left_join(add_meta, "Genotyping_ID")


# Make sure outdir exists
if(!file.exists(out)){
    dir.create(out)
}

# Group this: 
# Median scores per sample per category
msamp = ms %>% 
    select(-predicted_labels) %>% 
    group_by(category, sanger_sample_id) %>% 
    mutate(
        median_ISG = median(!!sym(modulename))
    ) %>% 
    select(-c(!!modulename, cell)) %>% 
    distinct()


# Plot
ggplot(msamp, aes(x = category, y = median_ISG, fill=category)) +
  geom_violin(trim = FALSE, color = "black") +
  theme_classic() +
  geom_text(
        data = msamp %>%
            group_by(category) %>%
            summarise(
                median_val = median(median_ISG),
                max_median= max(median_ISG)  # Calculate max for positioning text
            ),
        aes(
            x = category,
            y = max_median*1.2,  # Position text slightly above the violins
            label = sprintf("%.2f", median_val)
        ),
        inherit.aes = FALSE,
        hjust = 0.5,
        color = "black",
        size = 4
    ) + 
  scale_fill_manual(values=umap.category.palette) + 
  labs(title = "ISG expression across category",x="Major population") + 
  theme(legend.position="none")

ggsave(paste0(out, "ISG_per_per_category.pdf"),
         width = 10, height = 4,device = cairo_pdf)

# If we look at major lineages? 
msampl = ms %>% 
    mutate(
        lineage = case_when(
            category == "Mesenchymal" ~ "Mesenchymal",
            category %in% c("Colonocyte", "Stem", "Secretory", "Enterocyte") ~ "Epithelial",
            TRUE ~ "Immune"
        )
    ) %>% 
    group_by(lineage, sanger_sample_id) %>% 
    mutate(
        median_ISG_lineage = median(!!sym(modulename))
    )


ggplot(msampl, aes(x = lineage, y = median_ISG_lineage, fill=lineage)) +
  geom_violin(trim = FALSE, color = "black") +
  theme_classic() +
  geom_text(
        data = msampl %>%
            group_by(lineage) %>%
            summarise(
                median_val = median(median_ISG_lineage),
                max_median= max(median_ISG_lineage)  # Calculate max for positioning text
            ),
        aes(
            x = lineage,
            y = max_median*1.2,  # Position text slightly above the violins
            label = sprintf("%.2f", median_val)
        ),
        inherit.aes = FALSE,
        hjust = 0.5,
        color = "black",
        size = 4
    ) + 
  labs(title = "ISG expression across lineage",x="Lineage") + 
  theme(legend.position="none")

ggsave(paste0(out, "ISG_per_per_lineage.pdf"),
         width = 6, height = 4,device = cairo_pdf)

# Just in healthy? 
ggplot(filter(msamp, disease_status == "Healthy"), aes(x = category, y = median_ISG, fill=category)) +
  geom_violin(trim = FALSE, color = "black") +
  theme_classic() +
  geom_text(
        data = msamp %>%
            filter(
                disease_status == "Healthy"
            ) %>% 
            group_by(category) %>%
            summarise(
                median_val = median(median_ISG),
                max_median= max(median_ISG)  # Calculate max for positioning text
            ),
        aes(
            x = category,
            y = max_median*1.2,  # Position text slightly above the violins
            label = sprintf("%.2f", median_val)
        ),
        inherit.aes = FALSE,
        hjust = 0.5,
        color = "black",
        size = 4
    ) + 
  scale_fill_manual(values=umap.category.palette) + 
  labs(title = "ISG expression across category - Healthy only",x="Major population") + 
  theme(legend.position="none")

ggsave(paste0(out, "ISG_per_per_category_healthy_only.pdf"),
         width = 10, height = 4,device = cairo_pdf)

########################
# What if we look at the cell-type level (with same filters we use for eQTL mapping)
########################
ncellsdf = ms %>% 
    group_by(sanger_sample_id, label_new) %>% 
    summarise(count = n()) %>%
    mutate(enough_cells = count > 4) %>% # 5 or more
    filter(enough_cells) %>% 
    mutate(sample_cell = paste(sanger_sample_id, label_new, sep = "-"))

nsamplesdf = ncellsdf %>% 
    group_by(label_new) %>%
    summarise(count_nsamples = n()) %>%
    mutate(enough_samples = count_nsamples > 29) %>% # 30 or more
    filter(enough_samples)

ms.filt.group = ms %>% 
    mutate(sample_cell = paste(sanger_sample_id, label_new, sep = "-")) %>% 
    filter(
        sample_cell %in% unique(ncellsdf$sample_cell),
        label_new %in% nsamplesdf$label_new
    ) %>% 
    group_by(label_new, sanger_sample_id) %>%
    mutate(
        median_ISG = median(!!sym(modulename))
    ) %>% 
    select(-c(!!modulename, cell)) %>% 
    distinct() %>% 
    group_by(label_new) %>%
    mutate(median_median_ISG = median(median_ISG)) %>% 
    ungroup() %>%
    arrange(-median_ISG)

ordered_celltypes = ms.filt.group %>% 
    distinct(label_new, median_median_ISG) %>% 
    arrange(median_median_ISG) %>% 
    pull(label_new)

ms.filt.group = ms.filt.group %>%
    mutate(label_new = factor(label_new, levels = ordered_celltypes))

# Plot
ggplot(ms.filt.group, aes(x = label_new, y = median_ISG, fill=category)) +
  geom_violin(trim = FALSE, color = "black") +
  theme_classic() +
  scale_fill_manual(values=umap.category.palette) + 
  labs(title = "ISG expression across cell-types",x="Cell type") + 
  theme(
    legend.position="none", 
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(size=30, face="bold")
    )

ggsave(paste0(out, "ISG_per_per_celltype.pdf"),
         width = 30, height = 10,device = cairo_pdf)

#######################
# How does this vary by disease status?
#######################
# Use data just from the TI - so a seperate analysis
mstisamp = fread(tires) %>% 
    mutate(
        category = gsub("_ct", "", predicted_category),
        disease_status = factor(disease_status, levels = c("Healthy", "CD"))
    ) %>% 
    left_join(add_meta, "Genotyping_ID") %>%  
    mutate(
        disease_inflam = paste0(disease_status, "-", ses_inflamed),
        disease_inflam = factor(disease_inflam, levels=c("Healthy-uninflamed", "CD-uninflamed", "CD-inflamed"))
    ) %>% 
    select(-predicted_labels) %>% 
    group_by(category, sanger_sample_id) %>% 
    mutate(
        median_ISG = median(!!sym(modulename))
    ) %>% 
    select(-c(!!modulename, cell)) %>% 
    distinct()


groups = c("disease_status", "disease_inflam")
for(g in groups) {
    print(paste0("..Plotting - ", g))
    vals = distinct(as.data.frame(mstisamp[,g]))[,1]
    comparisons = combn(as.character(vals), 2, simplify = FALSE)

    p = ggplot(mstisamp, aes(x = .data[[g]], y = median_ISG, fill=.data[[g]])) +
        geom_violin(trim = FALSE, color = "black") +
        stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black", fatten = 0) +
        facet_wrap(~ category, scales = "fixed", nrow=1) +
        theme_classic() +
        labs(title = paste0("ISG expression across ", g, " by major population")) + 
        stat_compare_means(comparisons = comparisons, method = "wilcox.test") + 
        theme(legend.position = "none")

    if(g == "disease_inflam"){
        p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
    }

    p

    ggsave(paste0(out, "ISG_per_", g, "_per_category_persamp.pdf"),
         width = 12, height = 7,device = cairo_pdf)
}


ggplot(mstisamp, aes(x = disease_status, y = median_ISG, fill=disease_status)) +
  geom_violin(trim = FALSE, color = "black") +
  facet_wrap(~ category, scales = "fixed", nrow=1) +
  theme_classic() +
  labs(title = "ISG expression across disease status by major population") + 
  stat_compare_means(comparisons = comparisons, method = "wilcox.test") + 
  theme(legend.position = "none")


ggsave(paste0(out, "ISG_per_disease_status_per_category_persamp.pdf"),
         width = 12, height = 4,device = cairo_pdf)








# What if we group epithelial cells together, ignore the disease_status
msampgrouped = ms %>% 
    select(-predicted_labels) %>% 
    mutate(
        lineage = gsub("Colonocyte|Enterocyte|Secretory|Stem", "Epithelial", category)
    ) %>% 
    group_by(lineage, sanger_sample_id) %>% 
    mutate(
        median_ISG = median(ISG)
    ) %>% 
    select(-c(ISG, cell)) %>% 
    distinct()

ggplot(msampgrouped, aes(x = lineage, y = median_ISG, fill=lineage)) +
  geom_violin(trim = FALSE, color = "black") +
  theme_classic() +
  labs(title = "ISG expression per lineage - per sample")

ggsave(paste0(out, "ISG_per_lineage_persamp.pdf"),
         width = 8, height = 4,device = cairo_pdf)