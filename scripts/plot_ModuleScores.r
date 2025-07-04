# 
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
indir="results/"
out="plotout/"
scale=F
if(scale){
    indir = paste0(indir, "scaled-")
    out = paste0(out, "scaled-")
}
allres=paste0(indir, "schogginsISG-module_scores.txt.gz")
tires=paste0(indir, "schogginsISG-TI-module_scores.txt.gz")
modulename="schogginsISG"
umap.palette.df = read.csv("/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/IBDVerse-sc-eQTL-code/data/palette.csv")
umap.category.palette = deframe(dplyr::select(umap.palette.df, category, category_color))

#######################
# What is the overall expression of this signature like?
#######################
# Load
ms = fread(allres)

# Adjust
ms = ms %>% 
    mutate(
        category = gsub("_ct", "", predicted_category),
        disease_status = factor(disease_status, levels = c("Healthy", "CD"))
    )

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


#######################
# How does this vary by disease status?
#######################
# Use data just from the TI - so a seperate analysis
mstisamp = fread(tires) %>% 
    mutate(
        category = gsub("_ct", "", predicted_category),
        disease_status = factor(disease_status, levels = c("Healthy", "CD"))
    ) %>% 
    select(-predicted_labels) %>% 
    group_by(category, sanger_sample_id) %>% 
    mutate(
        median_ISG = median(!!sym(modulename))
    ) %>% 
    select(-c(!!modulename, cell)) %>% 
    distinct()


vals = distinct(as.data.frame(mstisamp$disease_status))[,1]
comparisons = combn(as.character(vals), 2, simplify = FALSE)

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