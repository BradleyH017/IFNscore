# 
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
inputf="results/379ISGs-module_scores.txt.gz"
out="plotout/"

# Load
ms = fread(inputf)

# Adjust
ms = ms %>% 
    mutate(
        category = gsub("_ct", "", predicted_category),
        disease_status = factor(disease_status, levels = c("Healthy", "CD"))
    ) %>% 
    rename(ISG=`379ISGs`)

# Make sure outdir exists
if(!file.exists(out)){
    dir.create(out)
}

# Plot
ggplot(ms, aes(x = disease_status, y = ISG, fill=disease_status)) +
  geom_violin(trim = FALSE, color = "black") +
  facet_wrap(~ category, scales = "free_y") +
  theme_classic() +
  labs(title = "ISG expression across disease status by category")

ggsave(paste0(out, "ISG_per_disease_status_per_category.pdf"),
         width = 10, height = 4,device = cairo_pdf)

# Group per sample
msamp = ms %>% 
    select(-predicted_labels) %>% 
    group_by(category, sanger_sample_id) %>% 
    mutate(
        median_ISG = median(ISG)
    ) %>% 
    select(-c(ISG, cell)) %>% 
    distinct()

vals = distinct(as.data.frame(msamp$disease_status))[,1]
comparisons = combn(as.character(vals), 2, simplify = FALSE)

ggplot(msamp, aes(x = disease_status, y = median_ISG, fill=disease_status)) +
  geom_violin(trim = FALSE, color = "black") +
  facet_wrap(~ category, scales = "fixed", nrow=1) +
  theme_classic() +
  labs(title = "ISG expression across disease status by category - per sample") + 
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")


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