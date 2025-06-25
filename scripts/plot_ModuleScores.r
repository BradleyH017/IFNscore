# 
library(data.table)
library(ggplot2)
library(dplyr)
inputf="results/ISG-module_scores.txt.gz"
out="plotout/"

# Load
ms = fread(inputf)

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
    group_by(category, sanger_sample_id) %>% 
    mutate(
        median_ISG = median(ISG)
    ) %>% 
    select(-c(ISG, cell)) %>% 
    distinct()

ggplot(msamp, aes(x = disease_status, y = median_ISG, fill=disease_status)) +
  geom_violin(trim = FALSE, color = "black") +
  facet_wrap(~ category, scales = "free_y") +
  theme_classic() +
  labs(title = "ISG expression across disease status by category - per sample")

ggsave(paste0(out, "ISG_per_disease_status_per_category_persamp.pdf"),
         width = 10, height = 4,device = cairo_pdf)