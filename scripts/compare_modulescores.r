# Bradley August 2025
# Comparing the module scores for different gene sets

# Packages
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)

# Options
scale=F
indir="results/"
ifnmodulename="schogginsISG_fixed"
add_meta_file="input/eQTL_interactions.tsv" 
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
umap.category.palette = deframe(dplyr::select(umap.palette.df, category, category_color))

# Define other variables
out=paste0("plotout/comparion/")
if(scale){
    indir = paste0(indir, "scaled-")
    out = paste0(out, "scaled-")
}

# Make sure outdir exists
if(!file.exists(out)){
    dir.create(out)
}

# Define load and grouping function
load_and_group_module_scores = function(indir, name, add_meta_file=NULL, group_by){
    message(paste0("..Running on ", name))
    # Construct path
    allres=paste0(indir, name, "-module_scores.txt.gz")

    # Load
    message("..Loading")
    ms = fread(allres) %>% 
        mutate( # Couple adjustment
            predicted_labels = gsub("_ct", "", predicted_labels),
            category = gsub("_ct", "", predicted_category),
            disease_status = factor(disease_status, levels = c("Healthy", "CD"))
        ) %>% # Add the proper names
        left_join(select(annot.mapping, c('predicted_labels', 'label_new')))

    # Add meta
    if(!is.null(add_meta_file)){
        message("..Adding meta")
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
    }

    # Group
    message("..Grouping")
    msamp = ms %>% 
        select(-predicted_labels) %>% 
        group_by(across(all_of(group_by))) %>% 
        mutate(
            median_score = median(!!sym(name))
        ) %>% 
        select(-c(!!name, cell)) %>% 
        distinct()

    return(msamp)
}

##############
# Prepping the IFN module
##############
ifnmodulescores = load_and_group_module_scores(indir=indir, name=ifnmodulename, add_meta_file=add_meta_file, group_b=c("category", "sanger_sample_id")) %>% 
    rename(median_ISG = median_score) # adjust the name manually

##################
# Load in the other modules, compute correlations
##################
other_module_names = c("IBD_dsRNA_coloc_SHORT", "magma_dsRNAs_LONG")
bmslist = vector("list", length = length(other_module_names))
comp_res = vector("list", length = length(other_module_names))
names(bmslist) = other_module_names # To store the module scores side by side
names(comp_res) = other_module_names # To store the results
for(other in other_module_names){
    print(paste0("..Working on ", other))
    
    # Load the module data
    othermodulescore = load_and_group_module_scores(indir=indir, name=other, add_meta_file=add_meta_file, group_by=c("category", "sanger_sample_id")) %>% 
        rename(!!paste0("median_", other) := median_score)

    # Merge with the results from the IFN module
    print("..Merging module scores")
    bms = ifnmodulescores %>% 
        select(sanger_sample_id, category, median_ISG) %>% 
        left_join(
            othermodulescore %>% 
                select(sanger_sample_id, category, !!paste0("median_", other)), 
            by=c("category", "sanger_sample_id")) %>% 
        distinct()

    bmslist[[which(other_module_names == other)]] <- bms

    # Compute correlations within categories
    print("..Computing correlations")
    cats = unique(bms$category)
    res = NULL
    for(cat in cats){
        dat = bms %>% 
            filter(category == !!cat)
        corres = cor.test(as.numeric(dat$median_ISG), as.numeric(unlist(dat[,paste0("median_", other)][,1])))
        corres = data.frame(var1="median_ISG", var2=paste0("median_", other), category = cat, cor = corres$estimate, p=corres$p.value)
        if(is.null(res)){
            res = corres
        } else {
            res = rbind(res, corres)
        }
    }

    # Save
    comp_res[[which(other_module_names == other)]] <- res
}
res = do.call(rbind, comp_res)

# Plot one with most significance
top = res %>% 
    slice_min(p)

dat = bmslist[[gsub("median_", "", top$var2)]] %>% 
    filter(category == top$category)

othername = top$var2
ggplot(dat, aes(x = median_ISG, y = .data[[othername]])) +
    geom_point() +
    theme_classic() +
    labs(title = top$category[1]) + 
    annotate(
        "text",
        x = -Inf, y = Inf,
        label = paste0("cor = ", signif(top$cor, 2),
                   "\np = ", signif(top$p, 2)),
        hjust = -0.1, vjust = 1.1,
        size = 5
    )

ggsave(paste0(out, "top_correlation.png"),
         width = 4, height = 4)



# Plot the least
bottom = res %>% 
    slice_max(p)

dat = bmslist[[gsub("median_", "", bottom$var2)]] %>% 
    filter(category == bottom$category)

othername = bottom$var2
ggplot(dat, aes(x = median_ISG, y = .data[[othername]])) +
    geom_point() +
    theme_classic() +
    labs(title = bottom$category[1]) + 
    annotate(
        "text",
        x = -Inf, y = Inf,
        label = paste0("cor = ", signif(bottom$cor, 2),
                   "\np = ", signif(bottom$p, 2)),
        hjust = -0.1, vjust = 1.1,
        size = 5
    )

ggsave(paste0(out, "bottom_correlation.png"),
         width = 4, height = 4)

# Save the results
res = res %>% arrange(p)
write.table(res, paste0(out, "correlations_results.csv"), quote=F, sep = "\t", row.names=F)