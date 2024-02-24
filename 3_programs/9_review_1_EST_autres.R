# FigS3
# A.Davias
# 23/02/2024


library(tidyverse)
library(broom)
library(forestplot)
library(bkmr)
library(fields)
library(vegan)
library(phyloseq)
library(broom)
library(forestplot)
library(bkmr)
library(fields)
library(psych)
library(expss)
library(SRS)
library(see)
library(writexl)
library(forcats)
library(egg)
library(corrplot)
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')  


bdd <- bdd_taxa %>% 
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(all_of(phenols_vec_ln), 
              all_of(pfas_vec_ln)) %>%
  na.omit()

colnames(bdd) <- colnames(bdd) %>%
  str_replace_all(c(
    "mo_" = "",
    "ch_" = "",
    "_total_i_cor_" = " ",
    "MEPA" = "Methylparaben", 
    "ETPA" = "Ethylparaben", 
    "BUPA" = "Butylparaben", 
    "PRPA" = "Propylparaben", 
    "OXBE" = "Benzophenone 3", 
    "TRCS" = "Triclosan", 
    "BPA" = "Bisphenol A", 
    "BPS" = "Bisphenol S", 
    "_i_cor" = "",
    "_cor" = "", 
    "PFHpS" = "PFHpS t2", 
    "PFHxS" = "PFHxS t2", 
    "PFOS" = "PFOs t2", 
    "PFNA" = "PFNA t2", 
    "PFOA" = "PFOA t2", 
    "PFDA" = "PFDA t2", 
    "PFUnDA" = "PFUnDA t2", 
    "_ln" = ""
  ))

bdd <- round(cor(bdd, 
          use = "pairwise.complete.obs", 
          method = "pearson"), 3)


tiff(filename = "4_output/heatmap_cor_phenols.tiff", units = "mm", width = 250, height = 250, res = 300)
corrplot(bdd, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         tl.cex = 0.8,
         number.digits = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()
