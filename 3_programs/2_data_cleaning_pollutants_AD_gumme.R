# 2_variables_creation_PE
# A. Davias
# 19/10/2021

## Chargement des packages
library(tidyverse)
library(GGally)
library(gtsummary)
library(summarytools)
library(patchwork)
library(ggpubr)
library(grid)
library(questionr)
library(Hmisc)
library(rmarkdown)
library(knitr)
library(labelled)
library(distill)
library(rmdformats)
theme_gtsummary_language("en", decimal.mark = ".", big.mark = " ")
theme_gtsummary_compact(set_theme = TRUE)

## Chargement des données
source("3_programs/1_data_cleaning_covariates_AD_gumme.R", encoding = 'UTF-8')




# PHENOLS ----
## Choix des variables ----
data_phenols <- metadata %>% select(
  
  ident, 
  
  mo_MEPA_total_i_cor_t2, 
  mo_MEPA_total_i_cor_t3, 
  ch_MEPA_total_i_cor_M2, 
  ch_MEPA_total_i_cor_Y1, 
  
  mo_ETPA_total_i_cor_t2,  
  mo_ETPA_total_i_cor_t3, 
  ch_ETPA_total_cat_M2,                                      
  ch_ETPA_total_i_cor_Y1, 
  
  mo_PRPA_total_i_cor_t2, 
  mo_PRPA_total_i_cor_t3, 
  ch_PRPA_total_cat_M2,                                      
  ch_PRPA_total_i_cor_Y1, 
  
  mo_BUPA_total_cat_t2,                                          
  mo_BUPA_total_cat_t3,                                          
  ch_BUPA_total_cat_M2,                                      
  ch_BUPA_total_cat_Y1,                                      
  
  mo_BPA_total_i_cor_t2, 
  mo_BPA_total_i_cor_t3, 
  ch_BPA_total_i_cor_M2, 
  ch_BPA_total_i_cor_Y1, 
  
  mo_BPS_total_cat_t2,                                          
  mo_BPS_total_cat_t3,                                          
  ch_BPS_total_cat_M2,                                      
  ch_BPS_total_cat_Y1,                                      

  mo_OXBE_total_i_cor_t2, 
  mo_OXBE_total_i_cor_t3, 
  ch_OXBE_total_i_cor_M2, 
  ch_OXBE_total_i_cor_Y1, 
  
  mo_TRCS_total_i_cor_t2, 
  mo_TRCS_total_i_cor_t3, 
  ch_TRCS_total_i_cor_M2,                                      
  ch_TRCS_total_i_cor_Y1                                       
)



# Création de vecteurs cat, num et ln 
phenols_vec <- data_phenols %>% select(-ident) %>% colnames()
phenols_vec_cat <- data_phenols %>% select(contains("cat")) %>% colnames()
phenols_vec_num <- data_phenols %>% select(contains("total_i_cor")) %>% colnames()
phenols_vec_ln <- paste(phenols_vec_num, "ln", sep = "_")
phenols_vec_ter <- paste(phenols_vec_num, "ter", sep = "_")

## Création phénols ln ----
metadata[, phenols_vec_num] <- lapply(metadata[, phenols_vec_num], as.numeric)
## les NAs "introduitsé proviennent des phenols codées en character --> numérique qui avaient des NAs
phenols_ln <- metadata %>% 
  select(ident, all_of(phenols_vec_num))%>% 
  select(ident, everything())
phenols_ln[,phenols_vec_num] <- log(phenols_ln[,phenols_vec_num])
colnames(phenols_ln) <- c("ident", paste(colnames(phenols_ln[, 2:23]), "ln", sep="_")) 
metadata <- left_join(metadata, 
                  phenols_ln, 
                  by ="ident")

## Création phénols tertiles ----
phenols_ter <- metadata %>% 
  select(ident, all_of(phenols_vec_num))%>% 
  select(ident, everything()) 
phenols_ter[, phenols_vec_num] <- lapply(phenols_ter[, phenols_vec_num], 
                                             quant.cut, 
                                             nbclass = 3, 
                                             include.lowest = TRUE,
                                             right = FALSE,
                                             dig.lab = 2)
colnames(phenols_ter) <- c("ident", paste(colnames(phenols_ter[, 2:23]), "ter", sep="_"))
metadata <- left_join(metadata, 
                  phenols_ter, 
                  by ="ident")

## Création étiquettes catégories ----
metadata[, phenols_vec_cat] <- lapply(metadata[, phenols_vec_cat], as.character) 
metadata[, phenols_vec_cat]  <- metadata[, phenols_vec_cat] %>% mutate_all(na_if, c("", "NA"))
metadata[, phenols_vec_cat] <- lapply(metadata[, phenols_vec_cat], as.factor) 



metadata <- metadata %>%
  mutate(
  mo_MEPA_total_i_cor_t2_ter = fct_recode(                                      # MEPA
    mo_MEPA_total_i_cor_t2_ter,
    "1st tertile" = "[1.1,6.6)",
    "2nd tertile" = "[6.6,23)",
    "3rd tertile" = "[23,9.8e+03]"
  ),
  mo_MEPA_total_i_cor_t3_ter = fct_recode(
    mo_MEPA_total_i_cor_t3_ter,
    "1st tertile" = "[0.93,6.7)",
    "2nd tertile" = "[6.7,26)",
    "3rd tertile" = "[26,1.8e+04]"
  ),
  ch_MEPA_total_i_cor_M2_ter = fct_recode(
    ch_MEPA_total_i_cor_M2_ter,
    "1st tertile" = "[0.03,0.48)",
    "2nd tertile" = "[0.48,4.8)",
    "3rd tertile" = "[4.8,5.1e+03]"
  ),
  ch_MEPA_total_i_cor_Y1_ter = fct_recode(
    ch_MEPA_total_i_cor_Y1_ter,
    "1st tertile" = "[0.17,16)",
    "2nd tertile" = "[16,2.1e+02)",
    "3rd tertile" = "[2.1e+02,2.7e+04]"
  ),
  mo_ETPA_total_i_cor_t2_ter = fct_recode(                                      # ETPA
    mo_ETPA_total_i_cor_t2_ter,
    "1st tertile" = "[0.018,0.62)",
    "2nd tertile" = "[0.62,1.3)",
    "3rd tertile" = "[1.3,6.5e+02]"
  ),
  mo_ETPA_total_i_cor_t3_ter = fct_recode(
    mo_ETPA_total_i_cor_t3_ter,
    "1st tertile" = "[0.025,0.62)",
    "2nd tertile" = "[0.62,1.4)",
    "3rd tertile" = "[1.4,1.6e+03]"
  ),
  ch_ETPA_total_cat_M2 = fct_recode(
    ch_ETPA_total_cat_M2,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  ch_ETPA_total_i_cor_Y1_ter = fct_recode(
    ch_ETPA_total_i_cor_Y1_ter,
    "1st tertile" = "[0.022,1.1)",
    "2nd tertile" = "[1.1,12)",
    "3rd tertile" = "[12,1.5e+03]"
  ),
  mo_PRPA_total_i_cor_t2_ter = fct_recode(                                      # PRPA
    mo_PRPA_total_i_cor_t2_ter,
    "1st tertile" = "[9.2e-05,0.1)",
    "2nd tertile" = "[0.1,1.6)",
    "3rd tertile" = "[1.6,1.2e+03]"
  ),
  mo_PRPA_total_i_cor_t3_ter = fct_recode(
    mo_PRPA_total_i_cor_t3_ter,
    "1st tertile" = "[4.8e-05,0.094)",
    "2nd tertile" = "[0.094,1.9)",
    "3rd tertile" = "[1.9,1.3e+03]"
  ),
  ch_PRPA_total_cat_M2 = fct_recode(
    ch_PRPA_total_cat_M2,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  ch_PRPA_total_i_cor_Y1_ter = fct_recode(
    ch_PRPA_total_i_cor_Y1_ter,
    "1st tertile" = "[0.00023,0.71)",
    "2nd tertile" = "[0.71,8.8)",
    "3rd tertile" = "[8.8,1.1e+04]"
  ),
  mo_BUPA_total_cat_t2 = fct_recode(                                            # BUPA
    mo_BUPA_total_cat_t2,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  mo_BUPA_total_cat_t3 = fct_recode(
    mo_BUPA_total_cat_t3,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  ch_BUPA_total_cat_M2 = fct_recode(
    ch_BUPA_total_cat_M2,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  ch_BUPA_total_cat_Y1 = fct_recode(
    ch_BUPA_total_cat_Y1,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ), 
  mo_BPA_total_i_cor_t2_ter = fct_recode(                                       # BPA
    mo_BPA_total_i_cor_t2_ter,
    "1st tertile" = "[0.028,1.4)",
    "2nd tertile" = "[1.4,2.5)",
    "3rd tertile" = "[2.5,71]"
  ),
  mo_BPA_total_i_cor_t3_ter = fct_recode(
    mo_BPA_total_i_cor_t3_ter,
    "1st tertile" = "[0.018,1.2)",
    "2nd tertile" = "[1.2,2.4)",
    "3rd tertile" = "[2.4,41]"
  ),
  ch_BPA_total_i_cor_M2_ter = fct_recode(
    ch_BPA_total_i_cor_M2_ter,
    "1st tertile" = "[0.018,0.53)",
    "2nd tertile" = "[0.53,1)",
    "3rd tertile" = "[1,7.5]"
  ),
  ch_BPA_total_i_cor_Y1_ter = fct_recode(
    ch_BPA_total_i_cor_Y1_ter,
    "1st tertile" = "[0.16,2.1)",
    "2nd tertile" = "[2.1,3.8)",
    "3rd tertile" = "[3.8,44]"
  ),
  mo_BPS_total_cat_t2 = fct_recode(                                             # BPS
    mo_BPS_total_cat_t2,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ), 
  mo_BPS_total_cat_t3 = fct_recode(
    mo_BPS_total_cat_t3,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  ch_BPS_total_cat_M2 = fct_recode(
    ch_BPS_total_cat_M2,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  ch_BPS_total_cat_Y1 = fct_recode(
    ch_BPS_total_cat_Y1,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"
  ),
  mo_OXBE_total_i_cor_t2_ter = fct_recode(                                      # OXBE
    mo_OXBE_total_i_cor_t2_ter,
    "1st tertile" = "[0.066,0.54)",
    "2nd tertile" = "[0.54,1.6)",
    "3rd tertile" = "[1.6,1.4e+03]"
  ),
  mo_OXBE_total_i_cor_t3_ter = fct_recode(
    mo_OXBE_total_i_cor_t3_ter,
    "1st tertile" = "[0.014,0.47)",
    "2nd tertile" = "[0.47,1.2)",
    "3rd tertile" = "[1.2,1.5e+03]"
  ),
  ch_OXBE_total_i_cor_M2_ter = fct_recode(
    ch_OXBE_total_i_cor_M2_ter,
    "1st tertile" = "[0.0075,0.22)",
    "2nd tertile" = "[0.22,0.45)",
    "3rd tertile" = "[0.45,1.7e+02]"
  ),
  ch_OXBE_total_i_cor_Y1_ter = fct_recode(
    ch_OXBE_total_i_cor_Y1_ter,
    "1st tertile" = "[0.039,0.28)",
    "2nd tertile" = "[0.28,0.52)",
    "3rd tertile" = "[0.52,1.1e+02]"
  ),
  mo_TRCS_total_i_cor_t2_ter = fct_recode(                                      # TRCS
    mo_TRCS_total_i_cor_t2_ter,
    "1st tertile" = "[0.014,0.57)",
    "2nd tertile" = "[0.57,1.6)",
    "3rd tertile" = "[1.6,2.7e+03]"
  ),
  mo_TRCS_total_i_cor_t3_ter = fct_recode(
    mo_TRCS_total_i_cor_t3_ter,
    "1st tertile" = "[0.011,0.58)",
    "2nd tertile" = "[0.58,1.4)",
    "3rd tertile" = "[1.4,8.4e+02]"
  ),
  ch_TRCS_total_i_cor_M2_ter = fct_recode(
    metadata$ch_TRCS_total_i_cor_M2_ter,
    "1st tertile" = "[0.0047,0.18)",
    "2nd tertile" = "[0.18,0.34)",
    "3rd tertile" = "[0.34,27]"
  ),
  ch_TRCS_total_i_cor_Y1_ter = fct_recode(
    ch_TRCS_total_i_cor_Y1_ter,
    "1st tertile" = "[0.01,0.19)",
    "2nd tertile" = "[0.19,0.35)",
    "3rd tertile" = "[0.35,63]"
  )
)


## Création de phénols à deux catégories ----
metadata <- metadata %>% mutate(
  ch_ETPA_total_cat_M2_2 = fct_recode(ch_ETPA_total_cat_M2, 
                                      "<LOQ" = "<LOD",
                                      "<LOQ" = "LOD-LOQ"),
  ch_PRPA_total_cat_M2_2 = fct_recode(ch_PRPA_total_cat_M2,
                                       ">LOD" = "LOD-LOQ",
                                       ">LOD" = ">LOQ"),
  ch_BUPA_total_cat_M2_2 = fct_recode(ch_BUPA_total_cat_M2,
                                       ">LOD" = "LOD-LOQ",
                                       ">LOD" = ">LOQ"),
  mo_BPS_total_cat_t2_2 = fct_recode(mo_BPS_total_cat_t2,
                                      ">LOD" = "LOD-LOQ",
                                      ">LOD" = ">LOQ"), 
  mo_BPS_total_cat_t3_2 = fct_recode(mo_BPS_total_cat_t3,
                                      ">LOD" = "LOD-LOQ",
                                      ">LOD" = ">LOQ"), 
  ch_BPS_total_cat_M2_2 = fct_recode(ch_BPS_total_cat_M2,
                                      ">LOD" = "LOD-LOQ",
                                      ">LOD" = ">LOQ"))
  

metadata <- metadata %>%
  mutate( 
    mo_BUPA_total_cat_t2 = fct_relevel(
      mo_BUPA_total_cat_t2,
      "<LOD", "LOD-LOQ", ">LOQ"), 
    mo_BUPA_total_cat_t3 = fct_relevel(
      mo_BUPA_total_cat_t3,
      "<LOD", "LOD-LOQ", ">LOQ"), 
    ch_BUPA_total_cat_Y1 = fct_relevel(
      ch_BUPA_total_cat_Y1,
      "<LOD", "LOD-LOQ", ">LOQ"), 
    ch_BPS_total_cat_Y1 = fct_relevel(
      ch_BPS_total_cat_Y1,
      "<LOD", "LOD-LOQ", ">LOQ")) 

# PFAS --------
## Choix des variables ----
metadata <- metadata %>%
  rename(
    mo_PFDA_i_cor = mo_pfda_i_cor,
    mo_PFUnDA_i_cor = mo_pfunda_i_cor,
    mo_PFHpS_i_cor = mo_pfhps_i_cor,
    
    mo_PFDoDa_cat = mo_PFDoDA_cat,
    mo_PFTrDa_cat = mo_PFTrDA_cat,
    
    mo_6_2diPAP_cat = mo__6_2diPAP_cat,
    mo_8_2diPAP_cat = mo__8_2diPAP_cat
  )



data_pfas <- metadata %>% select(ident, 
                                 mo_PFOA_cor, 
                                 mo_PFNA, 
                                 mo_PFHxS_cor,
                                 mo_PFOS_cor,
                                 
                                 mo_PFDA_i_cor, 
                                 mo_PFUnDA_i_cor, 
                                 mo_PFHpS_i_cor,
                                 
                                 mo_PFDoDa_cat, 
                                 mo_PFHxPA_cat, 
                                 mo_PFHpA_cat,
                                 mo_PFTrDa_cat, 
                                 mo_PFBS_cat, 
                                 mo_PFOSA_cat, 
                                 mo_6_2diPAP_cat,
                                 mo_8_2diPAP_cat)


pfas_vec_num <- c("mo_PFOA_cor", 
                  "mo_PFNA", 
                  "mo_PFHxS_cor",
                  "mo_PFOS_cor",
                  
                  "mo_PFDA_i_cor", 
                  "mo_PFUnDA_i_cor", 
                  "mo_PFHpS_i_cor")
pfas_vec_ln <- paste(pfas_vec_num, "ln", sep = "_")
pfas_vec_ter <- paste(pfas_vec_num, "ter", sep = "_")
pfas_vec_cat <- c("mo_PFDoDa_cat", 
                  "mo_PFHxPA_cat", 
                  "mo_PFHpA_cat",
                  "mo_PFTrDa_cat", 
                  "mo_PFBS_cat", 
                  "mo_PFOSA_cat", 
                  "mo_6_2diPAP_cat",
                  "mo_8_2diPAP_cat")






## Création pfas ln ------
pfas_ln <- metadata %>% select(ident, all_of(pfas_vec_num))
pfas_ln[,pfas_vec_num] <- log(pfas_ln[,pfas_vec_num])
colnames(pfas_ln) <- c("ident", paste(colnames(pfas_ln[, 2:8]), "ln", sep="_"))   
metadata <- left_join(metadata, 
                   pfas_ln, 
                   by ="ident")


## Création pfas tertiles ----
pfas_ter <- metadata %>% 
  select(ident, all_of(pfas_vec_num))%>% 
  select(ident, everything()) 
pfas_ter[, pfas_vec_num] <- lapply(pfas_ter[, pfas_vec_num], 
                                         quant.cut, 
                                         nbclass = 3, 
                                         include.lowest = TRUE,
                                         right = FALSE,
                                         dig.lab = 2)
colnames(pfas_ter) <- c("ident", paste(colnames(pfas_ter[, 2:8]), "ter", sep="_"))
metadata <- left_join(metadata, 
                      pfas_ter, 
                      by ="ident")

## Création étiquettes catégories ----
metadata[, pfas_vec_cat] <- lapply(metadata[, pfas_vec_cat], as.character) 
metadata[, pfas_vec_cat] <- lapply(metadata[, pfas_vec_cat], as.factor) 


metadata <- mutate(
  metadata, 
  mo_PFOS_cor_ter = fct_recode(mo_PFOS_cor_ter,
                           "1st tertile" = "[0.76,2.9)",
                           "2nd tertile" = "[2.9,4.5)",
                           "3rd tertile" = "[4.5,79]"),
  mo_PFOA_cor_ter =  fct_recode(mo_PFOA_cor_ter,
                            "1st tertile" = "[0.14,0.86)",
                            "2nd tertile" = "[0.86,1.3)",
                            "3rd tertile" = "[1.3,19]"), 
  mo_PFNA_ter = fct_recode(mo_PFNA_ter,
                           "1st tertile" = "[0.17,0.42)",
                           "2nd tertile" = "[0.42,0.59)",
                           "3rd tertile" = "[0.59,4.3]"), 
  mo_PFDA_i_cor_ter = fct_recode(mo_PFDA_i_cor_ter,
                           "1st tertile" = "[0.039,0.19)",
                           "2nd tertile" = "[0.19,0.26)",
                           "3rd tertile" = "[0.26,2.7]"), 
  mo_PFHxS_cor_ter = fct_recode(mo_PFHxS_cor_ter,
                            "1st tertile" = "[0.13,0.53)",
                            "2nd tertile" = "[0.53,0.77)",
                            "3rd tertile" = "[0.77,1e+02]"), 
  mo_PFHpS_i_cor_ter = fct_recode(mo_PFHpS_i_cor_ter,
                            "1st tertile" = "[0.013,0.071)",
                            "2nd tertile" = "[0.071,0.11)",
                            "3rd tertile" = "[0.11,6.3]"), 
  mo_PFUnDA_i_cor_ter = fct_recode(mo_PFUnDA_i_cor_ter,
                             "1st tertile" = "[0.026,0.12)",
                             "2nd tertile" = "[0.12,0.18)",
                             "3rd tertile" = "[0.18,0.73]"), 
  mo_PFDoDa_cat = fct_recode(mo_PFDoDa_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"), 
  mo_PFHxPA_cat = fct_recode(mo_PFHxPA_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"), 
  mo_PFHpA_cat = fct_recode(mo_PFHpA_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"), 
  mo_PFTrDa_cat = fct_recode(mo_PFTrDa_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"), 
  mo_PFBS_cat = fct_recode(mo_PFBS_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"), 
  mo_PFOSA_cat = fct_recode(mo_PFOSA_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"), 
  mo_6_2diPAP_cat = fct_recode(mo_6_2diPAP_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"), 
  mo_8_2diPAP_cat = fct_recode(mo_8_2diPAP_cat,
                             "<LOD" = "1",
                             "LOD-LOQ" = "2",
                             ">LOQ" = "3"))

## Création de pfas à deux catégories ----
metadata <- mutate(
  metadata, 
  mo_PFHpA_cat_2 = fct_recode(mo_PFHpA_cat,
                                        ">LOD" = "LOD-LOQ",
                                        ">LOD" = ">LOQ"),
  mo_PFBS_cat_2 = fct_recode(mo_PFBS_cat,
                                       ">LOD" = "LOD-LOQ",
                                       ">LOD" = ">LOQ"),
  mo_PFOSA_cat_2 = fct_recode(mo_PFOSA_cat,
                                        ">LOD" = "LOD-LOQ",
                                        ">LOD" = ">LOQ"),
  mo_6_2diPAP_cat_2 = fct_recode(mo_6_2diPAP_cat,
                                           ">LOD" = "LOD-LOQ",
                                           ">LOD" = ">LOQ"),
  mo_8_2diPAP_cat_2 = fct_recode(mo_8_2diPAP_cat,
                                  ">LOD" = "LOD-LOQ",
                                  ">LOD" = ">LOQ")
  )



metadata <- metadata %>%
  mutate( 
mo_PFDoDa_cat = fct_relevel(
  mo_PFDoDa_cat,
  "<LOD", "LOD-LOQ", ">LOQ"),
mo_PFHxPA_cat = fct_relevel(
  mo_PFHxPA_cat,
  "<LOD", "LOD-LOQ", ">LOQ"),
mo_PFTrDa_cat = fct_relevel(
  mo_PFTrDa_cat,
  "<LOD", "LOD-LOQ", ">LOQ"))



# PHTHALATES ----
## Choix des variables ----
data_phthalates <- metadata %>% 
  select(
    
    ident, 
    
    mo_cxMiNP_i_cor_t2,                                                         # Mono-4-methyl-7-carboxyoctyl phthalate 
    mo_cxMiNP_i_cor_t3, 
    ch_cxMiNP_i_cor_M2, 
    ch_cxMiNP_i_cor_Y1,    
    
    mo_MBzP_i_cor_t2,                                                           # Monobenzyl phthalate 
    mo_MBzP_i_cor_t3, 
    ch_MBzP_i_cor_M2, 
    ch_MBzP_i_cor_Y1,
    
    mo_MECPP_i_cor_t2,                                                          # Mono(2-ethyl-5-carboxypentyl) phthalate  
    mo_MECPP_i_cor_t3, 
    ch_MECPP_i_cor_M2, 
    ch_MECPP_i_cor_Y1,
    
    mo_MEHHP_i_cor_t2,                                                          # Mono(2-ethyl-5-hydroxyhexyl) phthalate  
    mo_MEHHP_i_cor_t3, 
    ch_MEHHP_i_cor_M2, 
    ch_MEHHP_i_cor_Y1, 
    
    mo_MEHP_i_cor_t2,                                                           # Mono(2-ethylhexyl) phthalate 
    mo_MEHP_i_cor_t3, 
    ch_MEHP_i_cor_M2, 
    ch_MEHP_i_cor_Y1, 
    
    mo_MEOHP_i_cor_t2,                                                          # Mono(2-ethyl-5oxohexyl)phthalate  
    mo_MEOHP_i_cor_t3, 
    ch_MEOHP_i_cor_M2, 
    ch_MEOHP_i_cor_Y1, 
    
    mo_MEP_i_cor_t2,                                                            # Monoethyl phthalate 
    mo_MEP_i_cor_t3,
    ch_MEP_i_cor_M2, 
    ch_MEP_i_cor_Y1,
    
    mo_MiBP_i_cor_t2,                                                           # Mono-isobutyl phthalate 
    mo_MiBP_i_cor_t3, 
    ch_MiBP_i_cor_M2, 
    ch_MiBP_i_cor_Y1, 
    
    mo_MMCHP_i_cor_t2,                                                          # Mono-2-carboxymethyl hexyl phthalate 
    mo_MMCHP_i_cor_t3, 
    ch_MMCHP_i_cor_M2, 
    ch_MMCHP_i_cor_Y1,
    
    mo_MnBP_i_cor_t2,                                                           # Mono-n-butyl phthalate 
    mo_MnBP_i_cor_t3, 
    ch_MnBP_i_cor_M2, 
    ch_MnBP_i_cor_Y1,
    
    mo_ohMiNP_i_cor_t2,                                                         # Mono-4-methyl-7-hydroxyoctyl phthalate 
    mo_ohMiNP_i_cor_t3, 
    ch_ohMiNP_i_cor_M2,
    ch_ohMiNP_i_cor_Y1,
    
    mo_ohMPHP_i_cor_t2,                                                         # 6-Hydroxy Monopropylheptylphthalate
    mo_ohMPHP_i_cor_t3, 
    ch_ohMPHP_cat_M2, 
    ch_ohMPHP_i_cor_Y1, 
    
    mo_ohMINCH_i_cor_t2, 
    mo_ohMINCH_i_cor_t3, 
    ch_ohMINCH_cat_M2, 
    ch_ohMINCH_i_cor_Y1,

    mo_oxoMINCH_i_cor_t2, 
    mo_oxoMINCH_i_cor_t3, 
    ch_oxoMINCH_cat_M2, 
    ch_oxoMINCH_i_cor_Y1,
    
    mo_oxoMiNP_i_cor_t2,                                                        # Mono-4-methyl-7-oxooctyl phthalate 
    mo_oxoMiNP_i_cor_t3, 
    ch_oxoMiNP_i_cor_M2, 
    ch_oxoMiNP_i_cor_Y1
)


phthalates_vec_num <- data_phthalates %>% select(contains("_i_cor_")) %>% colnames()
phthalates_vec_cat <- c("ch_ohMINCH_cat_M2", "ch_ohMPHP_cat_M2", "ch_oxoMINCH_cat_M2")
phthalates_vec_ln <- paste(phthalates_vec_num, "ln", sep = "_")
phthalates_vec_ter <- paste(phthalates_vec_num, "ter", sep = "_")


## Création phthalates ln ----
metadata[, phthalates_vec_num] <- lapply(metadata[, phthalates_vec_num], as.numeric) 
## les NAs "introduitsé proviennent des phenols codées en character --> numérique qui avaient des NAs
phthalates_ln <- metadata %>% 
  select(ident, all_of(phthalates_vec_num)) %>% 
  select(ident, everything())

phthalates_ln[,phthalates_vec_num] <- log(phthalates_ln[,phthalates_vec_num])
colnames(phthalates_ln) <- c("ident", paste(colnames(phthalates_ln[, 2:58]), "ln", sep="_")) 
metadata <- left_join(metadata,
                  phthalates_ln, 
                  by ="ident")


## Création phthalates tertiles ----
phthalates_ter <- metadata %>% 
  select(ident, all_of(phthalates_vec_num))%>% 
  select(ident, everything()) 
phthalates_ter[, phthalates_vec_num] <- lapply(phthalates_ter[, phthalates_vec_num], 
                                         quant.cut, 
                                         nbclass = 3, 
                                         include.lowest = TRUE,
                                         right = FALSE,
                                         dig.lab = 2)
colnames(phthalates_ter) <- c("ident", paste(colnames(phthalates_ter[, 2:58]), "ter", sep="_"))
metadata <- left_join(metadata, 
                      phthalates_ter, 
                      by ="ident")

## Créattion étiquettes catégories ----
metadata[, phthalates_vec_cat] <- lapply(metadata[, phthalates_vec_cat], as.character) 
metadata[, phthalates_vec_cat]  <- metadata[, phthalates_vec_cat] %>% mutate_all(na_if, c("", "NA"))
metadata[, phthalates_vec_cat] <- lapply(metadata[, phthalates_vec_cat], as.factor) 

metadata <- metadata %>% mutate(
  mo_cxMiNP_i_cor_t2_ter = fct_recode(                                          # cxMiNP
    mo_cxMiNP_i_cor_t2_ter,
    "1st tertile" = "[1.7,3.9)",
    "2nd tertile" = "[3.9,6)",
    "3rd tertile" = "[6,2.3e+02]"
  ),
  mo_cxMiNP_i_cor_t3_ter = fct_recode(
    mo_cxMiNP_i_cor_t3_ter,
    "1st tertile" = "[1.7,3.7)",
    "2nd tertile" = "[3.7,5.4)",
    "3rd tertile" = "[5.4,1.3e+02]"
  ),
  ch_cxMiNP_i_cor_M2_ter = fct_recode(
    ch_cxMiNP_i_cor_M2_ter,
    "1st tertile" = "[0.39,1.6)",
    "2nd tertile" = "[1.6,1.8)",
    "3rd tertile" = "[1.8,6.9]"
  ),
  ch_cxMiNP_i_cor_Y1_ter = fct_recode(
    ch_cxMiNP_i_cor_Y1_ter,
    "1st tertile" = "[0.35,3.2)",
    "2nd tertile" = "[3.2,4.6)",
    "3rd tertile" = "[4.6,1.1e+02]"
  ),
  
  
  mo_MBzP_i_cor_t2_ter = fct_recode(                                            # MBzP
    mo_MBzP_i_cor_t2_ter,
    "1st tertile" = "[0.64,3.5)",
    "2nd tertile" = "[3.5,5.9)",
    "3rd tertile" = "[5.9,1.6e+02]"
  ),
  mo_MBzP_i_cor_t3_ter = fct_recode(
    mo_MBzP_i_cor_t3_ter,
    "1st tertile" = "[0.53,3.2)",
    "2nd tertile" = "[3.2,5.8)",
    "3rd tertile" = "[5.8,1.4e+02]"
  ),
  ch_MBzP_i_cor_M2_ter = fct_recode(
    ch_MBzP_i_cor_M2_ter,
    "1st tertile" = "[0.066,0.57)",
    "2nd tertile" = "[0.57,1.1)",
    "3rd tertile" = "[1.1,5.3]"
  ),
  ch_MBzP_i_cor_Y1_ter = fct_recode(
    ch_MBzP_i_cor_Y1_ter,
    "1st tertile" = "[0.061,2.1)",
    "2nd tertile" = "[2.1,5.8)",
    "3rd tertile" = "[5.8,2.1e+02]"
  ),
  
  
  mo_MECPP_i_cor_t2_ter = fct_recode(                                           # MECPP
    mo_MECPP_i_cor_t2_ter,
    "1st tertile" = "[1.2,8.4)",
    "2nd tertile" = "[8.4,12)",
    "3rd tertile" = "[12,3.1e+02]"
  ),
  mo_MECPP_i_cor_t3_ter = fct_recode(
    mo_MECPP_i_cor_t3_ter,
    "1st tertile" = "[2.8,8.2)",
    "2nd tertile" = "[8.2,13)",
    "3rd tertile" = "[13,1e+02]"
  ),
  ch_MECPP_i_cor_M2_ter = fct_recode(
    ch_MECPP_i_cor_M2_ter,
    "1st tertile" = "[3.5,5.7)",
    "2nd tertile" = "[5.7,7.8)",
    "3rd tertile" = "[7.8,62]"
  ),
  ch_MECPP_i_cor_Y1_ter = fct_recode(
    ch_MECPP_i_cor_Y1_ter,
    "1st tertile" = "[0.56,7.9)",
    "2nd tertile" = "[7.9,16)",
    "3rd tertile" = "[16,2.1e+02]"
  ),
  
  
  mo_MEHHP_i_cor_t2_ter = fct_recode(                                           # MEHHP
    mo_MEHHP_i_cor_t2_ter,
    "1st tertile" = "[1.3,5.7)",
    "2nd tertile" = "[5.7,9.2)",
    "3rd tertile" = "[9.2,2.9e+02]"
  ),
  mo_MEHHP_i_cor_t3_ter = fct_recode(
    mo_MEHHP_i_cor_t3_ter,
    "1st tertile" = "[1.5,5.5)",
    "2nd tertile" = "[5.5,9.2)",
    "3rd tertile" = "[9.2,1e+02]"
  ),
  ch_MEHHP_i_cor_M2_ter = fct_recode(
    ch_MEHHP_i_cor_M2_ter,
    "1st tertile" = "[0.31,0.94)",
    "2nd tertile" = "[0.94,1.4)",
    "3rd tertile" = "[1.4,9]"
  ),
  ch_MEHHP_i_cor_Y1_ter = fct_recode(
    ch_MEHHP_i_cor_Y1_ter,
    "1st tertile" = "[0.38,3.8)",
    "2nd tertile" = "[3.8,7.7)",
    "3rd tertile" = "[7.7,1.5e+02]"
  ),
  
  
  mo_MEHP_i_cor_t2_ter = fct_recode(                                            # MEHP
    mo_MEHP_i_cor_t2_ter,
    "1st tertile" = "[0.25,1.7)",
    "2nd tertile" = "[1.7,3.2)",
    "3rd tertile" = "[3.2,1.3e+02]"
  ),
  mo_MEHP_i_cor_t3_ter = fct_recode(
    mo_MEHP_i_cor_t3_ter,
    "1st tertile" = "[0.2,1.3)",
    "2nd tertile" = "[1.3,2.6)",
    "3rd tertile" = "[2.6,23]"
  ),
  ch_MEHP_i_cor_M2_ter = fct_recode(
    ch_MEHP_i_cor_M2_ter,
    "1st tertile" = "[0.19,1.5)",
    "2nd tertile" = "[1.5,2.5)",
    "3rd tertile" = "[2.5,8.3]"
  ),
  ch_MEHP_i_cor_Y1_ter = fct_recode(
    ch_MEHP_i_cor_Y1_ter,
    "1st tertile" = "[0.12,1.5)",
    "2nd tertile" = "[1.5,3)",
    "3rd tertile" = "[3,27]"
  ),
  
  
  mo_MEOHP_i_cor_t2_ter = fct_recode(                                           # MEOHP
    mo_MEOHP_i_cor_t2_ter,
    "1st tertile" = "[0.92,4)",
    "2nd tertile" = "[4,6.7)",
    "3rd tertile" = "[6.7,2e+02]"
  ),
  mo_MEOHP_i_cor_t3_ter = fct_recode(
    mo_MEOHP_i_cor_t3_ter,
    "1st tertile" = "[0.93,4.1)",
    "2nd tertile" = "[4.1,6.6)",
    "3rd tertile" = "[6.6,72]"
  ),
  ch_MEOHP_i_cor_M2_ter = fct_recode(
    ch_MEOHP_i_cor_M2_ter,
    "1st tertile" = "[0.22,0.63)",
    "2nd tertile" = "[0.63,0.94)",
    "3rd tertile" = "[0.94,5.8]"
  ),
  ch_MEOHP_i_cor_Y1_ter = fct_recode(
    ch_MEOHP_i_cor_Y1_ter,
    "1st tertile" = "[0.33,2.6)",
    "2nd tertile" = "[2.6,5.8)",
    "3rd tertile" = "[5.8,96]"
  ),

  
  mo_MEP_i_cor_t2_ter = fct_recode(                                             # MEP
    mo_MEP_i_cor_t2_ter,
    "1st tertile" = "[2.3,15)",
    "2nd tertile" = "[15,37)",
    "3rd tertile" = "[37,1.3e+03]"
  ),
  mo_MEP_i_cor_t3_ter = fct_recode(
    mo_MEP_i_cor_t3_ter,
    "1st tertile" = "[1.4,13)",
    "2nd tertile" = "[13,32)",
    "3rd tertile" = "[32,1.4e+03]"
  ),
  ch_MEP_i_cor_M2_ter = fct_recode(
    ch_MEP_i_cor_M2_ter,
    "1st tertile" = "[0.16,3.3)",
    "2nd tertile" = "[3.3,5.9)",
    "3rd tertile" = "[5.9,78]"
  ),
  ch_MEP_i_cor_Y1_ter = fct_recode(
    ch_MEP_i_cor_Y1_ter,
    "1st tertile" = "[0.2,8.6)",
    "2nd tertile" = "[8.6,17)",
    "3rd tertile" = "[17,1.6e+02]"
  ),
  
  
  mo_MiBP_i_cor_t2_ter = fct_recode(                                            # MiBP
    mo_MiBP_i_cor_t2_ter,
    "1st tertile" = "[2,12)",
    "2nd tertile" = "[12,20)",
    "3rd tertile" = "[20,1.8e+02]"
  ),
  mo_MiBP_i_cor_t3_ter = fct_recode(
    mo_MiBP_i_cor_t3_ter,
    "1st tertile" = "[1.8,12)",
    "2nd tertile" = "[12,19)",
    "3rd tertile" = "[19,1.5e+02]"
  ),
  ch_MiBP_i_cor_M2_ter = fct_recode(
    ch_MiBP_i_cor_M2_ter,
    "1st tertile" = "[1.2,4.2)",
    "2nd tertile" = "[4.2,7.7)",
    "3rd tertile" = "[7.7,28]"
  ),
  ch_MiBP_i_cor_Y1_ter = fct_recode(
    ch_MiBP_i_cor_Y1_ter,
    "1st tertile" = "[0.36,10)",
    "2nd tertile" = "[10,18)",
    "3rd tertile" = "[18,3.6e+02]"
  ),
  
  
  mo_MMCHP_i_cor_t2_ter = fct_recode(                                           # MMCHP
    mo_MMCHP_i_cor_t2_ter,
    "1st tertile" = "[0.49,6.4)",
    "2nd tertile" = "[6.4,9.3)",
    "3rd tertile" = "[9.3,1.8e+02]"
  ),
  mo_MMCHP_i_cor_t3_ter = fct_recode(
    mo_MMCHP_i_cor_t3_ter,
    "1st tertile" = "[0.53,6.3)",
    "2nd tertile" = "[6.3,9.3)",
    "3rd tertile" = "[9.3,66]"
  ),
  ch_MMCHP_i_cor_M2_ter = fct_recode(
    ch_MMCHP_i_cor_M2_ter,
    "1st tertile" = "[0.63,3)",
    "2nd tertile" = "[3,3.6)",
    "3rd tertile" = "[3.6,9.4]"
  ),
  ch_MMCHP_i_cor_Y1_ter = fct_recode(
    ch_MMCHP_i_cor_Y1_ter,
    "1st tertile" = "[3.3,7.5)",
    "2nd tertile" = "[7.5,11)",
    "3rd tertile" = "[11,60]"
  ),
  
  
  mo_MnBP_i_cor_t2_ter = fct_recode(                                            # MnBP
    mo_MnBP_i_cor_t2_ter, 
    "1st tertile" = "[1.8,8.8)",
    "2nd tertile" = "[8.8,14)",
    "3rd tertile" = "[14,1.7e+02]"
  ),
  mo_MnBP_i_cor_t3_ter = fct_recode(
    mo_MnBP_i_cor_t3_ter,
    "1st tertile" = "[1.2,8.7)",
    "2nd tertile" = "[8.7,14)",
    "3rd tertile" = "[14,2e+02]"
  ),
  ch_MnBP_i_cor_M2_ter = fct_recode(
    ch_MnBP_i_cor_M2_ter,
    "1st tertile" = "[1.2,3.8)",
    "2nd tertile" = "[3.8,5.9)",
    "3rd tertile" = "[5.9,31]"
  ),
  ch_MnBP_i_cor_Y1_ter = fct_recode(
    ch_MnBP_i_cor_Y1_ter,
    "1st tertile" = "[0.49,8)",
    "2nd tertile" = "[8,16)",
    "3rd tertile" = "[16,1.5e+02]"
  ),
  
  
  mo_ohMINCH_i_cor_t2_ter = fct_recode(
    mo_ohMINCH_i_cor_t2_ter,
    "1st tertile" = "[0.36,1.3)",
    "2nd tertile" = "[1.3,2.4)",
    "3rd tertile" = "[2.4,1.1e+02]"),
  mo_ohMINCH_i_cor_t3_ter = fct_recode(
    mo_ohMINCH_i_cor_t3_ter,
    "1st tertile" = "[0.21,1.2)",
    "2nd tertile" = "[1.2,2.2)",
    "3rd tertile" = "[2.2,1.3e+02]"),
  ch_ohMINCH_cat_M2 = fct_recode(
    ch_ohMINCH_cat_M2,
    "<LOD" = "1",
    ">LOQ" = "3"),
  ch_ohMINCH_i_cor_Y1_ter = fct_recode(
    ch_ohMINCH_i_cor_Y1_ter,
    "1st tertile" = "[0.22,1.5)",
    "2nd tertile" = "[1.5,3)",
    "3rd tertile" = "[3,1.5e+02]"),
  
  
  mo_ohMiNP_i_cor_t2_ter = fct_recode(                                          # ohMiNP
    mo_ohMiNP_i_cor_t2_ter,
    "1st tertile" = "[0.46,3.5)",
    "2nd tertile" = "[3.5,8)",
    "3rd tertile" = "[8,1.4e+02]"
  ),
  mo_ohMiNP_i_cor_t3_ter = fct_recode(
    mo_ohMiNP_i_cor_t3_ter,
    "1st tertile" = "[0.79,3.2)",
    "2nd tertile" = "[3.2,6.9)",
    "3rd tertile" = "[6.9,2.6e+02]"
  ),
  ch_ohMiNP_i_cor_M2_ter = fct_recode(
    ch_ohMiNP_i_cor_M2_ter,
    "1st tertile" = "[0.26,0.62)",
    "2nd tertile" = "[0.62,1.1)",
    "3rd tertile" = "[1.1,1.5e+02]"
  ),
  ch_ohMiNP_i_cor_Y1_ter = fct_recode(
    ch_ohMiNP_i_cor_Y1_ter,
    "1st tertile" = "[0.042,1.8)",
    "2nd tertile" = "[1.8,4)",
    "3rd tertile" = "[4,66]"
  ),
  
  
  
  mo_ohMPHP_i_cor_t2_ter = fct_recode(                                          # ohMPHP
    mo_ohMPHP_i_cor_t2_ter,
    "1st tertile" = "[0.15,0.76)",
    "2nd tertile" = "[0.76,1)",
    "3rd tertile" = "[1,30]"
  ),
  mo_ohMPHP_i_cor_t3_ter = fct_recode(
    mo_ohMPHP_i_cor_t3_ter,
    "1st tertile" = "[0.17,0.7)",
    "2nd tertile" = "[0.7,0.95)",
    "3rd tertile" = "[0.95,36]"
  ),
  ch_ohMPHP_cat_M2 = fct_recode(
    ch_ohMPHP_cat_M2,
    "<LOD" = "1",
    ">LOQ" = "3"
  ),
  ch_ohMPHP_i_cor_Y1_ter = fct_recode(
    ch_ohMPHP_i_cor_Y1_ter,
    "1st tertile" = "[0.064,0.6)",
    "2nd tertile" = "[0.6,0.91)",
    "3rd tertile" = "[0.91,1.1e+02]"
  ),
  
  
  
  mo_oxoMiNP_i_cor_t2_ter = fct_recode(                                         # oxoMiNP
    mo_oxoMiNP_i_cor_t2_ter,
    "1st tertile" = "[0.14,1.6)",
    "2nd tertile" = "[1.6,2.8)",
    "3rd tertile" = "[2.8,1.2e+02]"
  ),
  mo_oxoMiNP_i_cor_t3_ter = fct_recode(
    mo_oxoMiNP_i_cor_t3_ter,
    "1st tertile" = "[0.42,1.6)",
    "2nd tertile" = "[1.6,2.8)",
    "3rd tertile" = "[2.8,64]"
  ),
  ch_oxoMiNP_i_cor_M2_ter = fct_recode(
    ch_oxoMiNP_i_cor_M2_ter,
    "1st tertile" = "[0.066,0.2)",
    "2nd tertile" = "[0.2,0.23)",
    "3rd tertile" = "[0.23,1.4]"
  ),
  ch_oxoMiNP_i_cor_Y1_ter = fct_recode(
    ch_oxoMiNP_i_cor_Y1_ter,
    "1st tertile" = "[0.064,0.94)",
    "2nd tertile" = "[0.94,1.8)",
    "3rd tertile" = "[1.8,63]"
  ), 
  mo_oxoMINCH_i_cor_t2_ter = fct_recode(
    mo_oxoMINCH_i_cor_t2_ter,
    "1st tertile" = "[0.08,1.2)",
    "2nd tertile" = "[1.2,2)",
    "3rd tertile" = "[2,69]"),
  mo_oxoMINCH_i_cor_t3_ter = fct_recode(
    mo_oxoMINCH_i_cor_t3_ter,
    "1st tertile" = "[0.13,1.1)",
    "2nd tertile" = "[1.1,2.1)",
    "3rd tertile" = "[2.1,94]"),
  ch_oxoMINCH_cat_M2 = fct_recode(
    ch_oxoMINCH_cat_M2,
    "<LOD" = "1",
    "LOD-LOQ" = "2",
    ">LOQ" = "3"),
  ch_oxoMINCH_i_cor_Y1_ter = fct_recode(
    ch_oxoMINCH_i_cor_Y1_ter,
    "1st tertile" = "[0.066,0.95)",
    "2nd tertile" = "[0.95,2.1)",
    "3rd tertile" = "[2.1,91]")
)



## Création de phthalates à deux catégories -----
metadata <- mutate(
  metadata, 
  ch_ohMINCH_cat_M2_2 = fct_recode(ch_ohMINCH_cat_M2,
                                    ">LOD" = ">LOQ"), 
  ch_ohMPHP_cat_M2_2 = fct_recode(ch_ohMPHP_cat_M2,
                                   ">LOD" = ">LOQ"), 
  ch_oxoMINCH_cat_M2_2 = fct_recode(ch_oxoMINCH_cat_M2,
                                     ">LOD" = "LOD-LOQ",
                                     ">LOD" = ">LOQ"))



write_labelled_csv(metadata, "2_final_data/metadata_gumme_clean_AD_220721.csv", 
                   row.names = FALSE)

