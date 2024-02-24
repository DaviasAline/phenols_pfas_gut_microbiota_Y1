# 2_data_selection
# A. Davias
# 19/10/2021


# Chargement des packages ----
library(tidyverse)
library(haven)
library(reshape2)
library(openxlsx)
library(GGally)
library(gtsummary)
library(summarytools)
library(patchwork)
library(AICcmodavg)
library(ggpubr)
library(grid)
library(questionr)
library(Hmisc)
library(rmarkdown)
library(knitr)
library(labelled)
library(distill)
library(rmdformats)
library(parameters)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(lazyeval)
library(car)
library(expss)
theme_gtsummary_language("en", decimal.mark = ".", big.mark = " ")
theme_gtsummary_compact(set_theme = TRUE)


# Chargement des données ----
metadata <- read_labelled_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/pollutants_gut_microbiota_Y1/2_final_data/metadata_gumme_clean_AD_220721.csv")
alphadiv_Y1 <- read_labelled_csv("0_source_data/alpha_diversity_ASVbased_Y1_labelled_AD_20220504_32.csv")
asv_taxa <- read_labelled_csv("0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") %>%
  select(ident, 
         ch_feces_ID_Y1, 
         starts_with("ch_feces_rel")) %>%
  select(!starts_with("ch_feces_rel_ASV"))
bdd_sg <- read_sas(
  "0_source_data/old/metadata_220901/base_aline_220901.sas7bdat",                  # base de données SEPAGES
  catalog_file = "0_source_data/formats.sas7bcat") %>%
  select(ident, 
         mo_pool_sg_T1, mo_pool_sg_T3, ch_pool_sg_M2, ch_pool_sg_Y1) 

bdd_outlier_pfas <- read_sas(
  "0_source_data/old/metadata_220909/base_aline_220909.sas7bdat",                  # base de données SEPAGES
  catalog_file = "0_source_data/formats.sas7bcat") %>%
  select(ident, 
         mo_PFHxS_cor_new = mo_PFHxS_cor,
         mo_PFOS_cor_new = mo_PFOS_cor,
         mo_PFHpS_i_cor_new = mo_pfhps_i_cor)

var_lab(metadata$ident) <- NULL
var_lab(metadata$ch_feces_ID_Y1) <- NULL
var_lab(alphadiv_Y1$ident) <- NULL
var_lab(alphadiv_Y1$ch_feces_ID_Y1) <- NULL
var_lab(asv_taxa$ident) <- NULL
var_lab(asv_taxa$ch_feces_ID_Y1) <- NULL


# Autres données importables (seulement si besoin)
#asv_rel <- read_labelled_csv("~/R/GUMME/pollutants_gut_microbiota_Y1_20220601/0_source_data/microbiota_Y1/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") %>%
#  select(ident, 
#         ch_feces_ID_Y1, 
#         starts_with("ch_feces_rel_ASV"))
#
#asv_raw <- read_labelled_csv("~/R/GUMME/pollutants_gut_microbiota_Y1_20220601/0_source_data/microbiota_Y1/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") %>%
#  select(ident, 
#         ch_feces_ID_Y1, 
#         starts_with("ch_feces_raw_ASV"))
#taxa <- read_csv("~/R/GUMME/pollutants_gut_microbiota_Y1_20220601/0_source_data/microbiota_Y1/taxa_table_ASVbased_Y1_AD_20220504_8.csv")


# Nettoyage gravité spécifique xx_pool_sg_xx ----
# A actualiser dans les codes de nettoyage de données au moment du changement de projet
sg_vec <- bdd_sg %>% select(mo_pool_sg_T1, mo_pool_sg_T3, ch_pool_sg_M2, ch_pool_sg_Y1) %>% colnames()
sg_vec_2 <- sg_vec %>% paste("2", sep = "_")
sg_vec_ter <- sg_vec_2 %>% paste("ter", sep = "_")
bdd_sg[, sg_vec] <- lapply(bdd_sg[, sg_vec], as.numeric)

bdd_sg <- bdd_sg %>% 
  mutate(ch_pool_sg_M2_2 = ifelse(ch_pool_sg_M2 >900, NA, ch_pool_sg_M2)) %>%   # Exclusion de ident outlier 
  
  mutate(mo_pool_sg_T1_2 = mo_pool_sg_T1*1000,                                  # Changement de l'unité de la gravité spécifique 
         mo_pool_sg_T3_2 = mo_pool_sg_T3*1000, 
         ch_pool_sg_M2_2 = ch_pool_sg_M2_2*1000, 
         ch_pool_sg_Y1_2 = ch_pool_sg_Y1*1000)



bdd_sg_ter <- bdd_sg %>%                          # Création variables sg sous formes de tertiles 
  select(ident, all_of(sg_vec_2))%>% 
  select(ident, everything()) 
bdd_sg_ter[, sg_vec_2] <- lapply(
  bdd_sg_ter[, sg_vec_2],
  quant.cut,
  nbclass = 3,
  include.lowest = TRUE,
  right = FALSE,
  dig.lab = 2)

colnames(bdd_sg_ter) <- c("ident", paste(colnames(bdd_sg_ter[, 2:5]), "ter", sep="_"))
bdd_sg_ter <- bdd_sg_ter %>%
  mutate(
    mo_pool_sg_T1_2_ter = fct_recode(mo_pool_sg_T1_2_ter,
                                     "1st tertile" = "[1006,1016)",
                                     "2nd tertile" = "[1016,1020)",
                                     "3rd tertile" = "[1020,1035]"), 
    mo_pool_sg_T3_2_ter = fct_recode(mo_pool_sg_T3_2_ter,
                                     "1st tertile" = "[1e+03,1.01e+03)",
                                     "2nd tertile" = "[1.01e+03,1.02e+03)",
                                     "3rd tertile" = "[1.02e+03,1.04e+03]"), 
    ch_pool_sg_M2_2_ter = fct_recode(ch_pool_sg_M2_2_ter,
                                     "1st tertile" = "[1000,1004)",
                                     "2nd tertile" = "[1004,1005)",
                                     "3rd tertile" = "[1005,1011]"), 
    ch_pool_sg_Y1_2_ter = fct_recode(ch_pool_sg_Y1_2_ter,
                                     "1st tertile" = "[1003,1010)",
                                     "2nd tertile" = "[1010,1013)",
                                     "3rd tertile" = "[1013,1029]"))
bdd_sg <- left_join(bdd_sg, 
                    bdd_sg_ter, 
                    by = "ident")
metadata <- left_join(metadata, 
                      bdd_sg, 
                      by = "ident")



# Nettoyage suspicion contamination BPA M2 et Y1, MEPA Y1 ----
# Pour ces 3 dosages, on utilise les concentrations conjuguées dans nos analyses principales à la place des concentrations totales
metadata[, c("ch_BPA_free_i_cor_M2", "ch_BPA_free_i_cor_Y1", "ch_MEPA_free_i_cor_Y1", 
             "ch_BPA_total_i_cor_M2", "ch_BPA_total_i_cor_Y1", "ch_MEPA_total_i_cor_Y1")] <- 
  lapply(metadata[, c("ch_BPA_free_i_cor_M2", "ch_BPA_free_i_cor_Y1", "ch_MEPA_free_i_cor_Y1", 
                      "ch_BPA_total_i_cor_M2", "ch_BPA_total_i_cor_Y1", "ch_MEPA_total_i_cor_Y1")], as.numeric)

metadata <- metadata %>%                         # create conjugated variables 
  mutate(
    ch_BPA_conj_i_cor_M2 = ch_BPA_total_i_cor_M2 - ch_BPA_free_i_cor_M2,
    ch_BPA_conj_i_cor_Y1 = ch_BPA_total_i_cor_Y1 - ch_BPA_free_i_cor_Y1,
    ch_MEPA_conj_i_cor_Y1 = ch_MEPA_total_i_cor_Y1 - ch_MEPA_free_i_cor_Y1)


## Problème des valeurs conjugées égales ou < 0 
metadata %>% 
  select(ch_BPA_conj_i_cor_M2, ch_BPA_total_i_cor_M2, ch_BPA_free_i_cor_M2) %>% 
  filter(ch_BPA_total_i_cor_M2 < ch_BPA_free_i_cor_M2) %>% 
  arrange(ch_BPA_conj_i_cor_M2, ch_BPA_total_i_cor_M2, ch_BPA_free_i_cor_M2)

metadata %>% 
  select(ch_BPA_conj_i_cor_Y1, ch_BPA_total_i_cor_Y1, ch_BPA_free_i_cor_Y1) %>% 
  filter(ch_BPA_total_i_cor_Y1 < ch_BPA_free_i_cor_Y1) %>% 
  arrange(ch_BPA_conj_i_cor_Y1, ch_BPA_total_i_cor_Y1, ch_BPA_free_i_cor_Y1)

metadata %>% 
  select(ch_MEPA_conj_i_cor_Y1, ch_MEPA_total_i_cor_Y1, ch_MEPA_free_i_cor_Y1) %>% 
  filter(ch_MEPA_total_i_cor_Y1 < ch_MEPA_free_i_cor_Y1) %>% 
  arrange(ch_MEPA_conj_i_cor_Y1, ch_MEPA_total_i_cor_Y1, ch_MEPA_free_i_cor_Y1)

metadata %>% 
  select(ch_BPA_conj_i_cor_M2, ch_BPA_conj_i_cor_Y1, ch_MEPA_conj_i_cor_Y1) %>% 
  filter(ch_BPA_conj_i_cor_M2 == 0 | ch_BPA_conj_i_cor_Y1 == 0 | ch_MEPA_conj_i_cor_Y1 == 0) %>% 
  arrange(ch_BPA_conj_i_cor_M2, ch_BPA_conj_i_cor_Y1, ch_MEPA_conj_i_cor_Y1)



metadata %>%                  # look at the first value above 0 
  select(ch_BPA_conj_i_cor_M2) %>%              
  filter(ch_BPA_conj_i_cor_M2>0) %>% 
  arrange(ch_BPA_conj_i_cor_M2) %>%
  slice(1)

metadata %>%                  # look at the first value above 0 
  select(ch_BPA_conj_i_cor_Y1) %>%              
  filter(ch_BPA_conj_i_cor_Y1>0) %>% 
  arrange(ch_BPA_conj_i_cor_Y1) %>%
  slice(1)

metadata %>%                  # look at the first value above 0 
  select(ch_MEPA_conj_i_cor_Y1) %>%              
  filter(ch_MEPA_conj_i_cor_Y1>0) %>% 
  arrange(ch_MEPA_conj_i_cor_Y1) %>%
  slice(1)


metadata %>%                  # look at the first value above 0 
  select(ch_BPA_conj_i_cor_M2) %>%              
  filter(ch_BPA_conj_i_cor_M2 < 0.001) %>% 
  arrange(ch_BPA_conj_i_cor_M2) 

metadata %>%                  # look at the first value above 0 
  select(ch_BPA_conj_i_cor_Y1) %>%              
  filter(ch_BPA_conj_i_cor_Y1 <0.01) %>% 
  arrange(ch_BPA_conj_i_cor_Y1)

metadata %>%                  # look at the first value above 0 
  select(ch_MEPA_conj_i_cor_Y1) %>%              
  filter(ch_MEPA_conj_i_cor_Y1 < 0.1) %>% 
  arrange(ch_MEPA_conj_i_cor_Y1) 


# BPA M2 conj: the smaller value above 0 is 0.005133234 --> then 0 values have to be replace by the nearest decimal value : 0.001
# BPA Y1 conj: the smaller value above 0 is 0.02123942 --> then 0 values have to be replace by the nearest decimal value : 0.01
# MEPA Y1 conj: the smaller value above 0 is 0.1665189 --> then 0 values have to be replace by the nearest decimal value : 0.1

metadata <- metadata %>%
  mutate(
    ch_BPA_conj_i_cor_M2_rec = if_else(ch_BPA_conj_i_cor_M2 < 0.001, 0.001, unvr(ch_BPA_conj_i_cor_M2)), 
    ch_BPA_conj_i_cor_Y1_rec = if_else(ch_BPA_conj_i_cor_Y1 <0.01, 0.01, unvr(ch_BPA_conj_i_cor_Y1)),
    ch_MEPA_conj_i_cor_Y1_rec = if_else(ch_MEPA_conj_i_cor_Y1 <0.1, 0.1, unvr(ch_MEPA_conj_i_cor_Y1))) %>%
  mutate(
    ch_BPA_conj_i_cor_M2_ln = log(ch_BPA_conj_i_cor_M2_rec),
    ch_BPA_conj_i_cor_Y1_ln = log(ch_BPA_conj_i_cor_Y1_rec),
    ch_MEPA_conj_i_cor_Y1_ln = log(ch_MEPA_conj_i_cor_Y1_rec),
    ch_BPA_free_i_cor_M2_ln = log(ch_BPA_free_i_cor_M2),
    ch_BPA_free_i_cor_Y1_ln = log(ch_BPA_free_i_cor_Y1),
    ch_MEPA_free_i_cor_Y1_ln = log(ch_MEPA_free_i_cor_Y1)
  ) 


# Nettoyage outlier ident 17673 pour PFAS (PFHpS, PFHxS et PFOS) ----
metadata <- metadata %>%
  mutate(
    mo_PFHpS_i_cor_ln_exclu17673 = ifelse(ident == 17673, NA, mo_PFHpS_i_cor_ln),
    mo_PFHxS_cor_ln_exclu17673 = ifelse(ident == 17673, NA, mo_PFHxS_cor_ln), 
    mo_PFOS_cor_ln_exclu17673 = ifelse(ident == 17673, NA, mo_PFOS_cor_ln))

bdd_outlier_pfas <- bdd_outlier_pfas %>%
  mutate(
    mo_PFHxS_cor_ln_new = log(mo_PFHxS_cor_new), 
    mo_PFOS_cor_ln_new = log(mo_PFOS_cor_new), 
    mo_PFHpS_i_cor_ln_new = log(mo_PFHpS_i_cor_new)
  ) %>%
  mutate(
    mo_PFHpS_i_cor_ln_new_exclu17673 = ifelse(ident == 17673, NA, mo_PFHpS_i_cor_ln_new),
    mo_PFHxS_cor_ln_new_exclu17673 = ifelse(ident == 17673, NA, mo_PFHxS_cor_ln_new), 
    mo_PFOS_cor_ln_new_exclu17673 = ifelse(ident == 17673, NA, mo_PFOS_cor_ln_new))
metadata <- left_join(metadata, 
                      bdd_outlier_pfas, 
                      by="ident")


# Nettoyage pour masses molaires phthalates ----
phthalates_vec_ms <- metadata %>%
  select(mo_DiNP_ms_i_cor_t2, mo_DiNP_ms_i_cor_t3, ch_DiNP_ms_i_cor_M2, ch_DiNP_ms_i_cor_Y1, 
         mo_DEHP_ms_i_cor_t2, mo_DEHP_ms_i_cor_t3, ch_DEHP_ms_i_cor_M2, ch_DEHP_ms_i_cor_Y1,
         mo_DINCH_ms_i_cor_t2, mo_DINCH_ms_i_cor_t3, ch_DINCH_ms_i_cor_Y1) %>% colnames()
phthalates_vec_ms_ln <- paste(phthalates_vec_ms, "ln", sep = "_")

## Création sommes molaires phthalates ln
metadata[, phthalates_vec_ms] <- lapply(metadata[, phthalates_vec_ms], as.numeric) 
## les NAs "introduitsé proviennent des phenols codées en character --> numérique qui avaient des NAs
metadata %>%
  select(mo_DiNP_ms_i_cor_t2, mo_DiNP_ms_i_cor_t3, ch_DiNP_ms_i_cor_M2, ch_DiNP_ms_i_cor_Y1, 
         mo_DEHP_ms_i_cor_t2, mo_DEHP_ms_i_cor_t3, ch_DEHP_ms_i_cor_M2, ch_DEHP_ms_i_cor_Y1,
         mo_DINCH_ms_i_cor_t2, mo_DINCH_ms_i_cor_t3, ch_DINCH_ms_i_cor_Y1) %>% describe()

phthalates_ms_ln <- metadata %>% 
  select(ident, all_of(phthalates_vec_ms)) %>% 
  select(ident, everything())

phthalates_ms_ln[,phthalates_vec_ms] <- log(phthalates_ms_ln[,phthalates_vec_ms])
colnames(phthalates_ms_ln) <- c("ident", paste(colnames(phthalates_ms_ln[, 2:12]), "ln", sep="_")) 
metadata <- left_join(metadata,
                      phthalates_ms_ln, 
                      by ="ident")

## Création sommes molaires phthalates ter
phthalates_ms_ter <- metadata %>% 
  select(ident, all_of(phthalates_vec_ms)) %>% 
  select(ident, everything())
phthalates_ms_ter[, phthalates_vec_ms] <- lapply(
  phthalates_ms_ter[, phthalates_vec_ms],
  quant.cut,
  nbclass = 3,
  include.lowest = TRUE,
  right = FALSE,
  dig.lab = 2)
colnames(phthalates_ms_ter) <- c("ident", paste(colnames(phthalates_ms_ter[, 2:12]), "ter", sep="_")) 

phthalates_ms_ter <- phthalates_ms_ter %>% mutate(
  mo_DiNP_ms_i_cor_t2_ter = fct_recode(
    mo_DiNP_ms_i_cor_t2_ter,
    "1st tertile" = "[0.0075,0.031)",
    "2nd tertile" = "[0.031,0.059)",
    "3rd tertile" = "[0.059,1.2]"
  ),
  mo_DiNP_ms_i_cor_t3_ter = fct_recode(
    mo_DiNP_ms_i_cor_t3_ter,
    "1st tertile" = "[0.0098,0.029)",
    "2nd tertile" = "[0.029,0.051)",
    "3rd tertile" = "[0.051,0.88]"
  ),
  ch_DiNP_ms_i_cor_M2_ter = fct_recode(
    ch_DiNP_ms_i_cor_M2_ter,
    "1st tertile" = "[0.003,0.0078)",
    "2nd tertile" = "[0.0078,0.0097)",
    "3rd tertile" = "[0.0097,0.48]"
  ),
  ch_DiNP_ms_i_cor_Y1_ter = fct_recode(
    ch_DiNP_ms_i_cor_Y1_ter,
    "1st tertile" = "[0.0015,0.02)",
    "2nd tertile" = "[0.02,0.035)",
    "3rd tertile" = "[0.035,0.66]"
  ),
  mo_DEHP_ms_i_cor_t2_ter = fct_recode(
    phthalates_ms_ter$mo_DEHP_ms_i_cor_t2_ter,
    "1st tertile" = "[0.018,0.089)",
    "2nd tertile" = "[0.089,0.14)",
    "3rd tertile" = "[0.14,3.7]"
  ),
  mo_DEHP_ms_i_cor_t3_ter = fct_recode(
    mo_DEHP_ms_i_cor_t3_ter,
    "1st tertile" = "[0.03,0.086)",
    "2nd etrtile" = "[0.086,0.13)",
    "3rd tertile" = "[0.13,1.2]"
  ),
  ch_DEHP_ms_i_cor_M2_ter = fct_recode(
    ch_DEHP_ms_i_cor_M2_ter,
    "1st tertile" = "[0.019,0.031)",
    "2nd tertile" = "[0.031,0.042)",
    "3rd tertile" = "[0.042,0.23]"
  ),
  ch_DEHP_ms_i_cor_Y1_ter = fct_recode(
    ch_DEHP_ms_i_cor_Y1_ter,
    "1st tertile" = "[0.005,0.059)",
    "2nd tertile" = "[0.059,0.11)",
    "3rd tertile" = "[0.11,1.6]"
  ),
  mo_DINCH_ms_i_cor_t2_ter = fct_recode(
    mo_DINCH_ms_i_cor_t2_ter,
    "1st tertile" = "[0.0023,0.0081)",
    "2nd tertile" = "[0.0081,0.014)",
    "3rd tertile" = "[0.014,0.57]"
  ),
  mo_DINCH_ms_i_cor_t3_ter = fct_recode(
    mo_DINCH_ms_i_cor_t3_ter,
    "1st tertile" = "[0.0016,0.0074)",
    "2nd tertile" = "[0.0074,0.013)",
    "3rd tertile" = "[0.013,0.69]"
  ),
  ch_DINCH_ms_i_cor_Y1_ter = fct_recode(
    ch_DINCH_ms_i_cor_Y1_ter,
    "1st tertile" = "[0.00096,0.0079)",
    "2nd tertile" = "[0.0079,0.017)",
    "3rd etrtile" = "[0.017,0.76]"
  )
)
metadata <- left_join(metadata,
                      phthalates_ms_ter, 
                      by ="ident")

# Création des vecteurs ----
## vecteurs phénols ----
phenols_vec <- metadata %>% select(
  mo_MEPA_total_i_cor_t2_ln,
  mo_MEPA_total_i_cor_t3_ln,
  ch_MEPA_total_i_cor_M2_ln,
  ch_MEPA_total_i_cor_Y1_ln,
  
  mo_ETPA_total_i_cor_t2_ln, 
  mo_ETPA_total_i_cor_t3_ln,
  ch_ETPA_total_cat_M2_2,                                      
  ch_ETPA_total_i_cor_Y1_ln,
  
  mo_PRPA_total_i_cor_t2_ln,
  mo_PRPA_total_i_cor_t3_ln,
  ch_PRPA_total_cat_M2_2,                                      
  ch_PRPA_total_i_cor_Y1_ln,
  
  mo_BUPA_total_cat_t2,                                          
  mo_BUPA_total_cat_t3,                                          
  ch_BUPA_total_cat_M2_2,                                      
  ch_BUPA_total_cat_Y1,                                      
  
  mo_BPA_total_i_cor_t2_ln,
  mo_BPA_total_i_cor_t3_ln,
  ch_BPA_total_i_cor_M2_ln,
  ch_BPA_total_i_cor_Y1_ln,
  
  mo_BPS_total_cat_t2_2,                                          
  mo_BPS_total_cat_t3_2,                                          
  ch_BPS_total_cat_M2_2,                                      
  ch_BPS_total_cat_Y1,                                      
  
  mo_OXBE_total_i_cor_t2_ln,
  mo_OXBE_total_i_cor_t3_ln,
  ch_OXBE_total_i_cor_M2_ln,
  ch_OXBE_total_i_cor_Y1_ln,
  
  mo_TRCS_total_i_cor_t2_ln,
  mo_TRCS_total_i_cor_t3_ln,
  ch_TRCS_total_i_cor_M2_ln,                                     
  ch_TRCS_total_i_cor_Y1_ln                                       
) %>% colnames()

phenols_vec_cat <- metadata %>% select(all_of(phenols_vec)) %>% select(contains("cat")) %>% colnames()
phenols_vec_ln  <- metadata %>% select(all_of(phenols_vec)) %>% select(!contains("cat")) %>% colnames()
phenols_vec_ter <- metadata %>% select(all_of(phenols_vec)) %>% select(!contains("cat")) %>% colnames() %>% str_replace("_ln", "_ter")
phenols_vec_num <- metadata %>% select(all_of(phenols_vec)) %>% select(!contains("cat")) %>% colnames() %>% str_replace("_ln", "")
phenols_vec_num_2 <- phenols_vec_num %>% str_replace_all(
  c("ch_BPA_total_i_cor_M2" = "ch_BPA_conj_i_cor_M2", 
    "ch_BPA_total_i_cor_Y1" = "ch_BPA_conj_i_cor_Y1", 
    "ch_MEPA_total_i_cor_Y1" = "ch_MEPA_conj_i_cor_Y1"))

phenols_vec_num_t2 <- metadata %>% select(all_of(phenols_vec_num)) %>% select(contains("t2")) %>% colnames()
phenols_vec_num_t3 <- metadata %>% select(all_of(phenols_vec_num)) %>% select(contains("t3")) %>% colnames()
phenols_vec_num_M2 <- metadata %>% select(all_of(phenols_vec_num)) %>% select(contains("M2")) %>% colnames()
phenols_vec_num_Y1 <- metadata %>% select(all_of(phenols_vec_num)) %>% select(contains("Y1")) %>% colnames()

phenols_vec_num_sg <- phenols_vec_num %>% str_replace_all("total_i_cor", "total_i_cor_sg")
phenols_vec_num_sg_ln <- paste(phenols_vec_num_sg, "ln", sep = "_")
phenols_vec_num_sg_t2 <- metadata %>% select(all_of(phenols_vec_num_sg)) %>% select(contains("t2")) %>% colnames()
phenols_vec_num_sg_t3 <- metadata %>% select(all_of(phenols_vec_num_sg)) %>% select(contains("t3")) %>% colnames()
phenols_vec_num_sg_M2 <- metadata %>% select(all_of(phenols_vec_num_sg)) %>% select(contains("M2")) %>% colnames()
phenols_vec_num_sg_Y1 <- metadata %>% select(all_of(phenols_vec_num_sg)) %>% select(contains("Y1")) %>% colnames()

phenols_vec_2 <- metadata %>% select(all_of(phenols_vec)) %>% colnames() %>% str_replace_all(
  c("ch_BPA_total_i_cor_M2_ln" = "ch_BPA_conj_i_cor_M2_ln", 
    "ch_BPA_total_i_cor_Y1_ln" = "ch_BPA_conj_i_cor_Y1_ln", 
    "ch_MEPA_total_i_cor_Y1_ln" = "ch_MEPA_conj_i_cor_Y1_ln"))



## vecteurs PFAS ----
pfas_vec <- metadata %>% select(mo_PFOA_cor_ln,
                                 mo_PFNA_ln,
                                 mo_PFHxS_cor_ln,
                                 mo_PFOS_cor_ln,
                                 
                                 mo_PFDA_i_cor_ln,
                                 mo_PFUnDA_i_cor_ln, 
                                 mo_PFHpS_i_cor_ln,
                                 
                                 mo_PFDoDa_cat, 
                                 mo_PFHxPA_cat, 
                                 mo_PFHpA_cat_2,
                                 mo_PFTrDa_cat, 
                                 mo_PFBS_cat_2, 
                                 mo_PFOSA_cat_2, 
                                 mo_6_2diPAP_cat_2,
                                 mo_8_2diPAP_cat_2) %>% colnames()

pfas_vec_cat <- metadata %>% select(all_of(pfas_vec)) %>% select(contains("cat")) %>% colnames()
pfas_vec_ln <- metadata %>% select(all_of(pfas_vec)) %>% select(!contains("cat")) %>% colnames()
pfas_vec_ter <- metadata %>% select(all_of(pfas_vec)) %>% select(!contains("cat")) %>% colnames() %>% str_replace("_ln", "_ter")
pfas_vec_num <- metadata %>% select(all_of(pfas_vec)) %>% select(!contains("cat")) %>% colnames() %>% str_replace("_ln", "")

pfas_vec_2 <- pfas_vec %>% str_replace_all(
  c("mo_PFHpS_i_cor_ln" = "mo_PFHpS_i_cor_ln_exclu17673", 
    "mo_PFHxS_cor_ln" = "mo_PFHxS_cor_ln_exclu17673", 
    "mo_PFOS_cor_ln" = "mo_PFOS_cor_ln_exclu17673"))
pfas_vec_3 <- pfas_vec %>% str_replace_all(
  c("mo_PFHpS_i_cor_ln" = "mo_PFHpS_i_cor_ln_new", 
    "mo_PFHxS_cor_ln" = "mo_PFHxS_cor_ln_new", 
    "mo_PFOS_cor_ln" = "mo_PFOS_cor_ln_new"))
pfas_vec_exclu17673 <- c("mo_PFHxS_cor_ln_exclu17673", 
                         "mo_PFOS_cor_ln_exclu17673", 
                         "mo_PFHpS_i_cor_ln_exclu17673", 
                         
                         "mo_PFOS_cor_ln_new_exclu17673", 
                         "mo_PFHxS_cor_ln_new_exclu17673", 
                         "mo_PFHpS_i_cor_ln_new_exclu17673")




## vecteur phthalates ----
phthalates_vec <- metadata %>% select(
  mo_ohMiNP_i_cor_t2_ln, mo_ohMiNP_i_cor_t3_ln, ch_ohMiNP_i_cor_M2_ln, ch_ohMiNP_i_cor_Y1_ln, 
  mo_oxoMiNP_i_cor_t2_ln, mo_oxoMiNP_i_cor_t3_ln, ch_oxoMiNP_i_cor_M2_ln, ch_oxoMiNP_i_cor_Y1_ln, 
  mo_cxMiNP_i_cor_t2_ln, mo_cxMiNP_i_cor_t3_ln, ch_cxMiNP_i_cor_M2_ln, ch_cxMiNP_i_cor_Y1_ln, 
  mo_DiNP_ms_i_cor_t2_ln, mo_DiNP_ms_i_cor_t3_ln, ch_DiNP_ms_i_cor_M2_ln, ch_DiNP_ms_i_cor_Y1_ln,   # Métabolites DiNP
  
  mo_MEOHP_i_cor_t2_ln, mo_MEOHP_i_cor_t3_ln, ch_MEOHP_i_cor_M2_ln, ch_MEOHP_i_cor_Y1_ln, 
  mo_MECPP_i_cor_t2_ln, mo_MECPP_i_cor_t3_ln, ch_MECPP_i_cor_M2_ln, ch_MECPP_i_cor_Y1_ln, 
  mo_MEHHP_i_cor_t2_ln, mo_MEHHP_i_cor_t3_ln, ch_MEHHP_i_cor_M2_ln, ch_MEHHP_i_cor_Y1_ln,    
  mo_MEHP_i_cor_t2_ln, mo_MEHP_i_cor_t3_ln, ch_MEHP_i_cor_M2_ln, ch_MEHP_i_cor_Y1_ln, 
  mo_MMCHP_i_cor_t2_ln, mo_MMCHP_i_cor_t3_ln, ch_MMCHP_i_cor_M2_ln, ch_MMCHP_i_cor_Y1_ln, 
  mo_DEHP_ms_i_cor_t2_ln, mo_DEHP_ms_i_cor_t3_ln, ch_DEHP_ms_i_cor_M2_ln, ch_DEHP_ms_i_cor_Y1_ln,   # Métabolites DEHP
  
  mo_MnBP_i_cor_t2_ln, mo_MnBP_i_cor_t3_ln, ch_MnBP_i_cor_M2_ln, ch_MnBP_i_cor_Y1_ln,               # Métabolite DBP
  mo_MiBP_i_cor_t2_ln, mo_MiBP_i_cor_t3_ln, ch_MiBP_i_cor_M2_ln, ch_MiBP_i_cor_Y1_ln,               # Métabolite DiBP
  mo_MBzP_i_cor_t2_ln, mo_MBzP_i_cor_t3_ln, ch_MBzP_i_cor_M2_ln, ch_MBzP_i_cor_Y1_ln,               # Metabolite MBzP
  mo_MEP_i_cor_t2_ln, mo_MEP_i_cor_t3_ln, ch_MEP_i_cor_M2_ln, ch_MEP_i_cor_Y1_ln,                   # Metabolite DEP
  mo_ohMPHP_i_cor_t2_ln, mo_ohMPHP_i_cor_t3_ln, ch_ohMPHP_cat_M2_2, ch_ohMPHP_i_cor_Y1_ln,          # Metabolite DEHP
  
  mo_ohMINCH_i_cor_t2_ln, mo_ohMINCH_i_cor_t3_ln, ch_ohMINCH_cat_M2_2, ch_ohMINCH_i_cor_Y1_ln, 
  mo_oxoMINCH_i_cor_t2_ln, mo_oxoMINCH_i_cor_t3_ln, ch_oxoMINCH_cat_M2_2, ch_oxoMINCH_i_cor_Y1_ln,
  mo_DINCH_ms_i_cor_t2_ln, mo_DINCH_ms_i_cor_t3_ln, ch_ohMINCH_cat_M2_2, ch_oxoMINCH_cat_M2_2, ch_DINCH_ms_i_cor_Y1_ln # Metabolites DINCH
) %>% colnames()


phthalates_vec_cat <- metadata %>% select(all_of(phthalates_vec)) %>% select(contains("cat")) %>% colnames()
phthalates_vec_ln  <- metadata %>% select(all_of(phthalates_vec)) %>% select(!contains("cat")) %>% colnames()
phthalates_vec_ter <- metadata %>% select(all_of(phthalates_vec)) %>% select(!contains("cat")) %>% colnames() %>% str_replace("_ln", "_ter")

phthalates_vec_num <- metadata %>% select(all_of(phthalates_vec)) %>% select(!contains("cat")) %>% colnames() %>% str_replace("_ln", "")
phthalates_vec_num_t2 <- metadata %>% select(all_of(phthalates_vec_num)) %>% select(contains("t2")) %>% colnames()
phthalates_vec_num_t3 <- metadata %>% select(all_of(phthalates_vec_num)) %>% select(contains("t3")) %>% colnames()
phthalates_vec_num_M2 <- metadata %>% select(all_of(phthalates_vec_num)) %>% select(contains("M2")) %>% colnames()
phthalates_vec_num_Y1 <- metadata %>% select(all_of(phthalates_vec_num)) %>% select(contains("Y1")) %>% colnames()


phthalates_vec_num_sg <- phthalates_vec_num %>% str_replace_all("_i_cor", "_i_cor_sg")
phthalates_vec_num_sg_t2 <- metadata %>% select(all_of(phthalates_vec_num_sg)) %>% select(contains("t2")) %>% colnames()
phthalates_vec_num_sg_t3 <- metadata %>% select(all_of(phthalates_vec_num_sg)) %>% select(contains("t3")) %>% colnames()
phthalates_vec_num_sg_M2 <- metadata %>% select(all_of(phthalates_vec_num_sg)) %>% select(contains("M2")) %>% colnames()
phthalates_vec_num_sg_Y1 <- metadata %>% select(all_of(phthalates_vec_num_sg)) %>% select(contains("Y1")) %>% colnames()

phthalates_vec_num_sg_ln <- paste(phthalates_vec_num_sg, "ln", sep = "_")

## vecteur covariables ----
covar_vec_i <- metadata %>% select(                       
  "ch_feces_RUN_Y1",                   
  "ch_feces_age_w_Y1_i",
  "po_delmod",
  "ch_food_intro_Y1_3cat_i",
  "ch_antibio_Y1_2cat_i",
  "mo_par_2cat",
  "mo_pets_i",
  "ch_sex",
  "mo_tob_gr_anyt_yn_n2_i",
  "Mo_ETS_anyT_yn1_opt_i",
  "ch_ETS_12m_opt36m",
  "mo_interpreg_3cat",
  "mo_dipl_3cat_i",
  "po_w_kg_3cat",
  "po_he_3cat_i",
  "ch_w_Y1_3cat_i",
  "ch_he_Y1_3cat_i",
  "po_gd",
  "mo_age",
  "mo_bmi_bepr_3cat_i",
  "bf_duration_till48w_4cat_i"
) %>% colnames()

covar_vec <- metadata %>% select(                       
  "ch_feces_RUN_Y1",                   
  "ch_feces_age_w_Y1",
  "po_delmod",
  "ch_food_intro_Y1_3cat",
  "ch_antibio_Y1_2cat",
  "mo_par_2cat",
  "mo_pets",
  "ch_sex",
  "mo_tob_gr_anyt_yn_n2",
  "Mo_ETS_anyT_yn1_opt",
  "ch_ETS_12m_opt36m",
  "mo_interpreg_3cat",
  "mo_dipl_3cat",
  "po_w_kg_3cat",
  "po_he_3cat",
  "ch_w_Y1_3cat",
  "ch_he_Y1_3cat",
  "po_gd",
  "mo_age",
  "mo_bmi_bepr_3cat",
  "bf_duration_till48w_4cat"
) %>% colnames()

covar_vec_num_i <- metadata %>% select(                       
  "ch_feces_age_w_Y1_i",
  "ch_antibio_Y1_i",
  "mo_par",
  "po_w_kg",
  "po_he_i",
  "ch_w_Y1_i",
  "ch_he_Y1_i",
  "po_gd",
  "mo_age",
  "mo_bmi_bepr_i",
  "bf_duration_till48w_i"
) %>% colnames()

covar_vec_num <- metadata %>% select(                       
  "ch_feces_age_w_Y1",
  "ch_antibio_Y1",
  "mo_par",
  "po_w_kg",
  "po_he",
  "ch_w_Y1",
  "ch_he_Y1",
  "po_gd",
  "mo_age",
  "mo_bmi_bepr",
  "bf_duration_till48w"
) %>% colnames()

covar_vec_cat_i <- metadata %>% select(                       
  "ch_feces_RUN_Y1",    
  "po_delmod",
  "ch_food_intro_Y1_3cat_i",
  "ch_antibio_Y1_2cat_i",
  "mo_par_2cat",
  "mo_pets_i",
  "ch_sex",
  "mo_tob_gr_anyt_yn_n2_i",
  "Mo_ETS_anyT_yn1_opt_i",
  "ch_ETS_12m_opt36m",
  "mo_interpreg_3cat",
  "mo_dipl_3cat_i",
  "po_w_kg_3cat",
  "po_he_3cat_i",
  "ch_w_Y1_3cat_i",
  "ch_he_Y1_3cat_i",
  "mo_bmi_bepr_3cat_i",
  "bf_duration_till48w_4cat_i"
) %>% colnames()

covar_vec_cat <- metadata %>% select(                       
  "ch_feces_RUN_Y1",    
  "po_delmod",
  "ch_food_intro_Y1_3cat",
  "ch_antibio_Y1_2cat",
  "mo_par_2cat",
  "mo_pets",
  "ch_sex",
  "mo_tob_gr_anyt_yn_n2",
  "Mo_ETS_anyT_yn1_opt",
  "ch_ETS_12m_opt36m",
  "mo_interpreg_3cat",
  "mo_dipl_3cat",
  "po_w_kg_3cat",
  "po_he_3cat",
  "ch_w_Y1_3cat",
  "ch_he_Y1_3cat",
  "mo_bmi_bepr_3cat",
  "bf_duration_till48w_4cat"
) %>% colnames()




## vecteur outcome ----
alpha_vec <- alphadiv_Y1 %>% select(-ident, -ch_feces_ID_Y1) %>% colnames()
taxa_vec <- c("ch_feces_rel_p1_Y1", "ch_feces_rel_p2_Y1", "ch_feces_rel_p3_Y1", "ch_feces_rel_p4_Y1", 
              "ch_feces_rel_g1_Y1", "ch_feces_rel_g2_Y1", "ch_feces_rel_g3_Y1", "ch_feces_rel_g4_Y1")





# Nettoyage gravité spécifique variables Matthieu xx_comp_i_cor_sg_xx ----
## Phénols ----
metadata[, phenols_vec_num_sg] <- lapply(metadata[, phenols_vec_num_sg], as.numeric)  # les NAs "introduitsé proviennent des phenols codées en character --> numérique qui avaient des NAs
phenols_ln_sg <- metadata %>% 
  select(ident, all_of(phenols_vec_num_sg))%>% 
  select(ident, everything())
phenols_ln_sg[,phenols_vec_num_sg] <- log(phenols_ln_sg[,phenols_vec_num_sg])
colnames(phenols_ln_sg) <- c("ident", paste(colnames(phenols_ln_sg[, 2:23]), "ln", sep="_")) 
metadata <- left_join(metadata, 
                       phenols_ln_sg, 
                       by = "ident")



## Phthalates ----
metadata[, phthalates_vec_num_sg] <- lapply(metadata[, phthalates_vec_num_sg], as.numeric) # les NAs introduit proviennent des phthalates codées en character --> numérique qui avaient des NAs
phthalates_ln_sg <- metadata %>% 
  select(ident, all_of(phthalates_vec_num_sg))%>% 
  select(ident, everything())
phthalates_ln_sg[,phthalates_vec_num_sg] <- log(phthalates_ln_sg[,phthalates_vec_num_sg])
colnames(phthalates_ln_sg) <- c("ident", paste(colnames(phthalates_ln_sg[, 2:69]), "ln", sep="_")) 
metadata <- left_join(metadata, 
                       phthalates_ln_sg, 
                       by = "ident")




# Préparation des données ----
## Classes des variables ----
metadata[, phenols_vec_cat] <- lapply(metadata[, phenols_vec_cat], as.factor) 
metadata[, pfas_vec_cat] <- lapply(metadata[, pfas_vec_cat], as.factor) 
metadata[, phthalates_vec_cat] <- lapply(metadata[, phthalates_vec_cat], as.factor) 
metadata[, covar_vec_cat] <- lapply(metadata[, covar_vec_cat], as.factor) 

## Ordre des catégories des variables catégorielles ----
metadata <- metadata %>% mutate(
  ch_food_intro_Y1_3cat = fct_relevel(
    ch_food_intro_Y1_3cat,
    "Between 0 and 6 months old", "Between 6 and 12 months old", "Not introduced at 12 months old"),
  ch_food_intro_Y1_3cat_i = fct_relevel(
    ch_food_intro_Y1_3cat_i,
    "Between 0 and 6 months old", "Between 6 and 12 months old", "Not introduced at 12 months old"),
  
  mo_par_2cat = fct_relevel(
    mo_par_2cat,
    "None", "1 child or more"),
  mo_interpreg_3cat = fct_relevel(
    mo_interpreg_3cat,
    "Under 2 years", "2 years and more", "Primiparous"),
  
  mo_dipl_3cat = fct_relevel(
    mo_dipl_3cat,
    "2years or less after graduation", "3-4years after graduation", ">=5years after graduation"),
  mo_dipl_3cat_i = fct_relevel(
    mo_dipl_3cat_i,
    "2years or less after graduation", "3-4years after graduation", ">=5years after graduation"),
  
  po_w_kg_3cat = fct_relevel(
    po_w_kg_3cat,
    "<3 Kg", "3-3.4 Kg", ">= 3.5 Kg"),
  po_he_3cat = fct_relevel(
    po_he_3cat,
    "<50 cm", "50-51 cm", ">= 52 cm"),
  ch_w_Y1_3cat = fct_relevel(
    ch_w_Y1_3cat,
    "<8.5 Kg", "8.5-9.9 Kg", ">=10 Kg"),
  ch_he_Y1_3cat = fct_relevel(
    ch_he_Y1_3cat_i,
    "<75 cm", "75-77.9 cm", ">=78 cm"),
  
  po_he_3cat_i = fct_relevel(
    po_he_3cat_i,
    "<50 cm", "50-51 cm", ">= 52 cm"),
  ch_w_Y1_3cat_i = fct_relevel(
    ch_w_Y1_3cat_i,
    "<8.5 Kg", "8.5-9.9 Kg", ">=10 Kg"),
  ch_he_Y1_3cat_i = fct_relevel(
    ch_he_Y1_3cat_i,
    "<75 cm", "75-77.9 cm", ">=78 cm"),
  
  mo_bmi_bepr_3cat = fct_relevel(
    mo_bmi_bepr_3cat,
    "<19 Kg/m2", "19-23.9 Kg/m2", ">=24 Kg/m2"),
  mo_bmi_bepr_3cat_i = fct_relevel(
    mo_bmi_bepr_3cat_i,
    "<19 Kg/m2", "19-23.9 Kg/m2", ">=24 Kg/m2"),
  
  bf_duration_till48w_4cat = fct_relevel(
    bf_duration_till48w_4cat,
    "Not breastfed", "<24 weeks", "24-47 weeks", "Still breastfeed at 48 weeks"), 
  bf_duration_till48w_4cat_i = fct_relevel(
    bf_duration_till48w_4cat_i,
    "Not breastfed", "<24 weeks", "24-47 weeks", "Still breastfeed at 48 weeks"))



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
      "<LOD", "LOD-LOQ", ">LOQ"), 
    
    mo_PFDoDa_cat = fct_relevel(
      mo_PFDoDa_cat,
      "<LOD", "LOD-LOQ", ">LOQ"),
    mo_PFHxPA_cat = fct_relevel(
      mo_PFHxPA_cat,
      "<LOD", "LOD-LOQ", ">LOQ"),
    mo_PFTrDa_cat = fct_relevel(
      mo_PFTrDa_cat,
      "<LOD", "LOD-LOQ", ">LOQ"))

## Labels de variables ----
metadata = modify(metadata,{
  var_lab(ident) = "SEPAGES dentity"
  var_lab(ch_feces_ID_Y1) = "Child age at the one-year stool collection"
  
  var_lab(ch_feces_age_w_Y1) = "Child age at stool collection (weeks)"    # variables dispo continues et catégorielles
  var_lab(ch_feces_age_w_Y1_4cat) = "Child age at stool collection"
  
  var_lab(ch_feces_age_w_Y1_i) = "Child age at stool collection (weeks)"   
  var_lab(ch_feces_age_w_Y1_4cat_i) = "Child age at stool collection"
  
  var_lab(po_gd) = "Gestational duration, completed weeks"    
  var_lab(po_gd_4cat) = "Gestational term, completed weeks" 
  
  var_lab(mo_age) = "Maternal age at conception (years)"                               
  var_lab(mo_age_4cat) = "Maternal age at conception"
  
  var_lab(mo_bmi_bepr) = "Maternal BMI before pregnancy (kg/m2)" 
  var_lab(mo_bmi_bepr_3cat) = "Maternal BMI before pregnancy" 
  var_lab(mo_bmi_bepr_i) = "Maternal BMI before pregnancy (kg/m2)" 
  var_lab(mo_bmi_bepr_3cat_i) = "Maternal BMI before pregnancy" 
  
  var_lab(mo_par) = "Maternal parity" 
  var_lab(mo_par_2cat) = "Maternal parity" 
  
  var_lab(po_w) = "Birth weight (g)" 
  var_lab(po_w_kg) = "Birth weight (Kg)"
  var_lab(po_w_kg_3cat) = "Birth weight" 
  
  var_lab(po_he) = "Birth length (cm)" 
  var_lab(po_he_3cat) = "Birth length" 
  var_lab(po_he_i) = "Birth length (cm)" 
  var_lab(po_he_3cat_i) = "Birth length" 
  
  var_lab(ch_w_Y1) = "Weight at one year (Kg)" 
  var_lab(ch_he_Y1) = "Length at one year (cm)"
  var_lab(ch_w_Y1_i) = "Weight at one year (Kg)" 
  var_lab(ch_he_Y1_i) = "Length at one year (cm)"
  
  var_lab(ch_w_Y1_3cat) = "Weight at one year"
  var_lab(ch_he_Y1_3cat) = "Length at one year"
  var_lab(ch_w_Y1_3cat_i) = "Weight at one year"
  var_lab(ch_he_Y1_3cat_i) = "Length at one year"
  
  var_lab(bf_duration_till48w) = "Breastfeeding duration (weeks)"
  var_lab(bf_duration_till48w_4cat) = "Breastfeeding duration"
  var_lab(bf_duration_till48w_i) = "Breastfeeding duration (weeks)"
  var_lab(bf_duration_till48w_4cat_i) = "Breastfeeding duration"
  
  var_lab(ch_antibio_Y1) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_3cat) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_2cat) = "Antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_i) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_3cat_i) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_2cat_i) = "Antibiotics use between 0-12 months old"
  
  var_lab(po_delmod) = "Delivery mode"              # variables forcément catégorielles 
  var_lab(ch_sex) = "Child sex"
  var_lab(mo_dipl) = "Maternal education" 
  var_lab(mo_dipl_3cat) = "Maternal education"
  var_lab(mo_dipl_3cat_i) = "Maternal education"
  var_lab(mo_pets) = "Presence of pets"  
  var_lab(mo_pets_i) = "Presence of pets"  
  var_lab(ch_food_intro_Y1) = "Period of introduction of solid food 0-1Y"
  var_lab(ch_food_intro_Y1_3cat) = "Period of introduction of solid food 0-1Y"
  var_lab(ch_food_intro_Y1_i) = "Period of introduction of solid food 0-1Y"
  var_lab(ch_food_intro_Y1_3cat_i) = "Period of introduction of solid food 0-1Y"
  var_lab(mo_tob_gr_anyt_yn_n2) = "Maternal active smoking during pregnancy"
  var_lab(Mo_ETS_anyT_yn1_opt) = "Maternal passive smoking during pregnancy"
  var_lab(ch_ETS_12m_opt36m) = "Child passive smoking during pregnancy"
  var_lab(mo_tob_gr_anyt_yn_n2_i) = "Maternal active smoking during pregnancy"
  var_lab(Mo_ETS_anyT_yn1_opt_i) = "Maternal passive smoking during pregnancy"
  var_lab(mo_interpreg) = "Interpregnancy interval (years)"
  var_lab(mo_interpreg_5cat) = "Interpregnancy interval"
  var_lab(mo_interpreg_3cat) = "Interpregnancy interval"
  
  var_lab(mo_ethnicity) = "Maternal ethnicity"                                    # varibales que l'on utilisera pas 
  var_lab(ch_sibling) = "Number of other children living at home (sibling or not)" 
  var_lab(ch_sibling_3cat) = "Number of other children living at home (sibling or not)" 
  var_lab(fa_dipl) = "Paternal education"  
  var_lab(mo_cats) = "Presence of cats" 
  var_lab(mo_dogs) = "Presence of dogs"
  var_lab(mo_birds) = "Presence of birds" 
  var_lab(mo_rodents) = "Presence of rodents" 
  var_lab(mo_other_pets) = "Presence of other pets"
})




var_label(metadata[, phenols_vec]) <- colnames(metadata[, phenols_vec]) %>%
  str_replace_all(
    c(
      "mo_" = "Maternal exposure to ",
      "ch_" = "Child exposure to ",
      "MEPA" = "methyparaben",
      "ETPA" = "ethyparaben",
      "PRPA" = "propylparaben",
      "BUPA" = "butylparaben",
      "BPA" = "bisphenol A",
      "BPS" = "bisphenol S",
      "OXBE" = "benzophenone 3",
      "TRCS" = "triclosan",
      "_total_i_cor_" = " ", 
      "_total_cat_" = " ",
      "t2" = "at trim.2",
      "t3" = "at trim.3",
      "M2" = "at 2 months",
      "Y1" = "at 12 months", 
      "_ter" = "", 
      "_ln" = "", 
      "_2" = ""
    )
  )



var_label(metadata[, pfas_vec]) <- colnames(metadata[, pfas_vec]) %>%
  str_replace_all(
    c(
      "mo_" = "Maternal exposure to ",
      "ch_" = "Child exposure to ",
      "_i_cor_ln" = "", 
      "_cor_ln" = "", 
      "_ln" = "",
      "_i_cor_ter" = "", 
      "_cat_2" = "",
      "_cat" = ""
    )
  )





var_label(metadata[, phthalates_vec]) <- colnames(metadata[, phthalates_vec]) %>%
  str_replace_all(
    c(
      "mo_" = "Maternal exposure to ",
      "ch_" = "Child exposure to ",
      "_i_cor_" = " ", 
      "_cat_" = " ",
      "t2" = "at trim.2",
      "t3" = "at trim.3",
      "M2" = "at 2 months",
      "Y1" = "at 12 months", 
      "_ter" = "", 
      "_ln" = "", 
      "_2" = ""
    )
  )


# Vérification codage ----
metadata %>% filter(statut == "inclu") %>% select(all_of(covar_vec))  %>% tbl_summary() 
metadata %>% filter(statut == "inclu") %>% select(all_of(phenols_vec))  %>% tbl_summary() 
metadata %>% filter(statut == "inclu") %>% select(all_of(pfas_vec))  %>% tbl_summary() 
metadata %>% filter(statut == "inclu") %>% select(all_of(phthalates_vec))  %>% tbl_summary() 


# Regroupement des bases de données ----
var_lab(metadata$ident) <- NULL
var_lab(metadata$ch_feces_ID_Y1) <- NULL
var_lab(alphadiv_Y1$ident) <- NULL
var_lab(alphadiv_Y1$ch_feces_ID_Y1) <- NULL
var_lab(asv_taxa$ident) <- NULL
var_lab(asv_taxa$ch_feces_ID_Y1) <- NULL
bdd_alpha <- left_join(alphadiv_Y1,
                       metadata[, c(
                         "ident",
                         "ch_feces_ID_Y1",
                         phenols_vec, 
                         phenols_vec_2, 
                         phenols_vec_ter,
                         phenols_vec_num_sg,
                         phenols_vec_num_sg_ln,
                         pfas_vec,
                         pfas_vec_2,
                         pfas_vec_3,
                         pfas_vec_ter,
                         pfas_vec_exclu17673,
                         phthalates_vec,
                         phthalates_vec_ter,
                         phthalates_vec_num_sg,
                         phthalates_vec_num_sg_ln,
                         covar_vec, 
                         covar_vec_i
                       )],
                       by = c("ch_feces_ID_Y1", "ident"),
                       all.y = TRUE)
bdd_alpha <- left_join(bdd_alpha,
                       bdd_sg,
                       by = "ident")

bdd_taxa <- left_join(asv_taxa,
                      metadata[, c(
                        "ident",
                        "ch_feces_ID_Y1",
                        phenols_vec, 
                        phenols_vec_2, 
                        phenols_vec_ter,
                        phenols_vec_num_sg,
                        phenols_vec_num_sg_ln,
                        pfas_vec,
                        pfas_vec_2,
                        pfas_vec_3,
                        pfas_vec_ter,
                        pfas_vec_exclu17673,
                        phthalates_vec,
                        phthalates_vec_ter,
                        phthalates_vec_num_sg,
                        phthalates_vec_num_sg_ln,
                        covar_vec,
                        covar_vec_i
                      )],
                      by = c("ch_feces_ID_Y1", "ident"),
                      all.y = TRUE)
bdd_taxa <- left_join(bdd_taxa,
                      bdd_sg,
                      by = "ident")

