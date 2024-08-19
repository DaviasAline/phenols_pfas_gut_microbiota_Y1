## Aline Davias
## Review 2 ES&T
## 28/05/2024


# Chargement des packages 
library(tidyverse)
library(vegan)
library(devtools)
library(pairwiseAdonis)

# 1 : tester la fonction pairwise.adonis ----
load("4_output/betadiv/betadiv.RData")
# on teste la fonction pairwise.adonis seulement pour les 3 expos qui étaient associées 

#___________________________________________________________________________________________________________
# Conclusion: sous réserve d'utiliser la même seed, le même nombre de permutation, les mêmes covariables 
# et le même groupe de référence (opposé), on obtient bien les mêmes résultats entre les 2 fonctions 
#___________________________________________________________________________________________________________

## OBXE M2 ----
set.seed(1996)
test_oxbe_M2_2 <- pairwise.adonis2(      # test de la suggestion du reviewer d'utiliser pairwise.adonis2
  all_dist_M2 ~
    ch_OXBE_total_i_cor_M2_ter +
    ch_feces_RUN_Y1 +
    ch_feces_age_w_Y1_i +
    po_delmod + ch_food_intro_Y1_3cat_i +
    ch_antibio_Y1_2cat_i +
    mo_par_2cat +
    mo_pets_i +
    ch_sex +
    mo_tob_gr_anyt_yn_n2_i +
    Mo_ETS_anyT_yn1_opt_i +
    ch_ETS_12m_opt36m +
    mo_dipl_3cat_i +
    po_w_kg_3cat +
    po_he_3cat_i +
    ch_w_Y1_3cat_i +
    ch_he_Y1_3cat_i +
    po_gd +
    mo_age +
    mo_bmi_bepr_3cat_i +
    bf_duration_till48w_4cat_i,
  data = bdd_final_M2,
  permutations = 999)


bdd_oxbe_M2_1_2 <- bdd_final_M2 %>%   # ancien code manuel 
  filter(ch_OXBE_total_i_cor_M2_ter %in% c("1st tertile", "2nd tertile"))
oxbe_to_keep_1_2 <- bdd_oxbe_M2_1_2$ident
bdd_oxbe_M2_1_2 <- bdd_oxbe_M2_1_2 %>%
  select(ident, all_of(covariables), ch_OXBE_total_i_cor_M2_ter, all_of(oxbe_to_keep_1_2))
all_dist_oxbe_M2_1_2 <- bdd_oxbe_M2_1_2 %>% select(all_of(oxbe_to_keep_1_2)) %>% as.dist()
set.seed(1996)
results_betadiv_multivar_bray_curtis_oxbe_M2_1st_2nd_ter <- 
  adonis2(all_dist_oxbe_M2_1_2 ~ 
            ch_OXBE_total_i_cor_M2_ter + 
            ch_feces_RUN_Y1 + 
            ch_feces_age_w_Y1_i + 
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i + 
            mo_par_2cat +
            mo_pets_i +
            ch_sex + 
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i, 
          data = bdd_oxbe_M2_1_2, 
          permutations = 999)

bdd_oxbe_M2_1_3 <- bdd_final_M2 %>%   # ancien code manuel 
  filter(ch_OXBE_total_i_cor_M2_ter %in% c("1st tertile", "3rd tertile"))
oxbe_to_keep_1_3 <- bdd_oxbe_M2_1_3$ident
bdd_oxbe_M2_1_3 <- bdd_oxbe_M2_1_3 %>%
  select(ident, all_of(covariables), ch_OXBE_total_i_cor_M2_ter, all_of(oxbe_to_keep_1_3))
all_dist_oxbe_M2_1_3 <- bdd_oxbe_M2_1_3 %>% select(all_of(oxbe_to_keep_1_3)) %>% as.dist()
set.seed(1996)
results_betadiv_multivar_bray_curtis_oxbe_M2_1st_3rd_ter <- 
  adonis2(all_dist_oxbe_M2_1_3 ~ 
            ch_OXBE_total_i_cor_M2_ter + 
            ch_feces_RUN_Y1 + 
            ch_feces_age_w_Y1_i + 
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i + 
            mo_par_2cat +
            mo_pets_i +
            ch_sex + 
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i, 
          data = bdd_oxbe_M2_1_3, 
          permutations = 999)


bdd_oxbe_M2_2_3 <- bdd_final_M2 %>%   # ancien code manuel 
  filter(ch_OXBE_total_i_cor_M2_ter %in% c("2nd tertile", "3rd tertile"))
oxbe_to_keep_2_3 <- bdd_oxbe_M2_2_3$ident
bdd_oxbe_M2_2_3 <- bdd_oxbe_M2_2_3 %>%
  select(ident, all_of(covariables), ch_OXBE_total_i_cor_M2_ter, all_of(oxbe_to_keep_2_3))
all_dist_oxbe_M2_2_3 <- bdd_oxbe_M2_2_3 %>% select(all_of(oxbe_to_keep_2_3)) %>% as.dist()
set.seed(1996)
results_betadiv_multivar_bray_curtis_oxbe_M2_2nd_3rd_ter <- 
  adonis2(all_dist_oxbe_M2_2_3 ~ 
            ch_OXBE_total_i_cor_M2_ter + 
            ch_feces_RUN_Y1 + 
            ch_feces_age_w_Y1_i + 
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i + 
            mo_par_2cat +
            mo_pets_i +
            ch_sex + 
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i, 
          data = bdd_oxbe_M2_2_3, 
          permutations = 999)


## MEPA Y1 ----
test_mepa_Y1 <-
  pairwise.adonis2(
    all_dist_Y1 ~
      ch_MEPA_total_i_cor_Y1_ter +
      ch_feces_RUN_Y1 +
      ch_feces_age_w_Y1_i +
      po_delmod + ch_food_intro_Y1_3cat_i +
      ch_antibio_Y1_2cat_i +
      mo_par_2cat +
      mo_pets_i +
      ch_sex +
      mo_tob_gr_anyt_yn_n2_i +
      Mo_ETS_anyT_yn1_opt_i +
      ch_ETS_12m_opt36m +
      mo_dipl_3cat_i +
      po_w_kg_3cat +
      po_he_3cat_i +
      ch_w_Y1_3cat_i +
      ch_he_Y1_3cat_i +
      po_gd +
      mo_age +
      mo_bmi_bepr_3cat_i +
      bf_duration_till48w_4cat_i,
    data = bdd_final_Y1,
    permutations = 999
  )

## PRPA Y1 ----
test_prpa_Y1 <-
  pairwise.adonis2(
    all_dist_Y1 ~
      ch_PRPA_total_i_cor_Y1_ter +
      ch_feces_RUN_Y1 +
      ch_feces_age_w_Y1_i +
      po_delmod + ch_food_intro_Y1_3cat_i +
      ch_antibio_Y1_2cat_i +
      mo_par_2cat +
      mo_pets_i +
      ch_sex +
      mo_tob_gr_anyt_yn_n2_i +
      Mo_ETS_anyT_yn1_opt_i +
      ch_ETS_12m_opt36m +
      mo_dipl_3cat_i +
      po_w_kg_3cat +
      po_he_3cat_i +
      ch_w_Y1_3cat_i +
      ch_he_Y1_3cat_i +
      po_gd +
      mo_age +
      mo_bmi_bepr_3cat_i +
      bf_duration_till48w_4cat_i,
    data = bdd_final_Y1,
    permutations = 999
  )

rm(list = ls())



# 2 : tester la fonction MaAsLin ----
library(tidyverse)
library(expss)
library(gtsummary)
library(Maaslin2)
library(labelled)
library(questionr)
library(sjlabelled)
library(openxlsx)
library(corrplot)
library(see)
library(psych)
library(compositions)
library(writexl)
library(broom)


## Vectors ----
covariates <- c("ch_feces_RUN_Y1",                   
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
                "bf_duration_till48w_4cat_i")
pollutants <- c("mo_MEPA_total_i_cor_t2_ln",
                "mo_MEPA_total_i_cor_t3_ln",
                "ch_MEPA_total_i_cor_M2_ln",
                "ch_MEPA_conj_i_cor_Y1_ln",   # conjugated concentration in the main analysis
                
                "mo_ETPA_total_i_cor_t2_ln", 
                "mo_ETPA_total_i_cor_t3_ln",
                "ch_ETPA_total_cat_M2_2",                                      
                "ch_ETPA_total_i_cor_Y1_ln",
                
                "mo_PRPA_total_i_cor_t2_ln",
                "mo_PRPA_total_i_cor_t3_ln",
                "ch_PRPA_total_cat_M2_2",                                      
                "ch_PRPA_total_i_cor_Y1_ln",
                
                "mo_BUPA_total_cat_t2",                                          
                "mo_BUPA_total_cat_t3",                                          
                "ch_BUPA_total_cat_M2_2",                                      
                "ch_BUPA_total_cat_Y1",                                      
                
                "mo_BPA_total_i_cor_t2_ln",
                "mo_BPA_total_i_cor_t3_ln",
                "ch_BPA_conj_i_cor_M2_ln",     # conjugated concentration in the main analysis
                "ch_BPA_conj_i_cor_Y1_ln",     # conjugated concentration in the main analysis
                
                "mo_BPS_total_cat_t2_2",                                          
                "mo_BPS_total_cat_t3_2",                                          
                "ch_BPS_total_cat_M2_2",                                      
                "ch_BPS_total_cat_Y1",                                      
                
                "mo_OXBE_total_i_cor_t2_ln",
                "mo_OXBE_total_i_cor_t3_ln",
                "ch_OXBE_total_i_cor_M2_ln",
                "ch_OXBE_total_i_cor_Y1_ln",
                
                "mo_TRCS_total_i_cor_t2_ln",
                "mo_TRCS_total_i_cor_t3_ln",
                "ch_TRCS_total_i_cor_M2_ln",                                     
                "ch_TRCS_total_i_cor_Y1_ln", 
                
                "mo_PFOA_cor_ln",
                "mo_PFNA_ln",
                "mo_PFHxS_cor_ln",
                "mo_PFOS_cor_ln",
                
                "mo_PFDA_i_cor_ln",
                "mo_PFUnDA_i_cor_ln", 
                "mo_PFHpS_i_cor_ln",
                
                "mo_PFDoDa_cat", 
                "mo_PFHxPA_cat", 
                "mo_PFHpA_cat_2",
                "mo_PFTrDa_cat", 
                "mo_PFBS_cat_2", 
                "mo_PFOSA_cat_2", 
                "mo_6_2diPAP_cat_2",
                "mo_8_2diPAP_cat_2")

## Metadata ----
load("2_final_data/metadata.RData")
input_metadata <- metadata %>%
  select(ident,                      
         all_of(covariates), 
         ch_MEPA_total_i_cor_Y1_ln, ch_BPA_total_i_cor_M2_ln, ch_BPA_total_i_cor_Y1_ln, # pour analyses de sensibilité
         ch_MEPA_free_i_cor_Y1_ln, ch_BPA_free_i_cor_M2_ln, ch_BPA_free_i_cor_Y1_ln,    # pour analyses de sensibilité
         all_of(pollutants)) %>%
  filter(!is.na(ch_feces_RUN_Y1)) %>%
  column_to_rownames("ident")
rm(metadata)

input_metadata <- input_metadata %>%
  mutate_at(vars(c("ch_feces_RUN_Y1",                   
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
                   "bf_duration_till48w_4cat_i")), as.factor)

input_metadata <- input_metadata %>%
  select(mo_MEPA_total_i_cor_t2_ln, 
         all_of(covariates))


## microbiota data ----
input_data <- 
  read_labelled_csv(
    "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
  select(
    ident, 
    starts_with("ch_feces_raw_g")) %>%
  filter(!is.na(ch_feces_raw_g1_Y1)) %>%
  column_to_rownames("ident")

var_label(input_data) <-  str_replace(
  var_label(input_data),                                    # set correct variable names                   
  "One year child feces relative abundance of ", "")
colnames(input_data) <- var_label(input_data)

input_data <- input_data %>%
  mutate_all(., ~ ifelse(. == 0, 1/5000, .))


test_maaslin <- 
  Maaslin2(input_data      = input_data,           # fichier d'abondance brutes des genres bactéries
           input_metadata  = input_metadata,       # fichier des varibales d'expo et des covariables 
           output          = "4_output/maaslin2",  # emplacement et nom du fichier de sortie
           normalization   = "TSS",                # normalisé en abondance relative
           transform       = "LOG",                # transformé en log
           analysis_method = "LM",                 # choisir la régression linéaire
           min_abundance   = 2/5000,               # the minimum abundance for each feature.
           min_prevalence  = 30,                   # The minimum % of samples for which a feature is detected at minimum abundance.
           
           fixed_effects   = "mo_MEPA_total_i_cor_t2_ln",  # explanatory variable 
           reference       = c("Between 0 and 6 months old,Under 2 years,2years or less after graduation,<3 Kg,<50 cm,<8.5 Kg,<75 cm,<19 Kg/m2,Not breastfed"))      # catégories de référence (pour les covariables catégorielles)


