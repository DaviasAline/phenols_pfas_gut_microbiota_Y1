## Aline Davias 
## 2021/10/08
## Metadonnées SEPAGES générés dans le cadre Projet GUMME 


# 0_data_reading ----
## Chargement des packages ----
library(haven)
library(tidyverse)
library(questionr)
library(lubridate)
library(mice)
library(labelled)
library(readr)
library(expss)
library(gtsummary)

## Chargement des données ----
source("3_programs/0_source_data_reading_AD_gumme.R", encoding = 'UTF-8')
age_feces_Y1 <- read_labelled_csv("0_source_data/age_feces_collection_Y1_labelled_AD_20220504_1.csv")
atb_Y1 <- read_labelled_csv("0_source_data/antibiotics_Y1_labelled_AD_20220314_9.csv")
solidfood_Y1 <- read_labelled_csv("0_source_data/solidfood_Y1_labelled_AD_20220302_7.csv")
weight_length_Y1 <- read_labelled_csv("0_source_data/weight_length_Y1_labelled_AD_20220301_2.csv")
metadata_microbiote <- 
  read_labelled_csv("0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") %>%
  select(!starts_with("ch_feces_rel")) %>%
  select(!starts_with("ch_feces_raw")) 
  

var_lab(metadata$ident) <- NULL
var_lab(metadata_microbiote$ident) <- NULL
var_lab(metadata_microbiote$ch_feces_ID_Y1) <- NULL
var_lab(age_feces_Y1$ident) <- NULL
var_lab(atb_Y1$ident) <- NULL
var_lab(solidfood_Y1$ident) <- NULL
var_lab(weight_length_Y1$ident) <- NULL

metadata <- left_join(metadata, 
                      metadata_microbiote, 
                      by = c("ident", "ch_feces_ID_Y1"))
metadata <- left_join(metadata, 
                      age_feces_Y1, 
                      by = "ident")
metadata <- left_join(metadata, 
                      atb_Y1, 
                      by = "ident")
metadata <- left_join(metadata, 
                      solidfood_Y1, 
                      by = "ident")
metadata <- left_join(metadata, 
                      weight_length_Y1, 
                      by = "ident")



# 1_data_cleaning ----
## Data cleaning ----
metadata <- metadata %>%
  mutate(
    statut = if_else(!is.na(ch_feces_OrderSeq_ASVbased_Y1), "inclu", "exclu"), 
    statut = as.factor(statut)
  ) %>%
  rename(
    mo_ethnicity = mt1saa1_q02,
    ch_sibling = mt1saa1_q04,
    fa_dipl = ft2sac1_q03,
    mo_pets = mt3eaf1_q01,       # Variables animaux (questionnaire MT3EAF1_V1)
    mo_dogs = mt3eaf1_q01p1,
    mo_cats = mt3eaf1_q01p5,
    mo_birds = mt3eaf1_q01p9,
    mo_rodents = mt3eaf1_q01p13,
    mo_other_pets =  mt3eaf1_q01p17
  )


covar_vec_cat <- c("ch_sex", "po_delmod","mo_par", "mo_dipl", "mo_ethnicity","ch_sibling","fa_dipl","mo_pets",
                   "mo_dogs","mo_cats","mo_birds","mo_rodents","mo_other_pets","mo_tob_gr_anyt_yn_n2","mo_tob_gr_anyt_yn1_n2",
                   "ch_ETS_12m_opt36m","Mo_ETS_anyT_yn1_opt")
metadata[, covar_vec_cat] <- lapply(metadata[, covar_vec_cat], as.character)
metadata[, covar_vec_cat] <- lapply(metadata[, covar_vec_cat], as.factor)

metadata <- metadata %>%
  mutate(
    mo_tob_gr_anyt_yn_n2 = fct_recode(mo_tob_gr_anyt_yn_n2,            # tabagisme actif de la mère pdt la grossesse
                                      "No" = "0",
                                      "Yes" = "1"), 
    mo_tob_gr_anyt_yn1_n2 = fct_recode(mo_tob_gr_anyt_yn1_n2,
                                       "No" = "0",
                                       "Yes" = "1"), 
    ch_ETS_12m_opt36m = fct_recode(ch_ETS_12m_opt36m,                  # tabagisme passif de l'enfant 
                                   "No" = "0",
                                   "Yes" = "1"), 
    Mo_ETS_anyT_yn1_opt = fct_recode(Mo_ETS_anyT_yn1_opt,              # tabagisme passif de la mère pdt la grossesse
                                     "No" = "0",
                                     "Yes" = "1"), 
    ch_sex = fct_recode(ch_sex,                                        # Labellisation de ch_sex
                        "Male" = "1", 
                        "Female" = "2"), 
    po_delmod = fct_recode(po_delmod,                                  # Labellisation de po_delmod
                           "Vaginal delivery" = "1", 
                           "C-section" = "2"), 
    mo_par_2cat = fct_recode(mo_par,                                        # Labellisation de mo_par
                             "None" = "0", 
                             "1 child or more" = "1", 
                             "1 child or more" = "2"), 
    ch_sibling_3cat = fct_recode(ch_sibling,                                # Labellisation de ch_sibling
                                 "None" = "0",
                                 "1 child" = "1",
                                 "2 child or more" = "2",
                                 "2 child or more" = "3"), 
    mo_ethnicity = fct_recode(mo_ethnicity,                            # Labellisation de mo_ethnicity
                              "Afrique" = "1",
                              "Amériques" = "2",
                              "Asie du Sud-Est" = "3",
                              "Europe" = "4",
                              "Méditéranée Orientale" = "5",
                              "Pacifique Occidental" = "6",            
                              "Autre" = "7", 
                              "Ne sait pas/ne souhaite pas répondre" = "99"), 
    mo_dipl = fct_recode(mo_dipl,                                      # Labellisation des variables éducation mère mo_dipl
                         "BEP/CAP/Highschool" = "2", 
                         "1-2years after graduation" = "3", 
                         "3-4years after graduation" = "4", 
                         ">=5years after graduation" = "5"), 
    fa_dipl = fct_recode(fa_dipl,                                      # Labellisation des variables éducation père fa_dipl              
                         "BEP/CAP/Highschool" = "1", 
                         "BEP/CAP/Highschool" = "2", 
                         "BEP/CAP/Highschool" = "3", 
                         "BEP/CAP/Highschool" = "4", 
                         "BEP/CAP/Highschool" = "5", 
                         "1-2years after graduation" = "6",
                         "3-4years after graduation" = "7",
                         ">=5years after graduation" = "8"), 
    mo_pets = fct_recode(mo_pets,                                      # Labellisation des variables animaux 
                         "No" = "0", 
                         "One or more" = "1"), 
    mo_dogs = fct_recode(mo_dogs, 
                         "No" = "0", 
                         "One or more" = "1"), 
    mo_cats = fct_recode(mo_cats, 
                         "No" = "0", 
                         "One or more" = "1"), 
    mo_birds = fct_recode(mo_birds, 
                          "No" = "0", 
                          "One or more" = "1"), 
    mo_rodents = fct_recode(mo_rodents, 
                            "No" = "0", 
                            "One or more" = "1"), 
    mo_other_pets = fct_recode(mo_other_pets, 
                               "No" = "0", 
                               "One or more" = "1")
  )



# warning message mo_ethnicity : Unknown levels in `f`: 6  : normal car absence d'observation pour cette modalité 
# warning message fa_dipl : Unknown levels in `f`: 1, 2  : normal car absence d'observation pour ces modalités 



## Création variable intervalle intergrossesse ----
metadata <- metadata %>% 
  rename(
    mo_number_previous_pregnancy = mt1haa1_q22,                   # Nombre de grossesse antérieures 
    mo_first_pregnancy = mt1haa1_q22p1,                           # Infos intervalle avec 1ère grossesse
    mo_month_first_pregnancy = mt1haa1_q22p3, 
    mo_year_first_pregnancy = mt1haa1_q22p4,        
    mo_second_pregnancy = mt1haa1_q22p19,                        # Infos intervalle avec 2nd grossesse
    mo_month_second_pregnancy = mt1haa1_q22p21, 
    mo_year_second_pregnancy = mt1haa1_q22p22, 
    mo_third_pregnancy = mt1haa1_q22p37,                         # Infos intervalle avec 3eme grossesse
    mo_month_third_pregnancy = mt1haa1_q22p39, 
    mo_year_third_pregnancy = mt1haa1_q22p40,
    mo_fourth_pregnancy = mt1haa1_q22p55,                        # Infos intervalle avec 4eme grossesse
    mo_month_fourth_pregnancy = mt1haa1_q22p57, 
    mo_year_fourth_pregnancy = mt1haa1_q22p58,
    mo_fifth_pregnancy = mt1haa1_q22p73,                         # Infos intervalle avec 5eme grossesse
    mo_month_fifth_pregnancy = mt1haa1_q22p75, 
    mo_year_fifth_pregnancy = mt1haa1_q22p76,
    mo_sixth_pregnancy = mt1haa1_q22p91,                         # Infos intervalle avec 6eme grossesse
    mo_month_sixth_pregnancy = mt1haa1_q22p93, 
    mo_year_sixth_pregnancy = mt1haa1_q22p94) %>%
  mutate(
    mo_number_previous_pregnancy = as.double(mo_number_previous_pregnancy), 
    mo_first_pregnancy = as.double(mo_first_pregnancy), 
    mo_month_first_pregnancy = as.double(mo_month_first_pregnancy), 
    mo_year_first_pregnancy = as.double(mo_year_first_pregnancy), 
    mo_second_pregnancy = as.double(mo_second_pregnancy), 
    mo_month_second_pregnancy = as.double(mo_month_second_pregnancy), 
    mo_year_second_pregnancy = as.double(mo_year_second_pregnancy),
    mo_third_pregnancy = as.double(mo_third_pregnancy), 
    mo_month_third_pregnancy = as.double(mo_month_third_pregnancy), 
    mo_year_third_pregnancy = as.double(mo_year_third_pregnancy), 
    mo_fourth_pregnancy = as.double(mo_fourth_pregnancy), 
    mo_month_fourth_pregnancy = as.double(mo_month_fourth_pregnancy), 
    mo_year_fourth_pregnancy = as.double(mo_year_fourth_pregnancy), 
    mo_fifth_pregnancy = as.double(mo_fifth_pregnancy), 
    mo_month_fifth_pregnancy = as.double(mo_month_fifth_pregnancy), 
    mo_year_fifth_pregnancy = as.double(mo_year_fifth_pregnancy), 
    mo_sixth_pregnancy = as.double(mo_sixth_pregnancy), 
    mo_month_sixth_pregnancy = as.double(mo_month_sixth_pregnancy), 
    mo_year_sixth_pregnancy = as.double(mo_year_sixth_pregnancy), 
  )


## Pour chaque variable mois codé en NA ou en 99 alors que l'on a une valeur pour l'année, on impute par la valeur 1
metadata <- metadata %>% 
  mutate(
    mo_month_first_pregnancy = if_else(mo_month_first_pregnancy < 13, mo_month_first_pregnancy, 1),
    mo_month_second_pregnancy = if_else(mo_month_second_pregnancy < 13, mo_month_second_pregnancy, 1),
    mo_month_third_pregnancy = if_else(mo_month_third_pregnancy < 13, mo_month_third_pregnancy, 1),
    mo_month_fourth_pregnancy = if_else(mo_month_fourth_pregnancy < 13, mo_month_fourth_pregnancy, 1),
    mo_month_fifth_pregnancy = if_else(mo_month_fifth_pregnancy < 13, mo_month_fifth_pregnancy, 1),
    mo_month_sixth_pregnancy = if_else(mo_month_sixth_pregnancy < 13, mo_month_sixth_pregnancy, 1)
  ) %>% 
  mutate(                                                       
    mo_month_first_pregnancy = if_else((!is.na(mo_year_first_pregnancy)) & (is.na(mo_month_first_pregnancy)),
                                       1, 
                                       mo_month_first_pregnancy), 
    mo_month_second_pregnancy = if_else((!is.na(mo_year_second_pregnancy)) & (is.na(mo_month_second_pregnancy)),
                                        1, 
                                        mo_month_second_pregnancy), 
    mo_month_third_pregnancy = if_else((!is.na(mo_year_third_pregnancy)) & (is.na(mo_month_third_pregnancy)),
                                       1, 
                                       mo_month_third_pregnancy), 
    mo_month_fourth_pregnancy = if_else((!is.na(mo_year_fourth_pregnancy)) & (is.na(mo_month_fourth_pregnancy)),
                                        1, 
                                        mo_month_fourth_pregnancy), 
    mo_month_fifth_pregnancy = if_else((!is.na(mo_year_fifth_pregnancy)) & (is.na(mo_month_fifth_pregnancy)),
                                       1, 
                                       mo_month_fifth_pregnancy), 
    mo_month_sixth_pregnancy = if_else((!is.na(mo_year_sixth_pregnancy)) & (is.na(mo_month_sixth_pregnancy)),
                                       1, mo_month_sixth_pregnancy)) %>%
  
  mutate(                               # Rassembler le mois/l'année pour chaque grossesse
    mo_date_first_pregnancy = str_c(mo_year_first_pregnancy, mo_month_first_pregnancy, sep="-"), 
    mo_date_second_pregnancy = str_c(mo_year_second_pregnancy, mo_month_second_pregnancy, sep="-"),
    mo_date_third_pregnancy = str_c(mo_year_third_pregnancy, mo_month_third_pregnancy, sep="-"), 
    mo_date_fourth_pregnancy = str_c(mo_year_fourth_pregnancy, mo_month_fourth_pregnancy, sep="-"), 
    mo_date_fifth_pregnancy = str_c(mo_year_fifth_pregnancy, mo_month_fifth_pregnancy, sep="-"), 
    mo_date_sixth_pregnancy = str_c(mo_year_sixth_pregnancy, mo_month_sixth_pregnancy, sep="-")) %>% 
  
  mutate(                                                      # Conversion en variable date 
    mo_date_first_pregnancy = ym(mo_date_first_pregnancy), 
    mo_date_second_pregnancy = ym(mo_date_second_pregnancy), 
    mo_date_third_pregnancy = ym(mo_date_third_pregnancy), 
    mo_date_fourth_pregnancy = ym(mo_date_fourth_pregnancy), 
    mo_date_fifth_pregnancy = ym(mo_date_fifth_pregnancy), 
    mo_date_sixth_pregnancy = ym(mo_date_sixth_pregnancy))




## Création de la variable date de la dernière grossesse ayant conduit à un enfant vivant 

metadata <- metadata %>%
  mutate(
    mo_date_last_pregnancy = NA, 
    mo_date_last_pregnancy = as.numeric(mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_sixth_pregnancy == 1, 
                                          mo_date_sixth_pregnancy, 
                                          mo_date_last_pregnancy))%>%
  mutate(mo_date_last_pregnancy = if_else(mo_fifth_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_fifth_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_fourth_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_fourth_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_third_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_third_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_second_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_second_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_first_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_first_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(
    mo_number_previous_pregnancy = as.factor(mo_number_previous_pregnancy),
    mo_interpreg = difftime(po_datedel, 
                            mo_date_last_pregnancy, 
                            units = "weeks"),
    mo_interpreg = as.numeric(mo_interpreg), 
    mo_interpreg = mo_interpreg / 52.1429,             # Transformation en année
    mo_interpreg_5cat = cut(mo_interpreg,              # création variable catégorielle 
                            include.lowest = TRUE,
                            right = FALSE,
                            dig.lab = 4,
                            breaks = c(1, 2, 3, 4, 11.9945106894433)), 
    mo_interpreg_5cat = fct_recode(mo_interpreg_5cat,
                                   "Under 2 years" = "[1,2)",
                                   "Between 2 and 3 years" = "[2,3)",
                                   "Between 3 and 4 years" = "[3,4)",
                                   "Over 4 years" = "[4,11.99]"), 
    mo_interpreg_5cat = as.character(mo_interpreg_5cat), 
    mo_par_2cat = as.character(mo_par_2cat), 
    mo_interpreg_5cat = if_else(mo_number_previous_pregnancy == 0 & (is.na(mo_interpreg_5cat)),   # Prise en compte des NA qui sont en fait des femmes primipares
                                "Primiparous", 
                                mo_interpreg_5cat), 
    mo_interpreg_5cat = if_else( mo_par_2cat == "None" & (is.na(mo_interpreg_5cat)), 
                                 "Primiparous", 
                                 mo_interpreg_5cat), 
    mo_interpreg_5cat = as.factor(mo_interpreg_5cat), 
    mo_par_2cat = as.factor( mo_par_2cat)) %>%
  mutate(
    mo_interpreg_5cat = fct_relevel(mo_interpreg_5cat,
                                    "Under 2 years", "Between 2 and 3 years", "Between 3 and 4 years", "Over 4 years", "Primiparous"),
    mo_interpreg_3cat = fct_recode(mo_interpreg_5cat,
                                   "2 years and more" = "Between 2 and 3 years",
                                   "2 years and more" = "Between 3 and 4 years",
                                   "2 years and more" = "Over 4 years"))





## Choix du codage des covariables ----
metadata <- metadata %>%                                    # création de variables catégorielles à partir de variables numériques 
  mutate(
    ch_feces_age_w_Y1 = as.numeric(ch_feces_age_w_Y1),
    ch_feces_age_w_Y1_4cat = cut(
      ch_feces_age_w_Y1,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(42, 51, 53, 55, 71)
    ),
    po_w_kg = po_w / 1000,
    po_w_kg_3cat = cut(po_w_kg,
                       include.lowest = TRUE,
                       right = FALSE,
                       dig.lab = 4,
                       breaks = c(0.9, 3, 3.5, 4.7)), 
    po_he_3cat = cut(po_he,
                     include.lowest = TRUE,
                     right = FALSE,
                     dig.lab = 4,
                     breaks = c(30, 50, 52, 60)),
    
    ch_w_Y1_3cat = cut(ch_w_Y1,
                       include.lowest = TRUE,
                       right = FALSE,
                       dig.lab = 4,
                       breaks = c(6.0, 8.5, 10, 14)),
    ch_he_Y1_3cat = cut(ch_he_Y1,
                        include.lowest = TRUE,
                        right = FALSE,
                        dig.lab = 4,
                        breaks = c(67.75, 75, 78, 83.75)),
    mo_bmi_bepr_3cat = cut(mo_bmi_bepr,
                           include.lowest = TRUE, 
                           right = FALSE,
                           dig.lab = 4,
                           breaks = c(16, 19, 24, 42)),
    po_gd_4cat = cut(po_gd,                                                 
                     include.lowest = TRUE,
                     right = FALSE,
                     dig.lab = 4,
                     breaks = c(28, 38, 40, 41, 42)), 
    mo_age_4cat = cut(mo_age,                                                
                      include.lowest = TRUE,
                      right = FALSE,
                      dig.lab = 4,
                      breaks = c(20, 27, 33, 36, 46)), 
    bf_duration_till48w_4cat = cut(bf_duration_till48w,                         
                                   include.lowest = TRUE,
                                   right = FALSE,
                                   dig.lab = 4,
                                   breaks = c(0, 1, 24, 47, 48))) %>%
  
  
  mutate(                                  # création d'étiquettes de catégories pour les nouvelles variables 
    ch_feces_age_w_Y1_4cat = fct_recode(
      ch_feces_age_w_Y1_4cat,
      "<51weeks" = "[42,51)",
      "51-52weeks" = "[51,53)",
      "53-54weeks" = "[53,55)",
      ">=55weeks" = "[55,71]"
    ),
    ch_antibio_Y1 = as.character(ch_antibio_Y1),
    ch_antibio_Y1_3cat = fct_recode(
      ch_antibio_Y1,
      "2 and more" = "2",
      "2 and more" = "3",
      "2 and more" = "4",
      "2 and more" = "5"
    ),
    ch_antibio_Y1_2cat = fct_recode(
      ch_antibio_Y1,
      "No" = "0",
      "Yes" = "1",
      "Yes" = "2",
      "Yes" = "3",
      "Yes" = "4",
      "Yes" = "5"
    ),
    ch_food_intro_Y1_3cat = fct_recode(ch_food_intro_Y1,
                                       "Between 0 and 6 months old" = "Between 3 and 6 months old",
                                       "Between 0 and 6 months old" = "Between 0 and 3 months old"
    ), 
    po_w_kg_3cat = fct_recode(po_w_kg_3cat,
                              "<3 Kg" = "[0.9,3)",
                              "3-3.4 Kg" = "[3,3.5)",
                              ">= 3.5 Kg" = "[3.5,4.7]"),
    po_he_3cat = fct_recode(po_he_3cat,
                            "<50 cm" = "[30,50)",
                            "50-51 cm" = "[50,52)",
                            ">= 52 cm" = "[52,60]"),
    ch_w_Y1_3cat = fct_recode(ch_w_Y1_3cat,
                              "<8.5 Kg" = "[6,8.5)",
                              "8.5-9.9 Kg" = "[8.5,10)",
                              ">=10 Kg" = "[10,14]"),
    ch_he_Y1_3cat = fct_recode(ch_he_Y1_3cat,
                               "<75 cm" = "[67.75,75)",
                               "75-77.9 cm" = "[75,78)",
                               ">=78 cm" = "[78,83.75]"),
    mo_bmi_bepr_3cat = fct_recode(mo_bmi_bepr_3cat,
                                  "<19 Kg/m2" = "[16,19)",
                                  "19-23.9 Kg/m2" = "[19,24)",
                                  ">=24 Kg/m2" = "[24,42]"),
    po_gd_4cat = fct_recode(po_gd_4cat,                                     
                            "<38 weeks of amenorrhea" = "[28,38)",
                            "38-39 weeks of amenorrhea" = "[38,40)",
                            "40 weeks of amenorrhea" = "[40,41)",
                            "> 40 weeks of amenorrhea" = "[41,42]"), 
    mo_age_4cat = fct_recode(mo_age_4cat,                                      
                             "20-26 years" = "[20,27)",
                             "27-32 years" = "[27,33)",
                             "33-35 years" = "[33,36)",
                             ">35 years" = "[36,46]"), 
    bf_duration_till48w_4cat = fct_recode(bf_duration_till48w_4cat,         
                                          "Not breastfed" = "[0,1)",
                                          "<24 weeks" = "[1,24)",
                                          "24-47 weeks" = "[24,47)",
                                          "Still breastfeed at 48 weeks" = "[47,48]"),                                       
    mo_interpreg_3cat = fct_relevel(mo_interpreg_3cat,
                                    "Primiparous", 
                                    "2 years and more", 
                                    "Under 2 years"), 
    mo_dipl = as.character(mo_dipl),
    mo_dipl_3cat = fct_recode(mo_dipl,
                              "2years or less after graduation" = "BEP/CAP/Highschool",
                              "2years or less after graduation" = "1-2years after graduation"))


## Imputation des valeurs manquantes covariables ----
metadata %>%
  filter(statut == "inclu") %>%
  select(                 
    ch_feces_age_w_Y1,     # 7 NA /356
    po_gd,                 # 0 NA /356
    mo_age,                # 0 NA /356
    mo_bmi_bepr,           # 1 NA /356
    mo_par,                # 0 NA /356
    po_w_kg,               # 0 NA /356
    po_he,                 # 1 NA /356
    ch_w_Y1,               # 5 NA /356
    ch_he_Y1,              # 13 NA /356
    bf_duration_till48w,   # 9 NA /356
    ch_antibio_Y1,         # 1 NA /356
    po_delmod,             # 0 NA /356
    ch_sex,                # 0 NA /356
    mo_dipl_3cat,          # 2 NA /356
    mo_pets,               # 45 NA /356
    mo_interpreg_3cat,     # 0 NA /356
    mo_tob_gr_anyt_yn_n2,  # 26 NA /356
    Mo_ETS_anyT_yn1_opt,   # 17 NA /356
    ch_ETS_12m_opt36m,     # 0 NA /356
    ch_food_intro_Y1       # 37 NA /356
  ) %>%
  tbl_summary()



### Imputations (package mice)
bdd_imput <- metadata %>%
  mutate(
    ch_antibio_Y1 = as.numeric(ch_antibio_Y1), 
    ch_food_intro_Y1 = as.factor(ch_food_intro_Y1)) %>%
  select(ident, 
         ch_feces_age_w_Y1,  # variables à imputer en continu
         mo_bmi_bepr,
         po_he, 
         ch_w_Y1, 
         ch_he_Y1,
         bf_duration_till48w,
         ch_antibio_Y1,
         
         mo_dipl_3cat,      # variables à imuter en continue
         mo_pets,
         mo_tob_gr_anyt_yn_n2,
         Mo_ETS_anyT_yn1_opt, 
         ch_food_intro_Y1)

str(bdd_imput)               # vérifier le codage des variables 

imput <-  mice(bdd_imput, m=1, maxit=50, seed=500)   # imputation

imput$method    # voir quelles méthodes ont été utilisés pour chaque variable à imputer

imput$imp$ch_feces_age_w_Y1     # voir les imputations 
imput$imp$mo_bmi_bepr
imput$imp$po_he 
imput$imp$ch_w_Y1
imput$imp$ch_he_Y1
imput$imp$bf_duration_till48w
imput$imp$ch_antibio_Y1
imput$imp$mo_dipl_3cat  
imput$imp$mo_pets
imput$imp$mo_tob_gr_anyt_yn_n2
imput$imp$Mo_ETS_anyT_yn1_opt 
imput$imp$ch_food_intro_Y1

bdd_imput <- complete(imput)
bdd_imput <- as.data.frame(bdd_imput)
colnames(bdd_imput) <- c("ident", paste(colnames(bdd_imput[, 2:13]), "i", sep="_"))  # les nouvelles variables imputées se terminent en _i

metadata <- left_join(metadata,    # merger les variables imputées au reste de la base de données
                      bdd_imput,
                      by = "ident")  

# on suprime les valeurs imputées créées chez des enfants qui n'ont pas eu de prélévement de selles à 1 an
metadata <- metadata %>%
  mutate(ch_feces_age_w_Y1_i = 
           if_else(statut == "inclu", ch_feces_age_w_Y1_i, ch_feces_age_w_Y1))   
# parmis les "exclu" de l'analyse, il y a quand même 4 valeurs pour la variable age au prélévelement de la selle 
# --> correspond à des échantillons avec une faible quantité d'adn donc exclu des analyses


metadata <- metadata %>%
  mutate(
    ch_feces_age_w_Y1_4cat_i = cut(
      ch_feces_age_w_Y1_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(42, 51, 53, 55, 71),
    ),
    mo_bmi_bepr_3cat_i = cut(
      mo_bmi_bepr_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(16, 19, 24, 42)
    ),
    po_he_3cat_i = cut(
      po_he_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(30, 50, 52, 60)
    ),
    ch_w_Y1_3cat_i = cut(
      ch_w_Y1_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(6.0, 8.5, 10, 14)
    ),
    ch_he_Y1_3cat_i = cut(
      ch_he_Y1_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(67.75, 75, 78, 83.75)
    ),
    bf_duration_till48w_4cat_i = cut(
      bf_duration_till48w_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(0, 1, 24, 47, 48)
    )
  ) %>%
  
  mutate(                                     # création d'étiquettes de catégories pour les nouvelles variables
    ch_feces_age_w_Y1_4cat_i = fct_recode(
      ch_feces_age_w_Y1_4cat_i,
      "<51weeks" = "[42,51)",
      "51-52weeks" = "[51,53)",
      "53-54weeks" = "[53,55)",
      ">=55weeks" = "[55,71]"
    ),
    mo_bmi_bepr_3cat_i = fct_recode(
      mo_bmi_bepr_3cat_i,
      "<19 Kg/m2" = "[16,19)",
      "19-23.9 Kg/m2" = "[19,24)",
      ">=24 Kg/m2" = "[24,42]"
    ),
    po_he_3cat_i = fct_recode(
      po_he_3cat_i,
      "<50 cm" = "[30,50)",
      "50-51 cm" = "[50,52)",
      ">= 52 cm" = "[52,60]"
    ),
    ch_w_Y1_3cat_i = fct_recode(
      ch_w_Y1_3cat_i,
      "<8.5 Kg" = "[6,8.5)",
      "8.5-9.9 Kg" = "[8.5,10)",
      ">=10 Kg" = "[10,14]"
    ),
    ch_he_Y1_3cat_i = fct_recode(
      ch_he_Y1_3cat_i,
      "<75 cm" = "[67.75,75)",
      "75-77.9 cm" = "[75,78)",
      ">=78 cm" = "[78,83.75]"
    ),
    bf_duration_till48w_4cat_i = fct_recode(
      bf_duration_till48w_4cat_i,
      "Not breastfed" = "[0,1)",
      "<24 weeks" = "[1,24)",
      "24-47 weeks" = "[24,47)",
      "Still breastfeed at 48 weeks" = "[47,48]"
    ),
    ch_antibio_Y1_3cat_i = fct_recode(
      as.factor(as.character(ch_antibio_Y1_i)),
      "2 and more" = "2",
      "2 and more" = "3",
      "2 and more" = "4",
      "2 and more" = "5"
    ),
    ch_antibio_Y1_2cat_i = fct_recode(
      as.factor(as.character(ch_antibio_Y1_i)),
      "No" = "0",
      "Yes" = "1",
      "Yes" = "2",
      "Yes" = "3",
      "Yes" = "4",
      "Yes" = "5"
    ),
    ch_food_intro_Y1_3cat_i = fct_recode(
      ch_food_intro_Y1_i,
      "Between 0 and 6 months old" = "Between 3 and 6 months old",
      "Between 0 and 6 months old" = "Between 0 and 3 months old"
    )
  )


metadata <- metadata %>% mutate(
  ch_food_intro_Y1_3cat_i = fct_relevel(
    ch_food_intro_Y1_3cat_i,
    "Between 0 and 6 months old", "Between 6 and 12 months old", "Not introduced at 12 months old"),
  mo_par_2cat = fct_relevel(
    mo_par_2cat,
    "None", "1 child or more"),
  mo_interpreg_3cat = fct_relevel(
    mo_interpreg_3cat,
    "Under 2 years", "2 years and more", "Primiparous"),
  mo_dipl_3cat_i = fct_relevel(
    mo_dipl_3cat_i,
    "2years or less after graduation", "3-4years after graduation", ">=5years after graduation"),
  
  po_w_kg_3cat = fct_relevel(
    po_w_kg_3cat,
    "<3 Kg", "3-3.4 Kg", ">= 3.5 Kg"),
  po_he_3cat_i = fct_relevel(
    po_he_3cat_i,
    "<50 cm", "50-51 cm", ">= 52 cm"),
  ch_w_Y1_3cat_i = fct_relevel(
    ch_w_Y1_3cat_i,
    "<8.5 Kg", "8.5-9.9 Kg", ">=10 Kg"),
  ch_he_Y1_3cat_i = fct_relevel(
    ch_he_Y1_3cat_i,
    "<75 cm", "75-77.9 cm", ">=78 cm"),
  mo_bmi_bepr_3cat_i = fct_relevel(
    mo_bmi_bepr_3cat_i,
    "<19 Kg/m2", "19-23.9 Kg/m2", ">=24 Kg/m2"),
  
  bf_duration_till48w_4cat_i = fct_relevel(
    bf_duration_till48w_4cat_i,
    "Not breastfed", "<24 weeks", "24-47 weeks", "Still breastfeed at 48 weeks"))

