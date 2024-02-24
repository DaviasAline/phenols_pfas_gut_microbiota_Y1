# Analyses  exposition aux polluants et microbiote 1 an 
# A. Davias
# 23/02/2023
# Tester un modèle de régression linéaire par fenetre d'exposition

# Chargement des données ----
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')       # fonctions
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)
library(lmtest)


# Préparation des données ----
bdd_alpha_t2 <- bdd_alpha %>% select(ident, 
                                     ch_feces_SpecRich_5000_ASV_Y1, 
                                     ch_feces_Shannon_5000_ASV_Y1, 
                                     ch_feces_Faith_5000_ASV_Y1,
                                     all_of(covar_vec_i),
                                     all_of(pollutant_vec_t2)) %>% na.omit()
bdd_alpha_t3 <- bdd_alpha %>% select(ident, 
                                     ch_feces_SpecRich_5000_ASV_Y1, 
                                     ch_feces_Shannon_5000_ASV_Y1, 
                                     ch_feces_Faith_5000_ASV_Y1,
                                     all_of(covar_vec_i),
                                     all_of(pollutant_vec_t3)) %>% na.omit()
bdd_alpha_M2 <- bdd_alpha %>% select(ident, 
                                     ch_feces_SpecRich_5000_ASV_Y1, 
                                     ch_feces_Shannon_5000_ASV_Y1, 
                                     ch_feces_Faith_5000_ASV_Y1,
                                     all_of(covar_vec_i),
                                     all_of(pollutant_vec_M2)) %>% na.omit()
bdd_alpha_Y1 <- bdd_alpha %>% select(ident, 
                                     ch_feces_SpecRich_5000_ASV_Y1, 
                                     ch_feces_Shannon_5000_ASV_Y1, 
                                     ch_feces_Faith_5000_ASV_Y1,
                                     all_of(covar_vec_i),
                                     all_of(pollutant_vec_Y1)) %>% na.omit()




# Fonctions ----
model_covar <- function(outcome, data) {
  lm(
    outcome ~
      
      ch_feces_RUN_Y1  +
      ch_feces_age_w_Y1_i +
      po_delmod +
      ch_food_intro_Y1_3cat_i +
      ch_antibio_Y1_2cat_i +
      mo_par_2cat +
      mo_pets_i +
      ch_sex +
      mo_tob_gr_anyt_yn_n2_i +
      Mo_ETS_anyT_yn1_opt_i +
      ch_ETS_12m_opt36m +
      mo_interpreg_3cat +
      mo_dipl_3cat_i +
      po_w_kg_3cat +
      po_he_3cat_i +
      ch_w_Y1_3cat_i +
      ch_he_Y1_3cat_i +
      po_gd +
      mo_age +
      mo_bmi_bepr_3cat_i +
      bf_duration_till48w_4cat_i,
    data = data
  ) 
}


mixture_t2_function <- function(outcome, data) {
  lm(
    outcome ~
      mo_MEPA_total_i_cor_t2_ln +
      mo_ETPA_total_i_cor_t2_ln +
      mo_PRPA_total_i_cor_t2_ln +
      mo_BUPA_total_cat_t2 +
      mo_BPA_total_i_cor_t2_ln +
      mo_BPS_total_cat_t2_2 +
      mo_OXBE_total_i_cor_t2_ln +
      mo_TRCS_total_i_cor_t2_ln +
      mo_PFOA_cor_ln +
      mo_PFNA_ln +
      mo_PFHxS_cor_ln +
      mo_PFOS_cor_ln +
      mo_PFDA_i_cor_ln +
      mo_PFUnDA_i_cor_ln +
      mo_PFHpS_i_cor_ln +
      mo_PFDoDa_cat +
      mo_PFHxPA_cat +
      mo_PFHpA_cat_2 +
      mo_PFTrDa_cat +
      mo_PFBS_cat_2 +
      mo_PFOSA_cat_2 +
      mo_6_2diPAP_cat_2 +
      mo_8_2diPAP_cat_2 +
      
      ch_feces_RUN_Y1  +
      ch_feces_age_w_Y1_i +
      po_delmod +
      ch_food_intro_Y1_3cat_i +
      ch_antibio_Y1_2cat_i +
      mo_par_2cat +
      mo_pets_i +
      ch_sex +
      mo_tob_gr_anyt_yn_n2_i +
      Mo_ETS_anyT_yn1_opt_i +
      ch_ETS_12m_opt36m +
      mo_interpreg_3cat +
      mo_dipl_3cat_i +
      po_w_kg_3cat +
      po_he_3cat_i +
      ch_w_Y1_3cat_i +
      ch_he_Y1_3cat_i +
      po_gd +
      mo_age +
      mo_bmi_bepr_3cat_i +
      bf_duration_till48w_4cat_i,
    data = data
  ) 
}


mixture_t3_function <- function(outcome, data) {
  lm(outcome ~
       mo_MEPA_total_i_cor_t3_ln +
       mo_ETPA_total_i_cor_t3_ln +
       mo_PRPA_total_i_cor_t3_ln +
       mo_BUPA_total_cat_t3 +
       mo_BPA_total_i_cor_t3_ln +
       mo_BPS_total_cat_t3_2 +
       mo_OXBE_total_i_cor_t3_ln +
       mo_TRCS_total_i_cor_t3_ln +
       
       ch_feces_RUN_Y1  +
       ch_feces_age_w_Y1_i +
       po_delmod +
       ch_food_intro_Y1_3cat_i +
       ch_antibio_Y1_2cat_i +
       mo_par_2cat +
       mo_pets_i +
       ch_sex +
       mo_tob_gr_anyt_yn_n2_i +
       Mo_ETS_anyT_yn1_opt_i +
       ch_ETS_12m_opt36m +
       mo_interpreg_3cat +
       mo_dipl_3cat_i +
       po_w_kg_3cat +
       po_he_3cat_i +
       ch_w_Y1_3cat_i +
       ch_he_Y1_3cat_i +
       po_gd +
       mo_age +
       mo_bmi_bepr_3cat_i +
       bf_duration_till48w_4cat_i,
     data = data) 
}

mixture_M2_function <- function(outcome, data) {
  lm(outcome ~
       ch_MEPA_total_i_cor_M2_ln +
       ch_ETPA_total_cat_M2_2 +
       ch_PRPA_total_cat_M2_2 +
       ch_BUPA_total_cat_M2_2 +
       ch_BPA_conj_i_cor_M2_ln +
       ch_BPS_total_cat_M2_2 +
       ch_OXBE_total_i_cor_M2_ln +
       ch_TRCS_total_i_cor_M2_ln +
       
       ch_feces_RUN_Y1  +
       ch_feces_age_w_Y1_i +
       po_delmod +
       ch_food_intro_Y1_3cat_i +
       ch_antibio_Y1_2cat_i +
       mo_par_2cat +
       mo_pets_i +
       ch_sex +
       mo_tob_gr_anyt_yn_n2_i +
       Mo_ETS_anyT_yn1_opt_i +
       ch_ETS_12m_opt36m +
       mo_interpreg_3cat +
       mo_dipl_3cat_i +
       po_w_kg_3cat +
       po_he_3cat_i +
       ch_w_Y1_3cat_i +
       ch_he_Y1_3cat_i +
       po_gd +
       mo_age +
       mo_bmi_bepr_3cat_i +
       bf_duration_till48w_4cat_i,
     data = data) 
}


mixture_Y1_function <- function(outcome, data) {
  lm(outcome ~
       ch_MEPA_conj_i_cor_Y1_ln +
       ch_ETPA_total_i_cor_Y1_ln +
       ch_PRPA_total_i_cor_Y1_ln +
       ch_BUPA_total_cat_Y1 +
       ch_BPA_conj_i_cor_Y1_ln +
       ch_BPS_total_cat_Y1 +
       ch_OXBE_total_i_cor_Y1_ln +
       ch_TRCS_total_i_cor_Y1_ln +
       
       ch_feces_RUN_Y1  +
       ch_feces_age_w_Y1_i +
       po_delmod +
       ch_food_intro_Y1_3cat_i +
       ch_antibio_Y1_2cat_i +
       mo_par_2cat +
       mo_pets_i +
       ch_sex +
       mo_tob_gr_anyt_yn_n2_i +
       Mo_ETS_anyT_yn1_opt_i +
       ch_ETS_12m_opt36m +
       mo_interpreg_3cat +
       mo_dipl_3cat_i +
       po_w_kg_3cat +
       po_he_3cat_i +
       ch_w_Y1_3cat_i +
       ch_he_Y1_3cat_i +
       po_gd +
       mo_age +
       mo_bmi_bepr_3cat_i +
       bf_duration_till48w_4cat_i,
     data = data) 
}

model <- function(data, 
                  outcome, 
                  exposure_vec, 
                  digit_beta_IC) {
  data %>%                  
    select(
      {{outcome}},
      all_of({{exposure_vec}}),
      all_of(covar_vec_i)) %>% 
    
    tbl_uvregression(
      method = lm ,
      y = {{outcome}},
      formula = "{y} ~ 
         {x} +
         ch_feces_RUN_Y1  +
         ch_feces_age_w_Y1_i +
         po_delmod +
         ch_food_intro_Y1_3cat_i +
         ch_antibio_Y1_2cat_i +
         mo_par_2cat +
         mo_pets_i +
         ch_sex +
         mo_tob_gr_anyt_yn_n2_i +
         Mo_ETS_anyT_yn1_opt_i +
         ch_ETS_12m_opt36m +
         mo_interpreg_3cat +
         mo_dipl_3cat_i +
         po_w_kg_3cat +
         po_he_3cat_i +
         ch_w_Y1_3cat_i +
         ch_he_Y1_3cat_i +
         po_gd +
         mo_age +
         mo_bmi_bepr_3cat_i +
         bf_duration_till48w_4cat_i",
      hide_n = TRUE,
      pvalue_fun = ~ style_pvalue(.x, digits = 2),
      estimate_fun = ~ style_sigfig(.x, digits = {{digit_beta_IC}})
    ) %>%
    add_global_p(keep = TRUE, singular.ok = TRUE) %>%
    bold_labels() }


# Résultats ----
## Alpha diversity resultats bruts ----
model_spec_covar_t2 <-
  model_covar(outcome = bdd_alpha_t2$ch_feces_SpecRich_5000_ASV_Y1,
                      data = bdd_alpha_t2)
model_sha_covar_t2 <-
  model_covar(outcome = bdd_alpha_t2$ch_feces_Shannon_5000_ASV_Y1,
                      data = bdd_alpha_t2)
model_fai_covar_t2 <-
  model_covar(outcome = bdd_alpha_t2$ch_feces_Faith_5000_ASV_Y1,
                      data = bdd_alpha_t2)

model_spec_covar_t3 <-
  model_covar(outcome = bdd_alpha_t3$ch_feces_SpecRich_5000_ASV_Y1,
              data = bdd_alpha_t3)
model_sha_covar_t3 <-
  model_covar(outcome = bdd_alpha_t3$ch_feces_Shannon_5000_ASV_Y1,
              data = bdd_alpha_t3)
model_fai_covar_t3 <-
  model_covar(outcome = bdd_alpha_t3$ch_feces_Faith_5000_ASV_Y1,
              data = bdd_alpha_t3)

model_spec_covar_M2 <-
  model_covar(outcome = bdd_alpha_M2$ch_feces_SpecRich_5000_ASV_Y1,
              data = bdd_alpha_M2)
model_sha_covar_M2 <-
  model_covar(outcome = bdd_alpha_M2$ch_feces_Shannon_5000_ASV_Y1,
              data = bdd_alpha_M2)
model_fai_covar_M2 <-
  model_covar(outcome = bdd_alpha_M2$ch_feces_Faith_5000_ASV_Y1,
              data = bdd_alpha_M2)

model_spec_covar_Y1 <-
  model_covar(outcome = bdd_alpha_Y1$ch_feces_SpecRich_5000_ASV_Y1,
              data = bdd_alpha_Y1)
model_sha_covar_Y1 <-
  model_covar(outcome = bdd_alpha_Y1$ch_feces_Shannon_5000_ASV_Y1,
              data = bdd_alpha_Y1)
model_fai_covar_Y1 <-
  model_covar(outcome = bdd_alpha_Y1$ch_feces_Faith_5000_ASV_Y1,
              data = bdd_alpha_Y1)

mixture_spec_t2 <-
  mixture_t2_function(outcome = bdd_alpha_t2$ch_feces_SpecRich_5000_ASV_Y1,
                      data = bdd_alpha_t2)
mixture_sha_t2 <-
  mixture_t2_function(outcome = bdd_alpha_t2$ch_feces_Shannon_5000_ASV_Y1,
                      data = bdd_alpha_t2)
mixture_fai_t2 <-
  mixture_t2_function(outcome = bdd_alpha_t2$ch_feces_Faith_5000_ASV_Y1,
                      data = bdd_alpha_t2)

mixture_spec_t3 <-
  mixture_t3_function(outcome = bdd_alpha_t3$ch_feces_SpecRich_5000_ASV_Y1,
                      data = bdd_alpha_t3)
mixture_sha_t3 <-
  mixture_t3_function(outcome = bdd_alpha_t3$ch_feces_Shannon_5000_ASV_Y1,
                      data = bdd_alpha_t3)
mixture_fai_t3 <-
  mixture_t3_function(outcome = bdd_alpha_t3$ch_feces_Faith_5000_ASV_Y1,
                      data = bdd_alpha_t3)

mixture_spec_M2 <-
  mixture_M2_function(outcome = bdd_alpha_M2$ch_feces_SpecRich_5000_ASV_Y1,
                      data = bdd_alpha_M2)
mixture_sha_M2 <-
  mixture_M2_function(outcome = bdd_alpha_M2$ch_feces_Shannon_5000_ASV_Y1,
                      data = bdd_alpha_M2)
mixture_fai_M2 <-
  mixture_M2_function(outcome = bdd_alpha_M2$ch_feces_Faith_5000_ASV_Y1,
                      data = bdd_alpha_M2)

mixture_spec_Y1 <-
  mixture_Y1_function(outcome = bdd_alpha_Y1$ch_feces_SpecRich_5000_ASV_Y1,
                      data = bdd_alpha_Y1)
mixture_sha_Y1 <-
  mixture_Y1_function(outcome = bdd_alpha_Y1$ch_feces_Shannon_5000_ASV_Y1,
                      data = bdd_alpha_Y1)
mixture_fai_Y1 <-
  mixture_Y1_function(outcome = bdd_alpha_Y1$ch_feces_Faith_5000_ASV_Y1,
                      data = bdd_alpha_Y1)

test_log_vraisemblance_spec_t2 <- lrtest(mixture_spec_t2, model_spec_covar_t2)
test_log_vraisemblance_sha_t2 <- lrtest(mixture_sha_t2, model_sha_covar_t2)
test_log_vraisemblance_fai_t2 <- lrtest(mixture_fai_t2, model_fai_covar_t2)

test_log_vraisemblance_spec_t3 <- lrtest(mixture_spec_t3, model_spec_covar_t3)
test_log_vraisemblance_sha_t3 <- lrtest(mixture_sha_t3, model_sha_covar_t3)
test_log_vraisemblance_fai_t3 <- lrtest(mixture_fai_t3, model_fai_covar_t3)

test_log_vraisemblance_spec_M2 <- lrtest(mixture_spec_M2, model_spec_covar_M2)
test_log_vraisemblance_sha_M2 <- lrtest(mixture_sha_M2, model_sha_covar_M2)
test_log_vraisemblance_fai_M2 <- lrtest(mixture_fai_M2, model_fai_covar_M2)

test_log_vraisemblance_spec_Y1 <- lrtest(mixture_spec_Y1, model_spec_covar_Y1)
test_log_vraisemblance_sha_Y1 <- lrtest(mixture_sha_Y1, model_sha_covar_Y1)
test_log_vraisemblance_fai_Y1 <- lrtest(mixture_fai_Y1, model_fai_covar_Y1)

## Alpha diversity resultats article ----
### Trim.2 ----
effectif_t2 <- bdd_alpha_t2 %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_SpecRich_5000_ASV_Y1))%>%
  select(all_of(pollutant_vec_t2)) %>%
  na.omit %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

uni_t2_specrich <- model(data = bdd_alpha_t2,
                         outcome = ch_feces_SpecRich_5000_ASV_Y1,
                         exposure_vec = pollutant_vec_t2,
                         digit_beta_IC = 1)

uni_t2_shannon <- model(data = bdd_alpha_t2,
                        outcome = ch_feces_Shannon_5000_ASV_Y1,
                        exposure_vec = pollutant_vec_t2,
                        digit_beta_IC = 2)

uni_t2_faith <- model(data = bdd_alpha_t2,
                      outcome = ch_feces_Faith_5000_ASV_Y1,
                      exposure_vec = pollutant_vec_t2,
                      digit_beta_IC = 1)

mixture_t2_specrich <-
  mixture_t2_function(
    outcome = bdd_alpha_t2$ch_feces_SpecRich_5000_ASV_Y1,
    data =  bdd_alpha_t2)

mixture_t2_shannon <-
  mixture_t2_function(
    outcome = bdd_alpha_t2$ch_feces_Shannon_5000_ASV_Y1,
    data =  bdd_alpha_t2)

mixture_t2_faith <- 
  mixture_t2_function(
    outcome = bdd_alpha_t2$ch_feces_Faith_5000_ASV_Y1,
    data =  bdd_alpha_t2)


mixture_t2_specrich <- mixture_t2_specrich %>%
  tbl_regression(
    include = all_of(pollutant_vec_t2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p() %>%
  bold_labels()

mixture_t2_shannon <- mixture_t2_shannon %>%
  tbl_regression(
    include = all_of(pollutant_vec_t2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()

mixture_t2_faith <- mixture_t2_faith %>%
  tbl_regression(
    include = all_of(pollutant_vec_t2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()

mixture_t2 <- tbl_merge(tbls = list(effectif_t2, 
                                    uni_t2_specrich, mixture_t2_specrich, 
                                    uni_t2_shannon, mixture_t2_shannon, 
                                    uni_t2_faith, mixture_t2_faith), 
                        tab_spanner = c("**Effectif**", 
                                        "**Uni-pollutant model, Specific richness**", "**Mixture model, Specific richness**", 
                                        "**Uni-pollutant model, Shannon diversity**", "**Mixture model, Shannon diversity**", 
                                        "**Uni-pollutant model, Faith phylogenetic diversity**", "**Mixture model, Faith phylogenetic diversity**"))


### Trim.3 ----
effectif_t3 <- bdd_alpha_t3 %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_SpecRich_5000_ASV_Y1))%>%
  select(all_of(pollutant_vec_t3)) %>%
  na.omit %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

uni_t3_specrich <- model(data = bdd_alpha_t3,
                         outcome = ch_feces_SpecRich_5000_ASV_Y1,
                         exposure_vec = pollutant_vec_t3,
                         digit_beta_IC = 1)

uni_t3_shannon <- model(data = bdd_alpha_t3,
                        outcome = ch_feces_Shannon_5000_ASV_Y1,
                        exposure_vec = pollutant_vec_t3,
                        digit_beta_IC = 2)

uni_t3_faith <- model(data = bdd_alpha_t3,
                      outcome = ch_feces_Faith_5000_ASV_Y1,
                      exposure_vec = pollutant_vec_t3,
                      digit_beta_IC = 1)

mixture_t3_specrich <- 
  mixture_t3_function(
    outcome = bdd_alpha_t3$ch_feces_SpecRich_5000_ASV_Y1, 
    data = bdd_alpha_t3) 

mixture_t3_shannon <- 
  mixture_t3_function(
    outcome = bdd_alpha_t3$ch_feces_Shannon_5000_ASV_Y1, 
    data = bdd_alpha_t3) 

mixture_t3_faith <- 
  mixture_t3_function(
    outcome = bdd_alpha_t3$ch_feces_Faith_5000_ASV_Y1, 
    data = bdd_alpha_t3) 

mixture_t3_specrich <- mixture_t3_specrich %>%
  tbl_regression(
    include = all_of(pollutant_vec_t3),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p() %>%
  bold_labels()

mixture_t3_shannon <- mixture_t3_shannon %>%
  tbl_regression(
    include = all_of(pollutant_vec_t3),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()

mixture_t3_faith <- mixture_t3_faith %>%
  tbl_regression(
    include = all_of(pollutant_vec_t3),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()


mixture_t3 <- tbl_merge(tbls = list(effectif_t3, 
                                    uni_t3_specrich, mixture_t3_specrich, 
                                    uni_t3_shannon, mixture_t3_shannon, 
                                    uni_t3_faith, mixture_t3_faith), 
                        tab_spanner = c("**Effectif**", 
                                        "**Uni-pollutant model, Specific richness**", "**Mixture model, Specific richness**", 
                                        "**Uni-pollutant model, Shannon diversity**", "**Mixture model, Shannon diversity**", 
                                        "**Uni-pollutant model, Faith phylogenetic diversity**", "**Mixture model, Faith phylogenetic diversity**"))

### 2 months ----
effectif_M2 <- bdd_alpha_M2 %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_SpecRich_5000_ASV_Y1))%>%
  select(all_of(pollutant_vec_M2)) %>%
  na.omit %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

uni_M2_specrich <- model(data = bdd_alpha_M2,
                         outcome = ch_feces_SpecRich_5000_ASV_Y1,
                         exposure_vec = pollutant_vec_M2,
                         digit_beta_IC = 1)

uni_M2_shannon <- model(data = bdd_alpha_M2,
                        outcome = ch_feces_Shannon_5000_ASV_Y1,
                        exposure_vec = pollutant_vec_M2,
                        digit_beta_IC = 2)

uni_M2_faith <- model(data = bdd_alpha_M2,
                      outcome = ch_feces_Faith_5000_ASV_Y1,
                      exposure_vec = pollutant_vec_M2,
                      digit_beta_IC = 1)


mixture_M2_specrich <- 
  mixture_M2_function(
    outcome = bdd_alpha_M2$ch_feces_SpecRich_5000_ASV_Y1, 
    data = bdd_alpha_M2)
mixture_M2_shannon <- 
  mixture_M2_function(
    outcome = bdd_alpha_M2$ch_feces_Shannon_5000_ASV_Y1, 
    data = bdd_alpha_M2)
mixture_M2_faith <- 
  mixture_M2_function(
    outcome = bdd_alpha_M2$ch_feces_Faith_5000_ASV_Y1, 
    data = bdd_alpha_M2)

mixture_M2_specrich <- mixture_M2_specrich %>%
  tbl_regression(
    include = all_of(pollutant_vec_M2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p() %>%
  bold_labels()

mixture_M2_shannon <- mixture_M2_shannon %>%
  tbl_regression(
    include = all_of(pollutant_vec_M2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()

mixture_M2_faith <- mixture_M2_faith %>%
  tbl_regression(
    include = all_of(pollutant_vec_M2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()


mixture_M2 <- tbl_merge(tbls = list(effectif_M2, 
                                    uni_M2_specrich, mixture_M2_specrich, 
                                    uni_M2_shannon, mixture_M2_shannon, 
                                    uni_M2_faith, mixture_M2_faith), 
                        tab_spanner = c("**Effectif**", 
                                        "**Uni-pollutant model, Specific richness**", "**Mixture model, Specific richness**", 
                                        "**Uni-pollutant model, Shannon diversity**", "**Mixture model, Shannon diversity**", 
                                        "**Uni-pollutant model, Faith phylogenetic diversity**", "**Mixture model, Faith phylogenetic diversity**"))


### 12 months ----
effectif_Y1 <- bdd_alpha_Y1 %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_SpecRich_5000_ASV_Y1))%>%
  select(all_of(pollutant_vec_Y1)) %>%
  na.omit %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

uni_Y1_specrich <- model(data = bdd_alpha_Y1,
                         outcome = ch_feces_SpecRich_5000_ASV_Y1,
                         exposure_vec = pollutant_vec_Y1,
                         digit_beta_IC = 1)

uni_Y1_shannon <- model(data = bdd_alpha_Y1,
                        outcome = ch_feces_Shannon_5000_ASV_Y1,
                        exposure_vec = pollutant_vec_Y1,
                        digit_beta_IC = 2)

uni_Y1_faith <- model(data = bdd_alpha_Y1,
                      outcome = ch_feces_Faith_5000_ASV_Y1,
                      exposure_vec = pollutant_vec_Y1,
                      digit_beta_IC = 1)

mixture_Y1_specrich <- 
  mixture_Y1_function(
    outcome = bdd_alpha_Y1$ch_feces_SpecRich_5000_ASV_Y1, 
    data = bdd_alpha_Y1) 

mixture_Y1_shannon <- 
  mixture_Y1_function(
    outcome = bdd_alpha_Y1$ch_feces_Shannon_5000_ASV_Y1, 
    data = bdd_alpha_Y1) 

mixture_Y1_faith <- 
  mixture_Y1_function(
    outcome = bdd_alpha_Y1$ch_feces_Faith_5000_ASV_Y1, 
    data = bdd_alpha_Y1) 

mixture_Y1_specrich <- mixture_Y1_specrich %>%
  tbl_regression(
    include = all_of(pollutant_vec_Y1),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p() %>%
  bold_labels()

mixture_Y1_shannon <- mixture_Y1_shannon %>%
  tbl_regression(
    include = all_of(pollutant_vec_Y1),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()


mixture_Y1_faith <- mixture_Y1_faith %>%
  tbl_regression(
    include = all_of(pollutant_vec_Y1),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE)%>%
  bold_p() %>%
  bold_labels()


mixture_Y1 <- tbl_merge(tbls = list(effectif_Y1, 
                                    uni_Y1_specrich, mixture_Y1_specrich, 
                                    uni_Y1_shannon, mixture_Y1_shannon, 
                                    uni_Y1_faith, mixture_Y1_faith), 
                        tab_spanner = c("**Effectif**", 
                                        "**Uni-pollutant model, Specific richness**", "**Mixture model, Specific richness**", 
                                        "**Uni-pollutant model, Shannon diversity**", "**Mixture model, Shannon diversity**", 
                                        "**Uni-pollutant model, Faith phylogenetic diversity**", "**Mixture model, Faith phylogenetic diversity**"))

## Taxonomy ----
### Trim.2 ----
effectif_t2 <- bdd_taxa %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_rel_g4_Y1))%>%
  select(all_of(pollutant_vec_t2)) %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

mixture_t2_g1 <- 
  mixture_t2_function(
    outcome = bdd_taxa$ch_feces_rel_g1_Y1, 
    data = bdd_taxa) 
mixture_t2_g2 <- 
  mixture_t2_function(
    outcome = bdd_taxa$ch_feces_rel_g2_Y1, 
    data = bdd_taxa) 
mixture_t2_g3 <- 
  mixture_t2_function(
    outcome = bdd_taxa$ch_feces_rel_g3_Y1, 
    data = bdd_taxa) 
mixture_t2_g4 <- 
  mixture_t2_function(
    outcome = bdd_taxa$ch_feces_rel_g4_Y1, 
    data = bdd_taxa) 

mixture_t2_g1 <- mixture_t2_g1%>%
  tbl_regression(
    include = all_of(pollutant_vec_t2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t2_g2 <- mixture_t2_g2 %>%
  tbl_regression(
    include = all_of(pollutant_vec_t2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t2_g3 <- mixture_t2_g3 %>%
  tbl_regression(
    include = all_of(pollutant_vec_t2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t2_g4 <- mixture_t2_g4 %>%
  tbl_regression(
    include = all_of(pollutant_vec_t2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t2 <- 
  tbl_merge(tbls = list(effectif_t2, mixture_t2_g1, mixture_t2_g2, mixture_t2_g3, mixture_t2_g4),
            tab_spanner = c("**Effectif**", 
                            "**Bifidobacterium**",
                            "**Bacteroides**",
                            "**Blautia",
                            "**Escherichia and Shigella genera**"))


### Trim.3 ----
effectif_t3 <- bdd_taxa %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_rel_g4_Y1))%>%
  select(all_of(pollutant_vec_t3)) %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

mixture_t3_g1 <- 
  mixture_t3_function(
    outcome = bdd_taxa$ch_feces_rel_g1_Y1, 
    data = bdd_taxa) 
mixture_t3_g2 <- 
  mixture_t3_function(
    outcome = bdd_taxa$ch_feces_rel_g2_Y1, 
    data = bdd_taxa) 
mixture_t3_g3 <- 
  mixture_t3_function(
    outcome = bdd_taxa$ch_feces_rel_g3_Y1, 
    data = bdd_taxa) 
mixture_t3_g4 <- 
  mixture_t3_function(
    outcome = bdd_taxa$ch_feces_rel_g4_Y1, 
    data = bdd_taxa) 

mixture_t3_g1 <- mixture_t3_g1%>%
  tbl_regression(
    include = all_of(pollutant_vec_t3),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t3_g2 <- mixture_t3_g2 %>%
  tbl_regression(
    include = all_of(pollutant_vec_t3),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t3_g3 <- mixture_t3_g3 %>%
  tbl_regression(
    include = all_of(pollutant_vec_t3),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t3_g4 <- mixture_t3_g4 %>%
  tbl_regression(
    include = all_of(pollutant_vec_t3),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_t3 <- 
  tbl_merge(tbls = list(effectif_t3, mixture_t3_g1, mixture_t3_g2, mixture_t3_g3, mixture_t3_g4),
            tab_spanner = c("**Effectif**", 
                            "**Bifidobacterium**",
                            "**Bacteroides**",
                            "**Blautia",
                            "**Escherichia and Shigella genera**"))



### 2 months ----
effectif_M2 <- bdd_taxa %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_rel_g4_Y1))%>%
  select(all_of(pollutant_vec_M2)) %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

mixture_M2_g1 <- 
  mixture_M2_function(
    outcome = bdd_taxa$ch_feces_rel_g1_Y1, 
    data = bdd_taxa) 
mixture_M2_g2 <- 
  mixture_M2_function(
    outcome = bdd_taxa$ch_feces_rel_g2_Y1, 
    data = bdd_taxa) 
mixture_M2_g3 <- 
  mixture_M2_function(
    outcome = bdd_taxa$ch_feces_rel_g3_Y1, 
    data = bdd_taxa) 
mixture_M2_g4 <- 
  mixture_M2_function(
    outcome = bdd_taxa$ch_feces_rel_g4_Y1, 
    data = bdd_taxa) 

mixture_M2_g1 <- mixture_M2_g1%>%
  tbl_regression(
    include = all_of(pollutant_vec_M2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_M2_g2 <- mixture_M2_g2 %>%
  tbl_regression(
    include = all_of(pollutant_vec_M2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_M2_g3 <- mixture_M2_g3 %>%
  tbl_regression(
    include = all_of(pollutant_vec_M2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_M2_g4 <- mixture_M2_g4 %>%
  tbl_regression(
    include = all_of(pollutant_vec_M2),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_M2 <- 
  tbl_merge(tbls = list(effectif_M2, mixture_M2_g1, mixture_M2_g2, mixture_M2_g3, mixture_M2_g4),
            tab_spanner = c("**Effectif**", 
                            "**Bifidobacterium**",
                            "**Bacteroides**",
                            "**Blautia",
                            "**Escherichia and Shigella genera**"))


### 12 months ----
effectif_Y1 <- bdd_taxa %>%                   # Création d'une colonne effectif
  filter(!is.na(ch_feces_rel_g4_Y1))%>%
  select(all_of(pollutant_vec_Y1)) %>%
  tbl_summary(missing = "no", 
              statistic = all_continuous() ~ "{N_nonmiss}") %>%
  bold_labels()

mixture_Y1_g1 <- 
  mixture_Y1_function(
    outcome = bdd_taxa$ch_feces_rel_g1_Y1, 
    data = bdd_taxa) 
mixture_Y1_g2 <- 
  mixture_Y1_function(
    outcome = bdd_taxa$ch_feces_rel_g2_Y1, 
    data = bdd_taxa) 
mixture_Y1_g3 <- 
  mixture_Y1_function(
    outcome = bdd_taxa$ch_feces_rel_g3_Y1, 
    data = bdd_taxa) 
mixture_Y1_g4 <- 
  mixture_Y1_function(
    outcome = bdd_taxa$ch_feces_rel_g4_Y1, 
    data = bdd_taxa) 

mixture_Y1_g1 <- mixture_Y1_g1%>%
  tbl_regression(
    include = all_of(pollutant_vec_Y1),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_Y1_g2 <- mixture_Y1_g2 %>%
  tbl_regression(
    include = all_of(pollutant_vec_Y1),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_Y1_g3 <- mixture_Y1_g3 %>%
  tbl_regression(
    include = all_of(pollutant_vec_Y1),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_Y1_g4 <- mixture_Y1_g4 %>%
  tbl_regression(
    include = all_of(pollutant_vec_Y1),
    pvalue_fun = ~ style_pvalue(.x, digits = 2),
    estimate_fun = ~ style_sigfig(.x, digits = 1)
  ) %>%
  add_global_p(type = "II", keep = TRUE) %>%
  bold_p(t=0.1) %>%
  bold_labels()

mixture_Y1 <- 
  tbl_merge(tbls = list(effectif_Y1, mixture_Y1_g1, mixture_Y1_g2, mixture_Y1_g3, mixture_Y1_g4),
            tab_spanner = c("**Effectif**", 
                            "**Bifidobacterium**",
                            "**Bacteroides**",
                            "**Blautia",
                            "**Escherichia and Shigella genera**"))


