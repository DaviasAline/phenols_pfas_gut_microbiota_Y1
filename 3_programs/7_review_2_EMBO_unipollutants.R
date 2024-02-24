# Review #2 for EMBO mol med
# A. Davias
# 04/07/2023

## Chargement des données et des fonctions existantes ----
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')  
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)
library(broom)
library(forestplot)
library(bkmr)
library(fields)
library(psych)


## Création de fonctions ----

# Fonction pour obtenir les résultats bruts
lm_func <- function(outcome, exposure, data){   
  model <- lm({{outcome}} ~
                exposure +
                ch_feces_RUN_Y1 +
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
  return(model)
}

# Fonction pour le nettoyage de résultats bruts
data_prep <- function(results_list, outcome_name) {
  results_list <- map(results_list, broom::tidy, conf.int = TRUE) %>%
    bind_rows(.id = "model_id") %>%
    filter(!term %in%
             c("(Intercept)",
               "ch_feces_RUN_Y1R3",
               "ch_feces_age_w_Y1_i",
               "po_delmodVaginal delivery",
               "ch_food_intro_Y1_3cat_iBetween 6 and 12 months old",
               "ch_food_intro_Y1_3cat_iNot introduced at 12 months old",
               "ch_antibio_Y1_2cat_iYes",
               "mo_par_2cat1 child or more",
               "mo_pets_iOne or more",
               "ch_sexMale",
               "mo_tob_gr_anyt_yn_n2_iYes",
               "Mo_ETS_anyT_yn1_opt_iYes",
               "ch_ETS_12m_opt36mYes",
               "mo_interpreg_3cat2 years and more",
               "mo_interpreg_3catPrimiparous",
               "mo_dipl_3cat_i3-4years after graduation",
               "mo_dipl_3cat_i>=5years after graduation",
               "po_w_kg_3cat3-3.4 Kg",
               "po_w_kg_3cat>= 3.5 Kg",
               "po_he_3cat_i50-51 cm",
               "po_he_3cat_i>= 52 cm",
               "ch_w_Y1_3cat_i8.5-9.9 Kg",
               "ch_w_Y1_3cat_i>=10 Kg",
               "ch_he_Y1_3cat_i75-77.9 cm",
               "ch_he_Y1_3cat_i>=78 cm",
               "po_gd",
               "mo_age",
               "mo_bmi_bepr_3cat_i19-23.9 Kg/m2",
               "mo_bmi_bepr_3cat_i>=24 Kg/m2",
               "bf_duration_till48w_4cat_i<24 weeks",
               "bf_duration_till48w_4cat_i24-47 weeks",
               "bf_duration_till48w_4cat_iStill breastfeed at 48 weeks"))%>%
    mutate(
      model_type = as.factor("adjusted"),
      outcome_name = outcome_name,
      term = str_replace_all(term, "exposure", "")) %>% 
    rename(outcome = outcome_name) %>%
    rename(exposure = model_id) %>%
    select(model_type,
           outcome,
           exposure,
           term,
           everything())
  
  return(results_list)
}

# Fonction pour obtenir les p-values globales pour les variables catégorielles dans les résultats bruts
global_p_cat <- function(results_list, outcome_name) {
  
  global_p_cat <- function(results_list) {
    Anova(
      mod = results_list, 
      type="III",
      singular.ok = TRUE)
  } 
  
  results <- lapply(results_list, global_p_cat)
  
  results <- map(results, broom::tidy, conf.int = TRUE) %>%
    bind_rows(.id = "model_id") %>%
    filter(term %in% "exposure") %>%
    rename(exposure = model_id, 
           p.value_anova = p.value, 
           sumsq_anova = sumsq, 
           F.values_anova = F.values) %>%
    select(outcome_name, exposure, sumsq_anova, F.values_anova, p.value_anova)
  
  return(results)
}

# Fonction pour tableaux de résultats article
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

effectif_column <- function(data, outcome, exposure_vec) {
  data %>%                   
    filter(!is.na({{outcome}}))%>%
    select(all_of({{exposure_vec}})) %>%
    tbl_summary(missing = "no", 
                statistic = all_continuous() ~ "{N_nonmiss}") %>%
    bold_labels()
}

# Fonction pour figures forestplots à partir du tableau bruts des résultats 
forestplot <- function(results_list, outcome_name) {
  results <- results_list %>%
    ggplot(aes(x = exposure, 
               y = estimate, 
               min = conf.low, 
               ymax = conf.high, 
               #color = interaction(exposure_window, term_2), 
               color = term_rec, 
               shape = p_value_shape)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.4) +
    labs(x = "Exposures", y = outcome_name) +
    theme_bw() +
    coord_flip()  +
    scale_shape_manual(values = c(19, 8),
                       name = "p-value") +
    guides(color = guide_legend(title = ""))+
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7), 
          legend.position = "bottom",
          legend.box = "vertical", 
          legend.justification = "right", 
          legend.spacing.y = unit(0, "cm"), 
          legend.spacing.x = unit(0, "cm"), 
          legend.box.margin = margin(0,0,0,0, "cm"), 
          legend.margin = margin(0,0,0,0, "cm"))
  
  return(results)
}

# Fonction pour voir la table de correlation avant de faire la correction pour compararaison multiple selon correlation entre les tests
cor_mixed <- function(data,method="spearman"){
  ## Purpose: 2 -  use the function mixedCor from psych to calculate a 
  #                correlation matrix on mixed type variables (continuous and categorical).
  #                Note: mixedCor consider categorical variable as ordered factors.
  ## Inputs: - data: dataframe (not mice) with only the variables to consider 
  ##                 (mixed type allowed)
  ## Output: - correlation matrix
  
  ## STEP 1 : TYPE OF VARIABLES
  continuous_var = which(sapply(data, class) == "numeric")
  names(continuous_var)=NULL
  # Categorical var
  categorical_var = which(sapply(data, class) == "factor")
  #  - 2 levels only
  binary_var = categorical_var[sapply(data[,categorical_var],nlevels)==2] 
  binary_var = binary_var[!is.na(binary_var)]
  names(binary_var)=NULL
  #  - More than 2 levels (but less than 8)
  poly_var = categorical_var[(sapply(data[,categorical_var],nlevels)>2 & sapply(data[,categorical_var],nlevels)<8)] %>% na.exclude()
  names(poly_var)=NULL
  
  ## STEP 2 : CORRELATION MATRIX USING MIXEDCOR FUNCTION (FROM PSYCH)
  # data converted in numeric (necessary)
  data[,] = lapply(data[,],as.numeric)
  # Correlation matrix
  cor = data  %>% 
    mixedCor(c=continuous_var,p=poly_var,d=binary_var,use="pairwise.complete.obs",method=method)%>% pluck('rho')
  return(cor)
}

# Fonction pour voir l'alpha corrigé pour la correction pour comparaison multiple 
alpha_corrected <- function(data, alpha=0.05) {
  # Purpose: Alpha-risk correction (FWER), inspired from https://www-ncbi-nlm-nih-gov.gate2.inist.fr/pmc/articles/PMC3325408/
  # Inputs
  # - data = dataset with all exposures to consider
  # - alpha = risk to correct
  # Output
  # - alpha corrected
  
  M <- ncol(data)
  
  ## FIRST STEP: CORRELATION MATRIX (see function "cor_mixed" in the path "general_functions.R")
  cor = cor_mixed(data) #handle mixed typed variables
  
  ## SECOND STEP : EIGENVALUES OF CORRELATION MATRIX
  lambdas <- base::eigen(cor)$values
  
  ## THIRD STEP: CORRECT ALPHA
  M0 <- M - sum( (lambdas > 1) * (lambdas - 1))
  M0 <- floor(M0) # always round down estimated value
  alpha_corrected <- alpha/M0
  
  return(alpha_corrected)
  
}

# Fonction pour voir le nombre de tests rééls après prise en compte de la corrélation entre les expositions 
# Peut aussi s'appliquer à la corrélation entre les outcomes si besoin 
# puis faire nombre d'expo x nombre d'outcome
# multiplier ce nombre à la p-value
# = cette correction pour comparaison multiple adapté de Li et al 2012 est une correction de Bonferroni qui prend en compte la corrélation entre les tests effectués
M0_corrected <- function(data, alpha=0.05) {
  # Purpose: Alpha-risk correction (FWER), inspired from https://www-ncbi-nlm-nih-gov.gate2.inist.fr/pmc/articles/PMC3325408/
  # Inputs
  # - data = dataset with all exposures to consider
  # - alpha = risk to correct
  # Output
  # - alpha corrected
  
  M <- ncol(data)
  
  ## FIRST STEP: CORRELATION MATRIX (see function "cor_mixed" in the path "general_functions.R")
  cor = cor_mixed(data, method = "pearson") #handle mixed typed variables
  
  ## SECOND STEP : EIGENVALUES OF CORRELATION MATRIX
  lambdas <- base::eigen(cor)$values
  
  ## THIRD STEP: CORRECT ALPHA
  M0 <- M - sum( (lambdas > 1) * (lambdas - 1))
  M0 <- floor(M0) # always round down estimated value
  alpha_corrected <- alpha/M0
  
  return(M0)
  
}



## Création des vecteurs ----
conf_vec <- bdd_alpha %>% 
  select(all_of(phenols_vec_2))%>% 
  select(contains(c("MEPA", "ETPA", "PRPA", "BUPA", "TRCS"))) %>% 
  colnames()
explo_vec <- bdd_alpha %>% 
  select(all_of(phenols_vec_2), all_of(pfas_vec)) %>% 
  select(!contains(c("MEPA", "ETPA", "PRPA", "BUPA", "TRCS"))) %>% 
  colnames()
pollutants_vec <- metadata %>% 
  select(all_of(phenols_vec_2), all_of(pfas_vec)) %>% 
  colnames()
alpha_vec <- bdd_alpha %>% 
  select(all_of(alpha_vec)) %>% 
  select(contains("5000")) %>% 
  colnames() 



## Results - Unipolluants analysis ----
### Tableaux bruts ----
#### Colonnes principales ----
results_list_1 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_SpecRich_5000_ASV_Y1, data = bdd_alpha)
results_list_2 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_Shannon_5000_ASV_Y1, data = bdd_alpha)
results_list_3 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_Faith_5000_ASV_Y1, data = bdd_alpha)

results_list_4 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p1_Y1, data = bdd_taxa)
results_list_5 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p2_Y1, data = bdd_taxa)
results_list_6 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p3_Y1, data = bdd_taxa)
results_list_7 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p4_Y1, data = bdd_taxa)

results_list_8 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g1_Y1, data = bdd_taxa)
results_list_9 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g2_Y1, data = bdd_taxa)
results_list_10 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g3_Y1, data = bdd_taxa)
results_list_11 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g4_Y1, data = bdd_taxa)

results_multi <- list(
  
  data_prep(results_list = results_list_1, outcome_name = "Specific richness"),
  data_prep(results_list = results_list_2, outcome_name = "Shannon diversity"),
  data_prep(results_list = results_list_3, outcome_name = "Faith phylogenetic diversity"),
  
  data_prep(results_list = results_list_4, outcome_name = "Firmicutes"),
  data_prep(results_list = results_list_5, outcome_name = "Actinobacteria"),
  data_prep(results_list = results_list_6, outcome_name = "Bacteroidetes"),
  data_prep(results_list = results_list_7, outcome_name = "Proteobacteria"),
  
  data_prep(results_list = results_list_8, outcome_name = "Bifidobacterium"),
  data_prep(results_list = results_list_9, outcome_name = "Bacteroides"),
  data_prep(results_list = results_list_10, outcome_name = "Blautia"),
  data_prep(results_list = results_list_11, outcome_name = "Escherichia and Shigella")) %>%
  
  bind_rows() %>%
  mutate(
    exposure = as.factor(exposure)) %>%
  mutate(
    exposure_type = case_when(
      exposure == "mo_BUPA_total_cat_t2" ~ "categorical_3",
      exposure == "mo_BUPA_total_cat_t3" ~ "categorical_3",
      exposure == "ch_BUPA_total_cat_Y1" ~ "categorical_3",
      exposure == "ch_BPS_total_cat_Y1" ~ "categorical_3",
      exposure == "mo_PFDoDa_cat" ~ "categorical_3",
      exposure == "mo_PFHxPA_cat" ~ "categorical_3",
      exposure == "mo_PFTrDa_cat" ~ "categorical_3",
      grepl("cat", exposure) ~ "categorical_2",
      TRUE ~ "continuous"), 
    exposure = str_replace_all(exposure,
                               c("mo_" = "",
                                 "ch_" = "",
                                 "_total_i_cor_" = " ",
                                 "_i_cor_" = " ", 
                                 "_ln" = "",
                                 "_conj" = "",
                                 "ln" = "",
                                 "_cor" = "",
                                 "MEPA" = "Methylparaben", 
                                 "ETPA" = "Ethylparaben", 
                                 "PRPA" = "Propylparaben", 
                                 "BUPA" = "Butylparaben", 
                                 "BPA" = "Bisphenol A", 
                                 "OXBE" = "Benzophenone 3", 
                                 "TRCS" = "Triclosan")), 
    exposure = fct_recode(exposure,
                          "6_2diPAP t2" = "6_2diPAP_cat_2",
                          "8_2diPAP t2" = "8_2diPAP_cat_2",
                          "Bisphenol S M2" = "BPS_total_cat_M2_2",
                          "Bisphenol S t2" = "BPS_total_cat_t2_2",
                          "Bisphenol S t3" = "BPS_total_cat_t3_2",
                          "Bisphenol S Y1" = "BPS_total_cat_Y1",
                          "Butylparaben M2" = "Butylparaben_total_cat_M2_2",
                          "Butylparaben t2" = "Butylparaben_total_cat_t2",
                          "Butylparaben t3" = "Butylparaben_total_cat_t3",
                          "Butylparaben Y1" = "Butylparaben_total_cat_Y1",
                          "Ethylparaben M2" = "Ethylparaben_total_cat_M2_2",
                          "PFBS t2" = "PFBS_cat_2",
                          "PFDA t2" = "PFDA ",
                          "PFDoDa t2" = "PFDoDa_cat",
                          "PFHpA t2" = "PFHpA_cat_2",
                          "PFHpS t2" = "PFHpS ",
                          "PFHxPA t2" = "PFHxPA_cat",
                          "PFHxS t2" = "PFHxS",
                          "PFNA t2" = "PFNA",
                          "PFOA t2" = "PFOA",
                          "PFOS t2" = "PFOS",
                          "PFOSA t2" = "PFOSA_cat_2",
                          "PFTrDa t2" = "PFTrDa_cat",
                          "PFUnDA t2" = "PFUnDA ",
                          "Propylparaben M2" = "Propylparaben_total_cat_M2_2"), 
    exposure = fct_relevel(exposure,
                           "8_2diPAP t2", "6_2diPAP t2", "PFOSA t2", "PFBS t2", "PFTrDa t2",
                           "PFHpA t2", "PFHxPA t2", "PFDoDa t2", "PFHpS t2", "PFUnDA t2",
                           "PFDA t2", "PFOS t2", "PFHxS t2", "PFNA t2", "PFOA t2", "Benzophenone 3 Y1",
                           "Benzophenone 3 M2", "Benzophenone 3 t3", "Benzophenone 3 t2",
                           "Bisphenol S Y1", "Bisphenol S M2", "Bisphenol S t3", "Bisphenol S t2",
                           "Bisphenol A Y1", "Bisphenol A M2", "Bisphenol A t3", "Bisphenol A t2",
                           "Triclosan Y1", "Triclosan M2", "Triclosan t3", "Triclosan t2",
                           "Butylparaben Y1", "Butylparaben M2", "Butylparaben t3", "Butylparaben t2",
                           "Propylparaben Y1", "Propylparaben M2", "Propylparaben t3", "Propylparaben t2",
                           "Ethylparaben Y1", "Ethylparaben M2", "Ethylparaben t3", "Ethylparaben t2",
                           "Methylparaben Y1", "Methylparaben M2", "Methylparaben t3", "Methylparaben t2"), 
    p_value_shape = case_when(p.value < 0.1 ~ "p.value <0.1",
                              p.value > 0.1 ~ "p.value >0.1"), 
    p_value_shape = fct_relevel(p_value_shape,
                                "p.value >0.1", 
                                "p.value <0.1"),
    exposure_window = case_when(grepl("t2", exposure) ~ "Trim.2", 
                                grepl("t3", exposure) ~ "Trim.3", 
                                grepl("M2", exposure) ~ "2 months", 
                                grepl("Y1", exposure) ~ "12 months", 
                                TRUE ~ "Trim.2"), 
    exposure_window = fct_relevel(exposure_window, 
                                  "Trim.2", 
                                  "Trim.3", 
                                  "2 months", 
                                  "12 months"),
    exposure_family = case_when(grepl("paraben", exposure) ~ "Phenols", 
                                grepl("Triclosan", exposure) ~ "Phenols", 
                                grepl("Bisphenol", exposure) ~ "Phenols", 
                                grepl("Benzophenone", exposure) ~ "Phenols", 
                                TRUE ~ "PFAS"), 
    analysis = case_when(grepl("paraben", exposure) ~ "confirmatory", 
                         grepl("Triclosan", exposure) ~ "confirmatory", 
                         TRUE ~ "exploratory"), 
    term_rec = case_when(exposure_type == "continuous" ~ "Continuous", 
                         exposure_type == "categorical_2" & term == ">LOQ" ~ ">LOQ, compared to <LOQ", 
                         exposure_type == "categorical_2" & term == ">LOD" ~ ">LOD, compared to <LOD", 
                         exposure_type == "categorical_3" & term == "LOD-LOQ" ~ "LOD-LOQ, compared to <LOD",
                         exposure_type == "categorical_3" & term == ">LOQ" ~ ">LOQ, compared to <LOD"),
    term_rec = fct_relevel(term_rec,
                           "Continuous", ">LOQ, compared to <LOQ", ">LOD, compared to <LOD",
                           "LOD-LOQ, compared to <LOD", ">LOQ, compared to <LOD"),
    row_to_choose = case_when(term == ">LOQ" & exposure == "Butylparaben t2" ~ "no",
                              term == ">LOQ" & exposure == "Butylparaben t3" ~ "no",
                              term == ">LOQ" & exposure == "Butylparaben Y1" ~ "no",
                              term == ">LOQ" & exposure == "Bisphenol S Y1" ~ "no",
                              term == ">LOQ" & exposure == "PFDoDa t2" ~ "no",
                              term == ">LOQ" & exposure == "PFHxPA t2" ~ "no",
                              term == ">LOQ" & exposure == "PFTrDa t2" ~ "no",
                              TRUE ~ "yes")) %>%
  select("model_type", 
         "analysis",
         "outcome", 
         "exposure",
         "exposure_type",
         "exposure_window",
         "exposure_family",
         "term",
         "term_rec",
         "row_to_choose",
         "estimate",
         "std.error",
         "statistic",
         "conf.low",
         "conf.high",
         "p.value",   
         "p_value_shape")


#### Correction pour comparaison multiple ----
##### Version avec tous les tests ----
# On regarde la table de correlation des expositions et des outcomes 
# la focntion MO_corrected estime le nombre rééls de tests réalisés après prise en compte de la corrélation
# puis on multiplie la p-value x le nombre d'outcome donné x le nombre d'expositions données 

test_expo <- bdd_alpha %>% 
  filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
  select(all_of(pollutants_vec)) %>% na.omit()

test_outcome <- bdd_alpha %>% 
  filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
  select(ident, 
         ch_feces_Shannon_5000_ASV_Y1, 
         ch_feces_Faith_5000_ASV_Y1, 
         ch_feces_SpecRich_5000_ASV_Y1) %>% na.omit()

test_outcome_2 <- bdd_taxa %>% 
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(ident, 
         ch_feces_rel_p1_Y1:ch_feces_rel_p4_Y1,
         ch_feces_rel_g1_Y1:ch_feces_rel_g4_Y1) %>% na.omit()

test_outcome <- left_join(test_outcome, test_outcome_2, by = "ident") %>% select(-ident)

var_lab(test_outcome$ch_feces_Shannon_5000_ASV_Y1) <- NULL 
var_lab(test_outcome$ch_feces_Faith_5000_ASV_Y1) <- NULL 
var_lab(test_outcome$ch_feces_SpecRich_5000_ASV_Y1) <- NULL 

var_lab(test_outcome$ch_feces_rel_p1_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p2_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p3_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p4_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g1_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g2_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g3_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g4_Y1) <- NULL 

cor_mixed_table_expo <- cor_mixed(data = test_expo, method = "pearson")
results_alpha_corrected_expo <- alpha_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo <- M0_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo
## on passe de 47 expositions à 29 expositions après prise en compte de leur corrélation

cor_mixed_table_outcome <- cor_mixed(data = test_outcome, method = "pearson")
results_alpha_corrected_outcome <- alpha_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome <- M0_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome
## on passe de 47 expositions à 29 expositions après prise en compte de leur corrélation


results_multi <- results_multi %>%
  mutate(
    FWER.p.value = ifelse(row_to_choose %in% "yes", p.value * 29 * 6, NA),
    FWER.p.value = ifelse(FWER.p.value > 1, 1, FWER.p.value)) 

##### Version avec que tests exploratoires----
# On regarde la table de correlation des expositions et des outcomes 
# la focntion MO_corrected estime le nombre rééls de tests réalisés après prise en compte de la corrélation
# puis on multiplie la p-value x le nombre d'outcome donné x le nombre d'expositions données 

test_expo <- bdd_alpha %>% 
  filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
  select(all_of(pollutants_vec)) %>% 
  select(-contains("MEPA"), 
         -contains("ETPA"),
         -contains("PRPA"),
         -contains("BUPA"), 
         -contains("TRCS")) %>% 
  na.omit()

test_outcome <- bdd_alpha %>% 
  filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
  select(ident, 
         ch_feces_Shannon_5000_ASV_Y1, 
         ch_feces_Faith_5000_ASV_Y1, 
         ch_feces_SpecRich_5000_ASV_Y1) %>% na.omit()

test_outcome_2 <- bdd_taxa %>% 
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(ident, 
         ch_feces_rel_p1_Y1:ch_feces_rel_p4_Y1,
         ch_feces_rel_g1_Y1:ch_feces_rel_g4_Y1) %>% na.omit()

test_outcome <- left_join(test_outcome, test_outcome_2, by = "ident") %>% select(-ident)

var_lab(test_outcome$ch_feces_Shannon_5000_ASV_Y1) <- NULL 
var_lab(test_outcome$ch_feces_Faith_5000_ASV_Y1) <- NULL 
var_lab(test_outcome$ch_feces_SpecRich_5000_ASV_Y1) <- NULL 

var_lab(test_outcome$ch_feces_rel_p1_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p2_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p3_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p4_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g1_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g2_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g3_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g4_Y1) <- NULL 

cor_mixed_table_expo <- cor_mixed(data = test_expo, method = "pearson")
results_alpha_corrected_expo <- alpha_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo <- M0_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo
## on passe de 27 expositions à 17 expositions après prise en compte de leur corrélation

cor_mixed_table_outcome <- cor_mixed(data = test_outcome, method = "pearson")
results_alpha_corrected_outcome <- alpha_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome <- M0_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome
## on passe de 11 outcomes à 6 outcomes après prise en compte de leur corrélation


results_multi <- results_multi %>%
  mutate(
    FWER.p.value = ifelse(row_to_choose %in% "yes", p.value * 17 * 6, NA),
    FWER.p.value = ifelse(FWER.p.value > 1, 1, FWER.p.value)) 


### Tableau pour Marion ----
bdd_marion <- results_multi %>%
  select(outcome, 
         Pollutant = exposure, 
         Variable_type = exposure_type, 
         Variable_categories = term_rec, 
         estimate, 
         p.value) %>%
  #separate(Pollutant, c("Pollutant", "Timing")) %>%
  mutate(
    #Timing = fct_recode(Timing, "T2" = "t2", "T3" = "t3"), 
    Variable_type = fct_recode(Variable_type, "cat_2" = "categorical_2"),
    ch_feces_specrich_5000_Y1_b = ifelse(outcome == "Specific richness", estimate, NA), 
    ch_feces_specrich_5000_Y1_p = ifelse(outcome == "Specific richness", p.value, NA), 
    ch_feces_shannon_5000_Y1_b = ifelse(outcome == "Shannon diversity", estimate, NA), 
    ch_feces_shannon_5000_Y1_p = ifelse(outcome == "Shannon diversity", p.value, NA), 
    ch_feces_faith_5000_Y1_b = ifelse(outcome == "Faith phylogenetic diversity", estimate, NA), 
    ch_feces_faith_5000_Y1_p = ifelse(outcome == "Faith phylogenetic diversity", p.value, NA), 
    
    ch_feces_firmicutes_5000_Y1_b = ifelse(outcome == "Firmicutes", estimate, NA), 
    ch_feces_firmicutes_5000_Y1_p = ifelse(outcome == "Firmicutes", p.value, NA), 
    ch_feces_actinobacteria_5000_Y1_b = ifelse(outcome == "Actinobacteria", estimate, NA), 
    ch_feces_actinobacteria_5000_Y1_p = ifelse(outcome == "Actinobacteria", p.value, NA),
    ch_feces_bacteroidetes_5000_Y1_b = ifelse(outcome == "Bacteroidetes", estimate, NA), 
    ch_feces_bacteroidetes_5000_Y1_p = ifelse(outcome == "Bacteroidetes", p.value, NA),
    ch_feces_proteobacteria_5000_Y1_b = ifelse(outcome == "Proteobacteria", estimate, NA), 
    ch_feces_proteobacteria_5000_Y1_p = ifelse(outcome == "Proteobacteria", p.value, NA),
    
    ch_feces_blautia_5000_Y1_b = ifelse(outcome == "Blautia", estimate, NA), 
    ch_feces_blautia_5000_Y1_p = ifelse(outcome == "Blautia", p.value, NA),
    ch_feces_bifidobacterium_5000_Y1_b = ifelse(outcome == "Bifidobacterium", estimate, NA), 
    ch_feces_bifidobacterium_5000_Y1_p = ifelse(outcome == "Bifidobacterium", p.value, NA),
    ch_feces_bacteroides_5000_Y1_b = ifelse(outcome == "Bacteroides", estimate, NA), 
    ch_feces_bacteroides_5000_Y1_p = ifelse(outcome == "Bacteroides", p.value, NA),
    ch_feces_escherichia_shigella_5000_Y1_b = ifelse(outcome == "Escherichia and Shigella", estimate, NA), 
    ch_feces_escherichia_shigella_5000_Y1_p = ifelse(outcome == "Escherichia and Shigella", p.value, NA)) %>%
  select(-outcome, -estimate, -p.value)

write_xlsx(
  bdd_marion, 
  "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations écrites/9. Figure_overview_SEPAGES/SEPAGES_overview_phenols_pfas.xlsx")


### Tableaux article - alpha diversité ----
# analyses confirmatoires sans correction pour les polluants avec hypothèses à priori
results_multi_alpha_conf <- tbl_merge(                         
  tbls = list(
    effectif_column(data = bdd_alpha, 
                    outcome = ch_feces_SpecRich_5000_ASV_Y1, 
                    exposure_vec = conf_vec),
    model(data = bdd_alpha, 
          outcome = ch_feces_SpecRich_5000_ASV_Y1,
          exposure_vec = conf_vec, 
          digit_beta_IC = 1),
    model(data = bdd_alpha,
          outcome = ch_feces_Shannon_5000_ASV_Y1,
          exposure_vec = conf_vec, 
          digit_beta_IC = 2),
    model(data = bdd_alpha, 
          outcome = ch_feces_Faith_5000_ASV_Y1,
          exposure_vec = conf_vec, 
          digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Specific richness**", 
                  "**Shannon diversity**", 
                  "**Faith's phylogenetic diversity**"))

# analyses exploratoires avec correction pour les polluants sans hypothèses à priori
results_multi_alpha_explo <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_alpha, 
      outcome = ch_feces_SpecRich_5000_ASV_Y1, 
      exposure_vec = explo_vec),
    model(
      data = bdd_alpha, 
      outcome = ch_feces_SpecRich_5000_ASV_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_alpha,
      outcome = ch_feces_Shannon_5000_ASV_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 2),
    model(
      data = bdd_alpha, 
      outcome = ch_feces_Faith_5000_ASV_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Specific richness**", 
                  "**Shannon diversity**", 
                  "**Faith's phylogenetic diversity**"))


### Tableaux article - taxonomie ----
# Analyses "explicatives" non corrigées pour comparaison multiple pour la taxonomie 
# pour les associations significatives pour l'alpha diversité
# pour les autres polluants non significatifs --> on met tout en analyse sup. 

results_multi_phyla_conf <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1, 
      exposure_vec = conf_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_p2_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p3_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p4_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Phylum Firmicutes**",
                  "**Phylum Actinobacteria**",
                  "**Phylum Bacteroidetes**",
                  "**Phylum Proteobacteria**"))

results_multi_genera_conf <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1, 
      exposure_vec = conf_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_g2_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g3_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g4_Y1,
      exposure_vec = conf_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Genus Bifidobacterium**",
                  "**Genus Bacteroides**",
                  "**Genus Blautia**",
                  "**Genera Eschericha / Shigella**"))


results_multi_phyla_explo <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1, 
      exposure_vec = explo_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_p2_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p3_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p4_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Phylum Firmicutes**",
                  "**Phylum Actinobacteria**",
                  "**Phylum Bacteroidetes**",
                  "**Phylum Proteobacteria**"))

results_multi_genera_explo <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1, 
      exposure_vec = explo_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_g2_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g3_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g4_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Genus Bifidobacterium**",
                  "**Genus Bacteroides**",
                  "**Genus Blautia**",
                  "**Genera Eschericha / Shigella**"))



### Tableaux article - mixture effects ----

### Figures article - forestplots ----
leg <- results_multi %>%
  filter(outcome == "Faith phylogenetic diversity") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Faith phylogenetic diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
leg <- get_legend(leg) %>% as_ggplot()

forestplot_alpha_confirmatory_1 <-
  results_multi %>%
  filter(outcome == "Specific richness") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Specific richness") +
  theme(legend.position = "none")

forestplot_alpha_confirmatory_2 <-
  results_multi %>%
  filter(outcome == "Shannon diversity") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_confirmatory_3 <- 
  results_multi %>%
  filter(outcome == "Faith phylogenetic diversity") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Faith phylogenetic diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_confirmatory <- 
  (forestplot_alpha_confirmatory_1 + forestplot_alpha_confirmatory_2 + forestplot_alpha_confirmatory_3) / leg + 
  plot_layout(heights = c(7.5, 1))

ggsave("4_output/phenols pfas/forestplot_alpha_confirmatory.tiff", 
       plot = forestplot_alpha_confirmatory, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 150,
       dpi = 300,
       limitsize = FALSE)

leg_exploratory <-  results_multi %>%
  filter(outcome == "Faith phylogenetic diversity") %>%
  filter(analysis == "exploratory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Faith phylogenetic diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
leg_exploratory <- get_legend(leg_exploratory) %>% as_ggplot()

forestplot_alpha_exploratory_1 <-
  results_multi %>%
  filter(outcome == "Specific richness") %>%
  filter(analysis == "exploratory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Specific richness") +
  theme(legend.position = "none")

forestplot_alpha_exploratory_2 <-
  results_multi %>%
  filter(outcome == "Shannon diversity") %>%
  filter(analysis == "exploratory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_exploratory_3 <- 
  results_multi %>%
  filter(outcome == "Faith phylogenetic diversity") %>%
  filter(analysis == "exploratory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Faith phylogenetic diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_exploratory <- 
  (forestplot_alpha_exploratory_1 + forestplot_alpha_exploratory_2 + forestplot_alpha_exploratory_3) / leg + 
  plot_layout(heights = c(10, 1))

ggsave("4_output/phenols pfas/forestplot_alpha_exploratory.tiff", 
       plot = forestplot_alpha_exploratory, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 180,
       dpi = 300,
       limitsize = FALSE)

### Figures article - Heatmap correlation pollutants and covariates  ----
pollutants_vec_num <- c(phenols_vec_num, pfas_vec_num)
test <- metadata %>% 
  select(all_of(covar_vec_num), 
         all_of(pollutants_vec_num))


cormat <- round(cor(test,
                    use = "pairwise.complete.obs",
                    method = "spearman"), 1)

cormat <- cormat %>%
  as.data.frame() %>%
  select(all_of(pollutants_vec_num)) %>%
  t() %>%
  as.data.frame() %>%
  select(all_of(covar_vec_num)) %>%
  as.matrix()

colnames(cormat) <- colnames(cormat) %>%
  str_replace_all(
    c("ch_feces_age_w_Y1" = "Child age (weeks)",
      "ch_antibio_Y1" = "Antibiotics use 0-12 months",
      "mo_par" = "Maternal parity",
      "po_w_kg" = "Birth weight (kg)",
      "po_he"= "Birth length (cm)", 
      "ch_w_Y1"="Weight at one year (Kg)", 
      "ch_he_Y1"="Length at one year (cm)", 
      "po_gd"= "Gestational age (weeks)", 
      "mo_age"="Maternal age before pregnancy", 
      "mo_bmi_bepr"="Maternal BMI before pregnancy", 
      "bf_duration_till48w"="Breastfeeding duration (weeks)"))

rownames(cormat) <- rownames(cormat) %>%
  str_replace_all(
    c("mo_" = "",
      "ch_" = "",
      "_total_i_cor_" = " ",
      "_i_cor_" = " ",
      "_i_cor" = "",
      "_cor" = "",
      "MEPA" = "Methylparaben", 
      "ETPA" = "Ethylparaben", 
      "PRPA" = "Propylparaben", 
      "BUPA" = "Butylparaben", 
      "BPA" = "Bisphenol A", 
      "OXBE" = "Benzophenone 3", 
      "TRCS" = "Triclosan"))

cormat_long <-
  melt(cormat, na.rm = TRUE)  # passer en df long rapidement

cormat_long$Var1 <- fct_relevel(
  cormat_long$Var1,
  "PFHpS", "PFUnDA", "PFDA", "PFOS", "PFHxS",
  "PFNA", "PFOA", "Triclosan Y1", "Triclosan M2", "Triclosan t3",
  "Triclosan t2", "Benzophenone 3 Y1", "Benzophenone 3 M2", "Benzophenone 3 t3",
  "Benzophenone 3 t2", "Bisphenol A Y1", "Bisphenol A M2", "Bisphenol A t3",
  "Bisphenol A t2", "Propylparaben Y1", "Propylparaben t3", "Propylparaben t2",
  "Ethylparaben Y1", "Ethylparaben t3", "Ethylparaben t2", "Methylparaben Y1",
  "Methylparaben M2", "Methylparaben t3", "Methylparaben t2"
)
heatmap <-                                                    # faire la heatmap
  ggplot(cormat_long, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Spearman\nCorrelation"
  ) +
  theme_minimal() +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value),
            color = "black",
            size = 2) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 8,
      hjust = 1), 
    axis.text.y = element_text(size =8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    #legend.justification = c(1, 0),
    #legend.position = c(0.4, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 5,
    barheight = 0.5,
    title.position = "top",
    title.hjust = 0.5
  ))
heatmap


ggsave(filename = "4_output/heatmap_cor_pollutants_covar.pdf",
       plot = heatmap, 
       units = "cm",
       width = 15, 
       height = 20, 
       limitsize = FALSE)

### Figures article - Heatmap correlation entre pollutants  ----
pollutants_vec_num <- c(phenols_vec_num, pfas_vec_num)
test <- metadata %>% 
  select(all_of(pollutants_vec_num))

cor.test(test$mo_PFNA, test$mo_PFDA_i_cor, method = "pearson")
cor.test(test$mo_PFNA, test$mo_PFOA_cor, method = "pearson")

cor.test(test$mo_PFOS_cor, test$mo_PFHxS_cor, method = "pearson")
cor.test(test$mo_PFOS_cor, test$mo_PFHpS_i_cor, method = "pearson")

cor.test(test$ch_MEPA_total_i_cor_Y1, test$ch_PRPA_total_i_cor_Y1, method = "pearson")



phenols_pfas_cormat <- metadata[, c(phenols_vec_num, pfas_vec_num)]
colnames(phenols_pfas_cormat) <- colnames(phenols_pfas_cormat) %>%
  str_replace_all(c(
    "mo_" = "",
    "ch_" = "",
    "_total_i_cor_" = " ",
    "_total_cat_" = " ", 
    "_i_cor" = "",
    "_cor" = "",
    "MEPA" = "Methylparaben", 
    "ETPA" = "Ethylparaben", 
    "BUPA" = "Butylparaben", 
    "PRPA" = "Propylparaben", 
    "OXBE" = "Benzophenone 3", 
    "TRCS" = "Triclosan", 
    "BPA" = "Bisphenol A", 
    "BPS" = "Bisphenol S", 
    "PFHpS" = "PFHpS t2", 
    "PFHxS" = "PFHxS t2", 
    "PFOS" = "PFOs t2", 
    "PFNA" = "PFNA t2", 
    "PFOA" = "PFOA t2", 
    "PFDA" = "PFDA t2", 
    "PFUnDA" = "PFUnDA t2"
  ))

# méhode 1
heatmap <- heatmap_cor(cormat = phenols_pfas_cormat)

# méthode 2
phenols_pfas_cormat <- round(cor(phenols_pfas_cormat,
                    use = "pairwise.complete.obs",
                    method = "pearson"), 1)
corrplot(phenols_pfas_cormat, method = "circle", type = "lower", diag = TRUE)
corrplot(phenols_pfas_cormat, method = "number", type = "lower", add = TRUE)

ggsave(filename = "4_output/phenols pfas/heatmap_cor_pollutants_covar.pdf",
       plot = heatmap, 
       units = "cm",
       width = 15, 
       height = 20, 
       limitsize = FALSE)

### Correlations entre indices d'alpha diversité  ----
cor.test(bdd_taxa$ch_feces_rel_p1_Y1, bdd_taxa$ch_feces_rel_g3_Y1, method = "pearson")
cor.test(bdd_taxa$ch_feces_rel_p2_Y1, bdd_taxa$ch_feces_rel_g1_Y1, method = "pearson")
cor.test(bdd_taxa$ch_feces_rel_p3_Y1, bdd_taxa$ch_feces_rel_g2_Y1, method = "pearson")
cor.test(bdd_taxa$ch_feces_rel_p4_Y1, bdd_taxa$ch_feces_rel_g4_Y1, method = "pearson")

heatmap_alpha <- bdd_alpha %>% 
  select(ch_feces_SpecRich_5000_ASV_Y1, 
         ch_feces_Shannon_5000_ASV_Y1, 
         ch_feces_Faith_5000_ASV_Y1) %>%
  rename("Specific richness" = ch_feces_SpecRich_5000_ASV_Y1, 
         "Shannon diversity" = ch_feces_Shannon_5000_ASV_Y1, 
         "Faith diversity" = ch_feces_Faith_5000_ASV_Y1)

cormat <- round(cor(cormat, 
                    use = "pairwise.complete.obs", 
                    method = "spearman"), 3)
cormat <- reorder_cormat(cormat)                 # réordonner les coef de cor
upper_tri <- get_upper_tri(cormat)               # obtenir que le triangle sup
melted_cormat <- melt(upper_tri, na.rm = TRUE)   # passer en df long rapidement 

heatmap <-                                       # heatmap
  ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Spearman\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    size = 12,
    hjust = 1
  )) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value),
            color = "black",
            size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.4, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7,
    barheight = 1,
    title.position = "top",
    title.hjust = 0.5
  ))

### Figures Poster SF DohaD - forestplots ----
forestplot <- function(results_list, outcome_name) {
  results <- results_list %>%
    ggplot(aes(x = exposure, 
               y = estimate, 
               min = conf.low, 
               ymax = conf.high, 
               #color = interaction(exposure_window, term_2), 
               color = term_rec, 
               shape = p_value_shape)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.6) +
    labs(x = "Exposures", y = outcome_name) +
    theme_bw() +
    coord_flip()  +
    scale_shape_manual(values = c(19, 8),
                       name = "p-value") +
    guides(color = guide_legend(title = ""))+
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 17),
          legend.title = element_blank(), 
          legend.position = "bottom",
          legend.box = "vertical", 
          legend.justification = "center", 
          legend.spacing.y = unit(0, "cm"), 
          legend.spacing.x = unit(0, "cm"), 
          legend.box.margin = margin(0,0,0,0, "cm"), 
          legend.margin = margin(0,0,0,0, "cm"), 
          #panel.background = element_rect(fill = "transparent", colour = NA), 
          plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
          legend.background = element_rect(fill = "transparent", colour = NA))
  
  return(results)
}

#### Figure 1 ----
leg <- results_multi %>%
  filter(outcome == "Faith phylogenetic diversity") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Faith phylogenetic diversity") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.background = element_rect(fill = "transparent", colour = NA))
leg <- get_legend(leg) %>% as_ggplot()

y_axis <-
  results_multi %>%
  filter(outcome == "Specific richness") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Specific richness")+
  theme(legend.position = "none")

forestplot_alpha_confirmatory_1 <-
  results_multi %>%
  filter(outcome == "Specific richness") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Specific richness") +
  theme(legend.position = "none")

forestplot_alpha_confirmatory_1 <-
  results_multi %>%
  filter(outcome == "Specific richness") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Specific richness") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_confirmatory_2 <-
  results_multi %>%
  filter(outcome == "Shannon diversity") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_confirmatory_3 <- 
  results_multi %>%
  filter(outcome == "Faith phylogenetic diversity") %>%
  filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Faith phylogenetic diversity") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())+
  theme(legend.position = "none")


ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/y_axis.tiff", 
       plot = y_axis, 
       device = "tiff",
       units = "mm",
       width = 130, 
       height = 280,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/forestplot_alpha_1.tiff", 
       plot = forestplot_alpha_confirmatory_1, 
       device = "tiff",
       units = "mm",
       width = 90, 
       height = 280,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/forestplot_alpha_2.tiff", 
       plot = forestplot_alpha_confirmatory_2, 
       device = "tiff",
       units = "mm",
       width = 90, 
       height = 280,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/forestplot_alpha_3.tiff", 
       plot = forestplot_alpha_confirmatory_3, 
       device = "tiff",
       units = "mm",
       width = 90, 
       height = 280,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/leg.tiff", 
       plot = leg, 
       device = "tiff",
       units = "mm",
       width = 378, 
       height = 20,
       dpi = 300,
       limitsize = FALSE)


#### Figure 2 ----
results_multi_red <- results_multi %>%
  filter(exposure %in% c("Ethylparaben t3", 
                         "Propylparaben Y1", 
                         "Butylparaben M2", 
                         "Bisphenol A t3", 
                         "Bisphenol S M2", 
                         "Benzophenone 3 t2")) 
leg <- results_multi_red %>%
  filter(outcome == "Firmicutes") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Firmicutes") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.background = element_rect(fill = "transparent", colour = NA))
leg <- get_legend(leg) %>% as_ggplot()

y_axis_taxa <-
  results_multi_red %>%
  filter(outcome == "Firmicutes") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Firmicutes") +
  theme(legend.position = "none")

forestplot_taxa_confirmatory_1 <-
  results_multi_red %>%
  filter(outcome == "Firmicutes") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Firmicutes") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_taxa_confirmatory_2 <-
  results_multi_red %>%
  filter(outcome == "Actinobacteria") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Actinobacteria") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_taxa_confirmatory_3 <- 
  results_multi_red %>%
  filter(outcome == "Bacteroidetes") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Bacteroidetes") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_taxa_confirmatory_4 <- 
  results_multi_red %>%
  filter(outcome == "Proteobacteria") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Proteobacteria") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())+
  theme(legend.position = "none")


ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/y_axis_taxa.tiff", 
       plot = y_axis_taxa, 
       device = "tiff",
       units = "mm",
       width = 110, 
       height = 100,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/forestplot_taxa_confirmatory_1.tiff", 
       plot = forestplot_taxa_confirmatory_1, 
       device = "tiff",
       units = "mm",
       width = 75, 
       height = 100,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/forestplot_taxa_confirmatory_2.tiff", 
       plot = forestplot_taxa_confirmatory_2, 
       device = "tiff",
       units = "mm",
       width = 75, 
       height = 100,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/forestplot_taxa_confirmatory_3.tiff", 
       plot = forestplot_taxa_confirmatory_3, 
       device = "tiff",
       units = "mm",
       width = 75, 
       height = 100,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/forestplot_taxa_confirmatory_4.tiff", 
       plot = forestplot_taxa_confirmatory_4, 
       device = "tiff",
       units = "mm",
       width = 80, 
       height = 100,
       dpi = 300,
       limitsize = FALSE)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/leg.tiff", 
       plot = leg, 
       device = "tiff",
       units = "mm",
       width = 378, 
       height = 20,
       dpi = 300,
       limitsize = FALSE)

