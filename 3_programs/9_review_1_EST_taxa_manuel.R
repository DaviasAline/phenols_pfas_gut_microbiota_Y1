## Aline Davias
## Révisions ES&T
## 15/01/2024
## Analyses de la taxonomie 


# Packages loading ----
library(tidyverse)
library(phyloseq)
library(expss)
library(gtsummary)
library(Maaslin2)
library(labelled)
library(questionr)
library(sjlabelled)
library(openxlsx)
source("3_programs/4_functions_AD_gumme.R")
rm(comp_effectifs, heatmap_cor_pairwise, model_covar, model_multi, model_summary, model_univ_multi, 
   table_cor, table_cor_sg, test_sensi_sg)
library(corrplot)
library(see)
library(psych)
library(compositions)
library(writexl)
library(broom)

# Data reading ----
## Gut microbiota ----
taxa_table <- read_csv("0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")
input_data <- 
  read_labelled_csv(
    "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
  select(
    ident, 
    starts_with("ch_feces_rel_g")) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  column_to_rownames("ident")

var_label(input_data) <-  str_replace(
  var_label(input_data),                                    # set correct variable names                   
  "One year child feces relative abundance of ", "")
colnames(input_data) <- var_label(input_data)



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


# Functions coding ----
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

# Fonction pour le choix du nombre de décimales
custom_pvalue_fun <- function(x) {
  sapply(x, function(p) {
    if (is.na(p)) {
      return(NA) # Retourner NA si p est NA
    } else if (p < 0.001) {
      # Pour p < 0.001, utiliser la notation scientifique pour afficher toutes les décimales
      return(format(p, scientific = TRUE))
    } else if (p >= 0.001 & p < 0.01) {
      # Pour 0.001 <= p < 0.01, afficher avec 3 décimales
      return(sprintf("%.3f", p))
    } else {
      # Pour p >= 0.01, afficher avec 2 décimales
      return(sprintf("%.2f", p))
    }
  })
}


# Genera selection ----
## Description des genres à étudier en régression linéaire (n=46) 
input_data %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  View()

genera_linear <- input_data %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()

descrip_genera_linear <- tbl_merge(
  tbls = 
    list(
      tbl_1 = 
        input_data %>% 
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_2 = 
        input_data %>% 
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
        mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
        set_label(genera_linear) %>%
        tbl_summary(
          type = list(everything() ~ "categorical"))), 
  tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))


## Description des genres à exclure (n=146) 
input_data %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
  View()

genera_excluded <- input_data %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
  colnames()

descrip_genera_excluded <- tbl_merge(
  tbls = 
    list(
      tbl_1 = 
        input_data %>% 
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_2 = 
        input_data %>% 
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
        mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
        set_label(genera_excluded) %>%
        tbl_summary(
          type = list(everything() ~ "categorical"))), 
  tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))


# Data cleaning ----
## Selection des genres à analyser
input_data_raw <- input_data %>% 
  select(all_of(genera_linear)) %>%
  mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%        # remplacement des valeurs 0 par 1/5000
  rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
  rownames_to_column(var = "ident")

input_data_log <- input_data %>% 
  select(all_of(genera_linear)) %>%
  mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%        # remplacement des valeurs 0 par 1/5000
  mutate_all(~ log(.)) %>%                              # transformation logarithmique
  rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
  rownames_to_column(var = "ident")

input_data_clr <- input_data %>% 
  select(all_of(genera_linear)) %>%
  mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%     # remplacement des valeurs 0 par 1/5000
  mutate_all(~ log(.))
rowMeans <- rowMeans(input_data_clr)
input_data_clr <- sweep(input_data_clr, 1, rowMeans, FUN = "-")   # on soustrait chaque valeur par la log mean de la lign
input_data_clr <- input_data_clr %>%
  rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
  rownames_to_column(var = "ident")
rm(rowMeans)

genera_linear <- str_replace_all(genera_linear, "genus ", "")

# Data description ----
## Covariates ----
descrip_covar <- 
  input_metadata %>% 
  select(all_of(covariates)) %>% 
  tbl_summary(type = list(all_categorical() ~ "categorical"))

## Included genera (n=46) ----
### Tables ----
descrip_genera_linear

comp_transfo <- tbl_merge(
  tbls = list(
    input_data_raw %>% select(-ident) %>% tbl_summary(),
    input_data_log %>% select(-ident) %>% tbl_summary(),
    input_data_clr %>% select(-ident) %>% tbl_summary()), 
  tab_spanner = c("**Raw**", "**Log**", "**Clr**"))
comp_transfo

### Density plots ----
densityplot(data = input_data_raw, vars = genera_linear[1:23])
densityplot(data = input_data_raw, vars = genera_linear[24:46])
densityplot(data = input_data_log, vars = genera_linear[1:23])
densityplot(data = input_data_log, vars = genera_linear[24:46])
densityplot(data = input_data_clr, vars = genera_linear[1:23])
densityplot(data = input_data_clr, vars = genera_linear[24:46])


### Heatmap of correlation ----
cormat_genera <- input_data %>% 
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  rename_with(~gsub("genus ", "", .), everything())
cormat_genera <- round(cor(cormat_genera, 
                           use = "pairwise.complete.obs", 
                           method = "pearson"), 1)
heatmap_genera <-                                       # heatmap
  reshape2::melt(cormat_genera, na.rm = TRUE) %>% # passer en df long rapidement 
  ggplot(aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    #size = 12,
    hjust = 1
  )) +
  coord_fixed() +
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
heatmap_genera

# Calcul du nombre de tests ----



# Running linear regressions ----
## Rassemblement des données en 1 seul dataframe 
input_metadata <- input_metadata %>% rownames_to_column(var = "ident")
data_df_log <- left_join(input_data_log, input_metadata, by = "ident")
data_df_raw <- left_join(input_data_raw, input_metadata, by = "ident")
data_df_clr <- left_join(input_data_clr, input_metadata, by = "ident")
corres <- 
  taxa_table %>% 
  select(Phyla_corres = ch_feces_phylum_ASVbased_Y1, 
         Outcome = ch_feces_genus_ASVbased_Y1) %>%
  filter(Outcome %in% genera_linear) %>%
  distinct(Outcome, .keep_all = TRUE)

## Version log ----
tbls_by_outcome_log <- vector("list", length(genera_linear))                        # création liste pour stocker les tableaux par outcome
names(tbls_by_outcome_log) <- genera_linear

for (outcome in genera_linear) {
  tbls_for_outcome_log <- vector("list", length(pollutants))
  names(tbls_for_outcome_log) <- pollutants
  
  for (exposure in pollutants) {                                                # formula setting
    terms <- c(exposure, covariates)
    formula <- reformulate(terms, response = outcome)
    model <- lm(formula, data = data_df_log)
    
    tbl <-                                                                      # running linear regression
      tbl_regression(
        model, 
        include = exposure,
        estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun,
        exponentiate = FALSE) %>%
      bold_p() %>%
      bold_labels() %>%
      add_global_p(include = exposure, singular.ok = TRUE, keep = TRUE)
    
    tbls_for_outcome_log[[exposure]] <- tbl
  }
  tbls_by_outcome_log[[outcome]] <- tbls_for_outcome_log
}

### Tableau pour l'article ----
stacked_tbls_by_outcome_log <- vector("list", length(genera_linear))            # Création liste pour stocker les tableaux empilés par outcome
names(stacked_tbls_by_outcome_log) <- genera_linear

for (outcome in names(tbls_by_outcome_log)) {                                   # Récupérer les tableaux de régression pour cet outcome
  tbls_for_this_outcome_log <- tbls_by_outcome_log[[outcome]]
  stacked_tbl_log <- do.call(tbl_stack, list(tbls = tbls_for_this_outcome_log)) # Empiler les tableaux en un seul tableau
  stacked_tbls_by_outcome_log[[outcome]] <- stacked_tbl_log                     # Ajouter le tableau empilé à la liste des tableaux empilés
}

results_tbl_log <- tbl_merge(tbls = stacked_tbls_by_outcome_log,                # Fusionner les tableaux empilés en un seul tableau
                             tab_spanner = genera_linear)


### Tableau pour générer des figures ----
table_log <- tibble()                                                           # Initialisation d'un tibble vide pour stocker les résultats finaux
for (i in seq_along(tbls_by_outcome_log)) {                                     # Nom de l'outcome pour cette itération
  outcome_name <- names(tbls_by_outcome_log)[i]
  
  for (j in seq_along(tbls_by_outcome_log[[i]])) {                              # Itération sur chaque tbl_regression dans la liste courante
    exposure_name <- names(tbls_by_outcome_log[[i]])[j]                         # Nom de la variable d'exposition pour cette itération
    tbl_data <- tbls_by_outcome_log[[i]][[j]] %>%                               # Extraction des données du tableau tbl_regression
      as_tibble() %>%
      mutate(Outcome = outcome_name, Exposure = exposure_name)
    
    table_log <- bind_rows(table_log, tbl_data)                                 # Ajout des données extraites au tibble final
  }
}
rm(terms, formula, model, 
   tbl, tbl_data, tbls_for_this_outcome_log, 
   stacked_tbl_log, tbls_for_outcome_log,
   stacked_tbls_by_outcome_log,
   exposure, exposure_name, i, j, outcome, outcome_name)

table_log <-                                                                    # Ajout variable de la correspondance en phyla
  left_join(table_log, corres, by = "Outcome") %>% 
  select(Phyla_corres, everything())

table_log <- table_log %>%
  select(Phyla_corres, 
         Outcome, 
         Pollutants = Exposure, 
         Beta = "**Beta**", 
         "95% CI" = "**95% CI**",
         "p-value" = "**p-value**", 
         "Characteristic" = "**Characteristic**") %>%
  mutate(
    Phyla_corres = as.factor(Phyla_corres), 
    Phyla_corres = fct_relevel(Phyla_corres,
                               "Firmicutes", "Actinobacteria", 
                               "Bacteroidetes", "Proteobacteria", 
                               "Verrucomicrobia", "Candidatus_Saccharibacteria"),
    Time_window = case_when(grepl("t2", Pollutants) ~ "Mother, pregnancy trim. 2", 
                            grepl("t3", Pollutants) ~ "Mother, pregnancy trim. 3",
                            grepl("M2", Pollutants) ~ "Child, 2 months",
                            grepl("Y1", Pollutants) ~ "Child, 12 months", 
                            .default = "Mother, pregnancy trim. 2"),
    Pollutants = str_replace_all(Pollutants,
                                 c(
                                   "mo_" = "",
                                   "ch_" = "",
                                   "MEPA" = "Methylparaben",
                                   "ETPA" = "Ethylparaben",
                                   "PRPA" = "Propylparaben",
                                   "BUPA" = "Butylparaben",
                                   "BPA" = "Bisphenol A",
                                   "BPS" = "Bisphenol S",
                                   "OXBE" = "Benzophenone 3",
                                   "TRCS" = "Triclosan",
                                   "_total_i_cor_" = "", 
                                   "_total_cat_" = "",
                                   "t2_2" = "",
                                   "t3_2" = "",
                                   "M2_2" = "",
                                   "t2" = "",
                                   "t3" = "",
                                   "M2" = "",
                                   "Y1" = "", 
                                   "_ter" = "", 
                                   "_ln" = "", 
                                   "_i_cor" = "", 
                                   "_cor" = "",
                                   "_cat_2" = "",
                                   "_cat" = ""
                                 )),
    Pollutants_Time_window = case_when(Time_window == "Mother, pregnancy trim. 2" ~ paste(Pollutants, "trim.2", sep = " "), 
                                       Time_window == "Mother, pregnancy trim. 3" ~ paste(Pollutants, "trim.3", sep = " "), 
                                       Time_window == "Child, 2 months" ~ paste(Pollutants, "2 months", sep = " "), 
                                       Time_window == "Child, 12 months" ~ paste(Pollutants, "12 months", sep = " ")),
    Pollutants_Time_window = 
      fct_relevel(Pollutants_Time_window, 
                  "8_2diPAP trim.2", "6_2diPAP trim.2", "PFOSA trim.2", "PFBS trim.2",
                  "PFTrDa trim.2", "PFHpA trim.2", "PFHxPA trim.2", "PFDoDa trim.2",
                  "PFHpS trim.2", "PFUnDA trim.2", "PFDA trim.2", "PFOS trim.2",
                  "PFHxS trim.2", "PFNA trim.2", "PFOA trim.2", "Benzophenone 3 12 months",
                  "Benzophenone 3 2 months", "Benzophenone 3 trim.3", "Benzophenone 3 trim.2",
                  "Bisphenol S 12 months", "Bisphenol S 2 months", "Bisphenol S trim.3",
                  "Bisphenol S trim.2", "Bisphenol A 12 months", "Bisphenol A 2 months",
                  "Bisphenol A trim.3", "Bisphenol A trim.2", "Triclosan 12 months",
                  "Triclosan 2 months", "Triclosan trim.3", "Triclosan trim.2",
                  "Butylparaben 12 months", "Butylparaben 2 months", "Butylparaben trim.3",
                  "Butylparaben trim.2", "Propylparaben 12 months", "Propylparaben 2 months",
                  "Propylparaben trim.3", "Propylparaben trim.2", "Ethylparaben 12 months",
                  "Ethylparaben 2 months", "Ethylparaben trim.3", "Ethylparaben trim.2",
                  "Methylparaben 12 months", "Methylparaben 2 months", "Methylparaben trim.3",
                  "Methylparaben trim.2"),
    Pollutants_Time_window_rec = str_replace_all(Pollutants_Time_window, 
                                                 c("trim.2" = "t2", 
                                                   "trim.3" = "t3", 
                                                   "12 months" = "Y1",
                                                   "2 months" = "M2", 
                                                   "Bisphenol " = "BP")), 
    Codage = ifelse(Characteristic %in% c("<LOQ", ">LOQ", "<LOD", ">LOD", "LOD-LOQ"), Characteristic, NA), 
    Codage_rec = case_when(
      Pollutants_Time_window_rec == "Ethylparaben M2" & Codage == ">LOQ" ~ ">LOQ, compared to <LOQ",
      Codage == ">LOD" ~ ">LOD, compared to <LOD",
      Codage == ">LOQ" ~ ">LOQ, compared to <LOD",
      Codage == "LOD-LOQ" ~ "LOD-LOQ, compared to <LOD" ), 
    Codage_rec = fct_relevel(Codage_rec,
                             ">LOQ, compared to <LOQ",
                             ">LOD, compared to <LOD",
                             "LOD-LOQ, compared to <LOD",
                             ">LOQ, compared to <LOD"),
    `p-value` = gsub("__", "", `p-value`),
    `p-value` = as.numeric(`p-value`),
    `q-value` = `p-value`/(29*31), 
    p_value_shape = ifelse(`p-value`<0.05, "p-value<0.05", "p-value≥0.05"),
    q_value_shape = ifelse(`q-value`<0.05, "q-value<0.05", "q-value≥0.05"), 
    antimicrobial = ifelse(str_detect(Pollutants_Time_window_rec, "Methylparaben|Ethylparaben|Propylparaben|Butylparaben|Triclosan"), "Yes", "No"), 
    antimicrobial = fct_relevel(antimicrobial, "Yes", "No"), 
    sens_beta = ifelse(Beta < 0, "Beta<0", "Beta≥0"), 
    sens_beta = fct_relevel(sens_beta, "Beta≥0", "Beta<0"))  %>% 
  separate(col = "95% CI", into = c("lower_CI", "upper_CI"), sep = ",", remove = FALSE) %>%
  mutate(
    lower_CI = as.numeric(lower_CI),
    upper_CI = as.numeric(upper_CI)
  ) %>%
  select(
    Phyla_corres,
    Outcome, 
    Pollutants, 
    Time_window, 
    Pollutants_Time_window, Pollutants_Time_window_rec, 
    Codage, Codage_rec,
    Beta, sens_beta, 
    "95% CI", lower_CI, upper_CI, 
    "p-value", p_value_shape, 
    "q-value", q_value_shape, 
    antimicrobial)


#### Gérer la variable sens_beta pour les variables à 2 niveaux d'expo ---- 
codage_lod_df <- table_log %>%
  filter(Codage_rec %in%  c(">LOD, compared to <LOD", ">LOQ, compared to <LOQ") & !is.na(sens_beta)) %>%
  select(Outcome, Pollutants_Time_window_rec, sens_beta) %>%
  distinct(Outcome, Pollutants_Time_window_rec, .keep_all = TRUE)               # S'assurer d'avoir des lignes uniques pour la jointure

table_log <- table_log %>%
  left_join(codage_lod_df, 
            by = c("Outcome", "Pollutants_Time_window_rec"),                    # Joindre cette information avec le dataframe original pour remplir les NA
            suffix = c("", "_fill")) %>%                                        # en utilisant group_by et summarise pour gérer les relations multiples
  group_by(Outcome, Pollutants_Time_window_rec) %>%
  mutate(sens_beta_rec = if_else(is.na(sens_beta), first(sens_beta_fill, default = NA), sens_beta)) %>%
  ungroup() %>%
  select(Phyla_corres,
         Outcome, 
         Pollutants, 
         Time_window, 
         Pollutants_Time_window, Pollutants_Time_window_rec, 
         Codage, Codage_rec,
         Beta, sens_beta, sens_beta_rec,
         "95% CI", lower_CI, upper_CI, 
         "p-value", p_value_shape, 
         "q-value", q_value_shape, 
         antimicrobial)  # Supprimer la colonne temporaire après le remplissage
rm(codage_lod_df)

#### Gérer la variable sens_beta pour les variables à 3 niveaux d'expo ---- 
replacement_value <- function(x) {
  unique_values <- na.omit(unique(x))
  if(length(unique_values) == 1) {
    return(unique_values)
  } else {
    return("Non linear")
  }
}
 
table_log <- table_log %>%                               # Applying the logic to the data frame
  group_by(Outcome, Pollutants_Time_window_rec) %>%
  mutate(consistent_sens_beta_rec = replacement_value(sens_beta_rec)) %>%
  ungroup() %>%
  mutate(
    sens_beta_rec = ifelse(is.na(Codage) & is.na(sens_beta_rec), consistent_sens_beta_rec, sens_beta_rec), 
    sens_beta_rec = case_when(
      sens_beta_rec == "1" ~ "Beta≥0", 
      sens_beta_rec == "2" ~ "Beta<0", 
      .default = sens_beta_rec)) %>%
  select(-consistent_sens_beta_rec) # remove the helper column 

table_log %>%                                                                        # vérification que ça a fonctionné
  filter(Pollutants_Time_window_rec %in% c("Butylparaben t2", 
                                           "Butylparaben t3", 
                                           "Butylparaben Y1", 
                                           "BPS Y1", 
                                           "PFDoDa t2", 
                                           "PFHxPA t2", 
                                           "PFTrDa t2")) %>% 
  select(-Phyla_corres, -Pollutants, -Time_window,  - Pollutants_Time_window, -"q-value", -"q_value_shape", -"p_value_shape", -antimicrobial) %>%
  View()

## Version clr ----
# tbls_by_outcome_clr <- vector("list", length(genera_linear))                        # création liste pour stocker les tableaux par outcome
# names(tbls_by_outcome_clr) <- genera_linear
# 
# for (outcome in genera_linear) {
#   tbls_for_outcome_clr <- vector("list", length(pollutants))
#   names(tbls_for_outcome_clr) <- pollutants
#   
#   for (exposure in pollutants) {                                                # formula setting
#     terms <- c(exposure, covariates)
#     formula <- reformulate(terms, response = outcome)
#     model <- lm(formula, data = data_df_clr)
#     
#     tbl <-                                                                      # running linear regression
#       tbl_regression(
#         model, 
#         include = exposure,
#         estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
#         pvalue_fun = custom_pvalue_fun,
#         exponentiate = FALSE) %>%
#       bold_p() %>%
#       bold_labels() %>%
#       add_global_p(include = exposure, singular.ok = TRUE, keep = TRUE)
#     
#     tbls_for_outcome_clr[[exposure]] <- tbl
#   }
#   tbls_by_outcome_clr[[outcome]] <- tbls_for_outcome_clr
# }
# 
# ### Tableau pour l'article ----
# stacked_tbls_by_outcome_clr <- vector("list", length(genera_linear))                # Création liste pour stocker les tableaux empilés par outcome
# names(stacked_tbls_by_outcome_clr) <- genera_linear
# 
# for (outcome in names(tbls_by_outcome_clr)) {                                       # Récupérer les tableaux de régression pour cet outcome
#   tbls_for_this_outcome <- tbls_by_outcome_clr[[outcome]]
#   stacked_tbl_clr <- do.call(tbl_stack, list(tbls = tbls_for_this_outcome))         # Empiler les tableaux en un seul tableau
#   stacked_tbls_by_outcome_clr[[outcome]] <- stacked_tbl_clr                             # Ajouter le tableau empilé à la liste des tableaux empilés
# }
# 
# results_tbl_clr <- tbl_merge(tbls = stacked_tbls_by_outcome_clr,                    # Fusionner les tableaux empilés en un seul tableau
#                              tab_spanner = genera_linear)
# 
# 
# ### Tableau pour générer des figures ----
# table_clr <- tibble()                                                          # Initialisation d'un tibble vide pour stocker les résultats finaux
# for (i in seq_along(tbls_by_outcome_clr)) {                                         # Nom de l'outcome pour cette itération
#   outcome_name <- names(tbls_by_outcome_clr)[i]
#   
#   for (j in seq_along(tbls_by_outcome_clr[[i]])) {                                  # Itération sur chaque tbl_regression dans la liste courante
#     exposure_name <- names(tbls_by_outcome_clr[[i]])[j]                             # Nom de la variable d'exposition pour cette itération
#     tbl_data <- tbls_by_outcome_clr[[i]][[j]] %>%                                   # Extraction des données du tableau tbl_regression
#       as_tibble() %>%
#       mutate(Outcome = outcome_name, Exposure = exposure_name)
#     
#     table_clr <- bind_rows(table_clr, tbl_data)                               # Ajout des données extraites au tibble final
#   }
# }
# 
# table_clr <- left_join(table_clr, corres, by = "Outcome") %>% select(Phyla_corres, everything())
# 
# table_clr <- table_clr %>%
#   select(Outcome, 
#          Pollutants = Exposure, 
#          Beta = "**Beta**", 
#          "95% CI" = "**95% CI**",
#          "p-value" = "**p-value**", 
#          "Characteristic" = "**Characteristic**") %>%
#   mutate(
#     Phyla_corres = as.factor(Phyla_corres), 
#     Phyla_corres = fct_relevel(Phyla_corres,
#                                "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria", "Verrucomicrobia", "Candidatus_Saccharibacteria"),
#     Time_window = case_when(grepl("t2", Pollutants) ~ "Mother, pregnancy trim. 2", 
#                             grepl("t3", Pollutants) ~ "Mother, pregnancy trim. 3",
#                             grepl("M2", Pollutants) ~ "Child, 2 months",
#                             grepl("Y1", Pollutants) ~ "Child, 12 months", 
#                             .default = "Mother, pregnancy trim. 2"),
#     Pollutants = str_replace_all(Pollutants,
#                                  c(
#                                    "mo_" = "",
#                                    "ch_" = "",
#                                    "MEPA" = "Methylparaben",
#                                    "ETPA" = "Ethylparaben",
#                                    "PRPA" = "Propylparaben",
#                                    "BUPA" = "Butylparaben",
#                                    "BPA" = "Bisphenol A",
#                                    "BPS" = "Bisphenol S",
#                                    "OXBE" = "Benzophenone 3",
#                                    "TRCS" = "Triclosan",
#                                    "_total_i_cor_" = "", 
#                                    "_total_cat_" = "",
#                                    "t2_2" = "",
#                                    "t3_2" = "",
#                                    "M2_2" = "",
#                                    "t2" = "",
#                                    "t3" = "",
#                                    "M2" = "",
#                                    "Y1" = "", 
#                                    "_ter" = "", 
#                                    "_ln" = "", 
#                                    "_i_cor" = "", 
#                                    "_cor" = "",
#                                    "_cat_2" = "",
#                                    "_cat" = ""
#                                  )),
#     Pollutants_Time_window = case_when(Time_window == "Mother, pregnancy trim. 2" ~ paste(Pollutants, "trim.2", sep = " "), 
#                                        Time_window == "Mother, pregnancy trim. 3" ~ paste(Pollutants, "trim.3", sep = " "), 
#                                        Time_window == "Child, 2 months" ~ paste(Pollutants, "2 months", sep = " "), 
#                                        Time_window == "Child, 12 months" ~ paste(Pollutants, "12 months", sep = " ")),
#     Pollutants_Time_window = 
#       fct_relevel(Pollutants_Time_window, 
#                   "8_2diPAP trim.2", "6_2diPAP trim.2", "PFOSA trim.2", "PFBS trim.2",
#                   "PFTrDa trim.2", "PFHpA trim.2", "PFHxPA trim.2", "PFDoDa trim.2",
#                   "PFHpS trim.2", "PFUnDA trim.2", "PFDA trim.2", "PFOS trim.2",
#                   "PFHxS trim.2", "PFNA trim.2", "PFOA trim.2", "Benzophenone 3 12 months",
#                   "Benzophenone 3 2 months", "Benzophenone 3 trim.3", "Benzophenone 3 trim.2",
#                   "Bisphenol S 12 months", "Bisphenol S 2 months", "Bisphenol S trim.3",
#                   "Bisphenol S trim.2", "Bisphenol A 12 months", "Bisphenol A 2 months",
#                   "Bisphenol A trim.3", "Bisphenol A trim.2", "Triclosan 12 months",
#                   "Triclosan 2 months", "Triclosan trim.3", "Triclosan trim.2",
#                   "Butylparaben 12 months", "Butylparaben 2 months", "Butylparaben trim.3",
#                   "Butylparaben trim.2", "Propylparaben 12 months", "Propylparaben 2 months",
#                   "Propylparaben trim.3", "Propylparaben trim.2", "Ethylparaben 12 months",
#                   "Ethylparaben 2 months", "Ethylparaben trim.3", "Ethylparaben trim.2",
#                   "Methylparaben 12 months", "Methylparaben 2 months", "Methylparaben trim.3",
#                   "Methylparaben trim.2"),
#     Pollutants_Time_window_rec = str_replace_all(Pollutants_Time_window, 
#                                                  c("trim.2" = "t2", 
#                                                    "trim.3" = "t3", 
#                                                    "12 months" = "Y1",
#                                                    "2 months" = "M2", 
#                                                    "Bisphenol " = "BP")), 
#     Codage = ifelse(Characteristic %in% c("<LOQ", ">LOQ", "<LOD", ">LOD", "LOD-LOQ"), Characteristic, NA), 
#     Codage_rec = fct_recode(Codage,
#                             ">LOD, compared to <LOD" = ">LOD",
#                             ">LOQ, compared to <LOD" = ">LOQ",
#                             "LOD-LOQ, compared to <LOD" = "LOD-LOQ"), 
#     Codage_rec = fct_relevel(Codage_rec,
#                              ">LOD, compared to <LOD",
#                              "LOD-LOQ, compared to <LOD",
#                              ">LOQ, compared to <LOD"),
#     `p-value` = gsub("__", "", `p-value`),
#     `p-value_rec` = as.numeric(`p-value`),
#     `q-value` = `p-value_rec`/(29*31), 
#     p_value_shape = ifelse(`p-value_rec`<0.05, "p-value<0.05", "p-value≥0.05"),
#     q_value_shape = ifelse(`q-value`<0.05, "q-value<0.05", "q-value≥0.05"), 
#     antimicrobial = ifelse(str_detect(Pollutants_Time_window_rec, "Methylparaben|Ethylparaben|Propylparaben|Butylparaben|Triclosan"), "Yes", "No"), 
#     antimicrobial = fct_relevel(antimicrobial, "Yes", "No"), 
#     sens_beta = ifelse(Beta < 0, "Beta<0", "Beta≥0"), 
#     sens_beta = fct_relevel(sens_beta, "Beta≥0", "Beta<0"))  %>% 
#   select(
#     Phyla_corres,
#     Outcome, 
#     Pollutants, 
#     Time_window, 
#     Pollutants_Time_window, Pollutants_Time_window_rec, 
#     Codage, Codage_rec,
#     Beta, sens_beta, 
#     "95% CI", 
#     "p-value", "p-value_rec", p_value_shape, 
#     "q-value", q_value_shape, 
#     antimicrobial)
# 



# Figures ----
### Mahatan plot final ----
mahatan_plot <- table_log  %>%
  filter(is.na(Codage)) %>%
  mutate(
    Pollutants_Time_window_rec_2 = ifelse(
      Outcome == "Romboutsia" & Pollutants_Time_window_rec %in% c("BPS M2", "Benzophenone 3 t2"), 
      "BPS M2, Benzophenone 3 t2", 
      Pollutants_Time_window_rec)) %>%
  filter(!(Outcome == "Romboutsia" & Pollutants_Time_window_rec == "BPS M2")) %>%
  ggplot(aes(x = -log10(`p-value`), y = Outcome)) +
  geom_point(aes(shape = sens_beta_rec), size = 2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05/(29*31)), linetype = "dashed", color = "blue") +
  theme_lucid() +
  labs(x = "-log10(P-value)", y = "Genera", shape = "") +
  geom_text(aes(label = ifelse(`p-value` < 0.0035, as.character(Pollutants_Time_window_rec_2), "")), hjust = 0, vjust = -0.3, angle = 35, size = 3.5) +
  scale_shape_manual(values = c("Beta<0" = 15, "Beta≥0" = 17, "Non linear" = 18)) +# 15: carré plein, 17: triangle plein
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right")
mahatan_plot
ggsave("4_output/taxa manual/manhattan_plot.tiff", 
       mahatan_plot, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 20, 
       width = 40)

### Forestplot final ----
forest_plot <- table_log %>% 
  filter(`p-value`<0.0035) %>% 
  mutate(Beta = as.numeric(Beta), 
         Codage_rec = if_else(is.na(Codage_rec), "Continuous", Codage_rec),
         Codage_rec = fct_relevel(Codage_rec,
                                  "Continuous", ">LOD, compared to <LOD", "LOD-LOQ, compared to <LOD",
                                  ">LOQ, compared to <LOD")) %>%
  ggplot(aes(x = Outcome, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = Pollutants_Time_window_rec,
             shape = Codage_rec)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Genera", y = "") +
  theme_lucid() +
  coord_flip()  +
  # facet_wrap(vars(Phyla_corres), scales = "free_y", ncol = 1) +
  # scale_shape_manual(values = c(19, 8),
  #                    name = "p-value") +
  guides(color = guide_legend(title = "Exposure and time window"), 
         shape = guide_legend(title = "Exposure coding"))+
  theme(
    # axis.title = element_text(size = 7),
    #     axis.text = element_text(size = 6),
    #     legend.text = element_text(size = 7),
    #     legend.title = element_text(size = 7), 
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right", 
    # legend.spacing.y = unit(0, "cm"), 
    # legend.spacing.x = unit(0, "cm"), 
    # legend.box.margin = margin(0,0,0,0, "cm"), 
    # legend.margin = margin(0,0,0,0, "cm")
  ) 

forest_plot
ggsave("4_output/taxa manual/forest_plot.tiff", 
       forest_plot, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)

# Strongest associations ----
## p<0.0035 ----

table_log %>%
  filter(is.na(Codage_rec)) %>%
  filter(`p-value`<0.0035)%>%
  View()

table_log %>%
  filter((Outcome == "Escherichia_Shigella" & Pollutants_Time_window_rec == "Benzophenone 3 M2") |
           (Outcome == "Clostridium_XlVa" & Pollutants_Time_window_rec == "Benzophenone 3 t2")|
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "BPA Y1") | 
           (Outcome == "Anaerostipes" & Pollutants_Time_window_rec == "Propylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Ethylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Enterobacter" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Collinsella" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "Benzophenone 3 t2") |
           (Outcome == "Klebsiella" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Coprococcus" & Pollutants_Time_window_rec == "BPA Y1") |
           (Outcome == "Lactococcus" & Pollutants_Time_window_rec == "BPA M2") |
           (Outcome == "Anaerotruncus" & Pollutants_Time_window_rec == "Butylparaben Y1")) %>%
  View()


results_signi <- tbl_merge(
  tbls = list(
    tbls_by_outcome_log$Escherichia_Shigella$ch_OXBE_total_i_cor_M2_ln,
    tbls_by_outcome_log$Clostridium_XlVa$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_log$Lachnospiracea_incertae_sedis$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_log$Lachnospiracea_incertae_sedis$ch_BPA_total_i_cor_Y1_ln,
    tbls_by_outcome_log$Anaerostipes$mo_PRPA_total_i_cor_t2_ln,
    tbls_by_outcome_log$Enterococcus$mo_ETPA_total_i_cor_t2_ln,
    tbls_by_outcome_log$Enterococcus$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_log$Enterobacter$ch_ETPA_total_i_cor_Y1_ln,
    tbls_by_outcome_log$Collinsella$ch_BPS_total_cat_M2_2,
    tbls_by_outcome_log$Romboutsia$ch_BPS_total_cat_M2_2,
    tbls_by_outcome_log$Romboutsia$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_log$Klebsiella$ch_ETPA_total_i_cor_Y1_ln,
    tbls_by_outcome_log$Coprococcus$ch_BPA_total_i_cor_Y1_ln,
    tbls_by_outcome_log$Lactococcus$ch_BPA_total_i_cor_M2_ln,
    tbls_by_outcome_log$Anaerotruncus$ch_BUPA_total_cat_Y1), 
  tab_spanner = c("**Escherichia_Shigella**", 
                  "**Clostridium_XlVa**",
                  "**Lachnospiracea_incertae_sedis**",
                  "**Lachnospiracea_incertae_sedis**",
                  "**Anaerostipes**",
                  "**Enterococcus**",
                  "**Enterococcus**",
                  "**Enterobacter**",
                  "**Collinsella**",
                  "**Romboutsia**",
                  "**Romboutsia**",
                  "**Klebsiella**",
                  "**Coprococcus**",
                  "**Lactococcus**",
                  "**Anaerotruncus**"))

write_xlsx(
  x = taxa_table %>%
    select(-"ch_feces_ASV_ID_Y1", -"ch_feces_domain_ASVbased_Y1", -"ch_feces_TAX_ASVbased_Y1") %>%
    filter(ch_feces_genus_ASVbased_Y1 %in% c("Escherichia_Shigella", 
                                             "Clostridium_XlVa",
                                             "Lachnospiracea_incertae_sedis",
                                             "Anaerostipes",
                                             "Enterococcus",
                                             "Enterobacter",
                                             "Collinsella",
                                             "Romboutsia",
                                             "Klebsiella",
                                             "Coprococcus",
                                             "Lactococcus",
                                             "Anaerotruncus")) %>% 
    mutate_if(is.character, as.factor) %>%
    arrange( ch_feces_phylum_ASVbased_Y1, 
             ch_feces_class_ASVbased_Y1,
             ch_feces_order_ASVbased_Y1,  
             ch_feces_family_ASVbased_Y1 ) %>%
    distinct(), 
  path = "4_output/taxa manual/correspondnace taxa des principales asso.xlsx")


## p<0.001 ----
table_log %>%
  filter(is.na(Codage_rec)) %>%
  filter(`p-value`<0.001)%>%
  View()


table_log %>%
  filter((Outcome == "Escherichia_Shigella" & Pollutants_Time_window_rec == "Benzophenone 3 M2") |
           (Outcome == "Clostridium_XlVa" & Pollutants_Time_window_rec == "Benzophenone 3 t2")|
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "BPA Y1") | 
           (Outcome == "Anaerostipes" & Pollutants_Time_window_rec == "Propylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Ethylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Enterobacter" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Collinsella" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "Benzophenone 3 t2") |
           (Outcome == "Klebsiella" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Coprococcus" & Pollutants_Time_window_rec == "BPA Y1") |
           (Outcome == "Lactococcus" & Pollutants_Time_window_rec == "BPA M2") |
           (Outcome == "Anaerotruncus" & Pollutants_Time_window_rec == "Butylparaben Y1")) %>%
  View()


results_signi_0.001 <- tbl_merge(
  tbls = list(
    # tbls_by_outcome_log$Escherichia_Shigella$ch_OXBE_total_i_cor_M2_ln,
    # tbls_by_outcome_log$Clostridium_XlVa$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_log$Lachnospiracea_incertae_sedis$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_log$Lachnospiracea_incertae_sedis$ch_BPA_total_i_cor_Y1_ln,
    # tbls_by_outcome_log$Anaerostipes$mo_PRPA_total_i_cor_t2_ln,
    # tbls_by_outcome_log$Enterococcus$mo_ETPA_total_i_cor_t2_ln,
    tbls_by_outcome_log$Enterococcus$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_log$Enterobacter$ch_ETPA_total_i_cor_Y1_ln,
    tbls_by_outcome_log$Collinsella$ch_BPS_total_cat_M2_2,
    # tbls_by_outcome_log$Romboutsia$ch_BPS_total_cat_M2_2,
    # tbls_by_outcome_log$Romboutsia$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_log$Klebsiella$ch_ETPA_total_i_cor_Y1_ln,
    # tbls_by_outcome_log$Coprococcus$ch_BPA_total_i_cor_Y1_ln,
    tbls_by_outcome_log$Lactococcus$ch_BPA_total_i_cor_M2_ln
    # tbls_by_outcome_log$Anaerotruncus$ch_BUPA_total_cat_Y1
  ), 
  tab_spanner = c(
    # "**Escherichia_Shigella**", 
    # "**Clostridium_XlVa**",
    "**Lachnospiracea_incertae_sedis**",
    "**Lachnospiracea_incertae_sedis**",
    # "**Anaerostipes**",
    # "**Enterococcus**",
    "**Enterococcus**",
    "**Enterobacter**",
    "**Collinsella**",
    # "**Romboutsia**",
    # "**Romboutsia**",
    "**Klebsiella**",
    # "**Coprococcus**",
    "**Lactococcus**"
    # "**Anaerotruncus**"
  ))


