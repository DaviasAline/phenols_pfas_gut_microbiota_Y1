# Analyses BKMR exposition aux polluants et microbiote 1 an 
# A. Davias
# 23/02/2023

load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')
library(bkmr)
library(fields)
library(future)
library(future.apply)
library(patchwork)
library(writexl)
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)



# Nettoyage des données  ----

## variables outcomes ----
bdd_outcomes_alpha <- 
  bdd_alpha %>% 
  select(
    ident, 
    ch_feces_SpecRich_5000_ASV_Y1,   
    ch_feces_Shannon_5000_ASV_Y1,
    ch_feces_Faith_5000_ASV_Y1) %>%
  na.omit()

bdd_outcomes_taxa <- bdd_taxa %>%
  select(ident, 
         all_of(taxa_vec)) %>%
  na.omit

## variables d'ajustement ----
bdd_covariates_i <-          # on force les variables catégorielles en variables continues 
  metadata %>%
  select(
    ident, 
    all_of(covar_vec), 
    all_of(covar_vec_i), 
    all_of(covar_vec_num_i)) %>%
  select(
    ident,  
    ch_feces_RUN_Y1,              # covariable à 2 catégories : ok 0/1
    ch_feces_age_w_Y1_i,  
    po_delmod,                    # covariable à 2 catégories : ok 0/1
    ch_food_intro_Y1_3cat_i,      # covariable à transformer 
    ch_antibio_Y1_i,         
    mo_par,  
    mo_pets_i,                    # covariable à 2 catégories : ok 0/1
    ch_sex,                       # covariable à 2 catégories : ok 0/1
    mo_tob_gr_anyt_yn_n2_i,       # covariable à 2 catégories : ok 0/1
    Mo_ETS_anyT_yn1_opt_i,        # covariable à 2 catégories : ok 0/1
    ch_ETS_12m_opt36m,            # covariable à 2 catégories : ok 0/1
    mo_interpreg_3cat,            # covariable à transformer 
    mo_dipl_3cat_i,               # covariable à transformer
    po_w_kg,
    po_he_i,
    ch_w_Y1_i,
    ch_he_Y1_i,
    po_gd,
    mo_age, 
    mo_bmi_bepr_i,
    bf_duration_till48w_i) %>%
  mutate(
    ch_feces_RUN_Y1 = as.numeric(fct_recode(ch_feces_RUN_Y1,
                                            "0" = "R2",
                                            "1" = "R3")), 
    po_delmod = as.numeric(fct_recode(po_delmod,
                                      "0" = "C-section",
                                      "1" = "Vaginal delivery")),
    mo_pets_i = as.numeric(fct_recode(mo_pets_i,
                                      "0" = "No",
                                      "1" = "One or more")),
    ch_sex = as.numeric(fct_recode(ch_sex,
                                   "0" = "Female",
                                   "1" = "Male")),
    mo_tob_gr_anyt_yn_n2_i = as.numeric(fct_recode(mo_tob_gr_anyt_yn_n2_i,
                                                   "0" = "No",
                                                   "1" = "Yes")),
    Mo_ETS_anyT_yn1_opt_i = as.numeric(fct_recode(Mo_ETS_anyT_yn1_opt_i,
                                                  "0" = "No",
                                                  "1" = "Yes")), 
    ch_ETS_12m_opt36m = as.numeric(fct_recode(ch_ETS_12m_opt36m,
                                              "0" = "No",
                                              "1" = "Yes")), 
    
    mo_interpreg_3cat_2years_and_more = as.numeric(fct_recode(mo_interpreg_3cat,
                                                              "0" = "Under 2 years",
                                                              "1" = "2 years and more",
                                                              "0" = "Primiparous")),
    mo_interpreg_3cat_primiparous = as.numeric(fct_recode(mo_interpreg_3cat,
                                                          "0" = "Under 2 years",
                                                          "0" = "2 years and more",
                                                          "1" = "Primiparous")), 
    
    ch_food_intro_Y1_3cat_i_6_12m = as.numeric(fct_recode(ch_food_intro_Y1_3cat_i,
                                                          "0" = "Between 0 and 6 months old",
                                                          "1" = "Between 6 and 12 months old",
                                                          "0" = "Not introduced at 12 months old")),
    ch_food_intro_Y1_3cat_i_not_intro = as.numeric(fct_recode(ch_food_intro_Y1_3cat_i,
                                                              "0" = "Between 0 and 6 months old",
                                                              "0" = "Between 6 and 12 months old",
                                                              "1" = "Not introduced at 12 months old")),
    mo_dipl_3cat_i_3_4y = as.numeric(fct_recode(mo_dipl_3cat_i,
                                                "0" = "2years or less after graduation",
                                                "1" = "3-4years after graduation",
                                                "0" = ">=5years after graduation")), 
    mo_dipl_3cat_i_5y = as.numeric(fct_recode(mo_dipl_3cat_i,
                                              "0" = "2years or less after graduation",
                                              "0" = "3-4years after graduation",
                                              "1" = ">=5years after graduation"))) %>%
  select(-c("mo_interpreg_3cat", "mo_dipl_3cat_i", "ch_food_intro_Y1_3cat_i"))
keep_covar <- bdd_covariates_i[, 2:25] %>% colnames()

## variables d'exposition ----
# possibilité 1 : on ne prend que les covariables qui sont numériques (version Nicolas)

bdd_exposure_1 <- metadata %>%          
  select(
    ident,
    all_of(phenols_vec_2), 
    all_of(pfas_vec))%>%           # variables d'exposition (déjà logtransformé)
  select(
    ident, 
    where(is.numeric))%>%               # conserver que les expositions codées en numérique
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
keep_expo_1 <- colnames(bdd_exposure_1[, 2:30])
bdd_exposure_1[, keep_expo_1] <- lapply(bdd_exposure_1[, keep_expo_1], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_exposure_1) <- c("ident", keep_expo_1)


# Préparation des matrices ----
## possibilité 1 : fenetres d'exposition confondues, sans variables catégorielles 
bdd_bkmr_alpha_1 <- 
  list(bdd_outcomes_alpha, bdd_covariates_i, bdd_exposure_1) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_1 <- bdd_bkmr_alpha_1 %>% select(all_of(keep_expo_1)) %>% as.matrix()
covariates_1 <- bdd_bkmr_alpha_1 %>% select(all_of(keep_covar)) %>% as.matrix()

outcome_specrich <- bdd_bkmr_alpha_2 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon <- bdd_bkmr_alpha_2 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()
outcome_faith <- bdd_bkmr_alpha_2 %>% select(ch_feces_Faith_5000_ASV_Y1) %>% as.matrix()

outcome_vec <- c("outcome_specrich", "outcome_shannon", "outcome_faith")



# Création des fonctions ----
bkmr_hierar_bis <- function(outcome) {
  set.seed(111)
  
  results_hierar_num_window_bis <- kmbayes(    # modèle sans les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_num, 
    X = covariates, 
    iter = 50000,          # mettre 50 000
    groups = c(1,2,3,4,    # methylparaben continue
               1,2,4,      # ethyparaben continue sauf M2
               1,2,4,      # propylparaben continue sauf M2
               # butylparaben (cat)
               1,2,3,4,    # bisphenol A continue
               # bisphenol S (cat)
               1,2,3,4,     # benzophenone 3 continue
               1,2,3,4,    # Triclosan continue
               5,5,5,5,5,5,5),   # 7 PFAS t2
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  
  results_hierar_all_window_bis <- kmbayes(    # modèle avec les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_all, 
    X = covariates, 
    iter = 50000,          # mettre 50 000
    groups = c(1,2,3,4,     # methylparaben 
               1,2,3,4,       # ethyparaben 
               1,2,3,4,       # propylparaben 
               1,2,3,4,       # butylparaben 
               1,2,3,4,       # bisphenol A 
               1,2,3,4,       # bisphenol S 
               1,2,3,4,       # benzophenone 3 
               1,2,3,4,      # Triclosan 
               5,5,5,5,5,5,5,5,5,5,5,5,5,5,5),
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  
  results_hierar_t2_num <- kmbayes(    # modèle sans les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_t2_num, 
    X = covariates, 
    iter = 50000,          # mettre 50 000
    groups = c(1,1,1,1,1,1,     # phenols
               2,2,2,2,2,2,2),  # pfas 
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  
  results_hierar_t2_all <- kmbayes(    # modèle avec les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_t2_all, 
    X = covariates, 
    iter = 50000,          # mettre 50 000
    groups = c(1,1,1,1,1,1,1,1,     # phenols
               2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),  # pfas 
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  
  
  results <- list(results_hierar_num_window_bis, 
                  results_hierar_all_window_bis, 
                  results_hierar_t2_num, 
                  results_hierar_t2_all)
  
  return(results)
}




TracePlot_group <- function(model_specrich, model_shannon, model_faith, titre){
  par(mfrow=c(3,3))
  TracePlot(fit = model_specrich, par = "beta", sel = TRUE) 
  TracePlot(fit = model_specrich, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_specrich, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_shannon, par = "beta", sel = TRUE) 
  TracePlot(fit = model_shannon, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_shannon, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_faith, par = "beta", sel = TRUE) 
  TracePlot(fit = model_faith, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_faith, par = "r", comp = 1, sel = TRUE) 
  mtext(titre, outer=TRUE, line = -1.5,font=2, cex=1, padj = 0)
}


pip_results_hierar <- function(bkmr_specrich, bkmr_shannon, bkmr_faith) {
  bkmr_pip_spechrich <- 
    ExtractPIPs(bkmr_specrich) %>% 
    as.data.frame() %>% 
    rename(groupPIP_specrich = groupPIP, 
           condPIP_specrich = condPIP) 
  
  bkmr_pip_shannon <- 
    ExtractPIPs(bkmr_shannon) %>% 
    as.data.frame() %>% 
    rename(groupPIP_shannon = groupPIP, 
           condPIP_shannon = condPIP) 
  
  bkmr_pip_faith <- 
    ExtractPIPs(bkmr_faith) %>% 
    as.data.frame() %>% 
    rename(groupPIP_faith = groupPIP, 
           condPIP_faith = condPIP) 
  
  results <- left_join(bkmr_pip_spechrich, bkmr_pip_shannon, bkmr_pip_faith, by = c("variable", "group"))
  results <- left_join(results, bkmr_pip_faith, by = c("variable", "group"))
  return(results)
  
}


risks_overall <- function(fit, y, Z) {
  results <- OverallRiskSummaries(
    fit = fit, 
    y = y,
    Z = Z,
    X = covariates,
    qs = seq(0.25, 0.95, by = 0.10),
    q.fixed = 0.25,
    method = "approx")
  return(results)
}     

risks_singvar <- function(fit, y, Z) {
  results <- SingVarRiskSummaries(
    fit = fit, 
    y = y,
    Z = Z,   
    X = covariates,
    qs.diff = c(0.25, 0.75),
    q.fixed = c(0.25, 0.50, 0.75),
    method = "approx")
  return(results)
}


plot_risks.overall <- function(risks.overall_sperich, 
                               risks.overall_shannon, 
                               risks.overall_faith
                               #, plot_title
                               ) {
  
  plot_risks.overall_sperich <-
    ggplot(risks.overall_sperich,
           aes(
             quantile,
             est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd
           )) +
    geom_pointrange() +
    labs(y = "Specific richness") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() +
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA))
  
  plot_risks.overall_shannon <-
    ggplot(risks.overall_shannon,
           aes(
             quantile,
             est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd
           )) +
    geom_pointrange() + 
    labs(y = "Shannon diversity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() + 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA))
  
  plot_risks.overall_faith <-
    ggplot(risks.overall_faith,
           aes(
             quantile,
             est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd
           )) +
    geom_pointrange() +
    labs(y = "Faith phylogenetic diversity")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() + 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA))
  
  plot_risks.overall <- 
    plot_risks.overall_sperich + 
    plot_risks.overall_shannon + 
    plot_risks.overall_faith + 
    plot_layout(ncol = 3) 
  #+ plot_annotation(title = plot_title)
  
  return(plot_risks.overall)
}



plot_risks.overall_sperich <- function(risks.overall_sperich, titre){
  ggplot(risks.overall_sperich,
         aes(
           quantile,
           est,
           ymin = est - 1.96 * sd,
           ymax = est + 1.96 * sd
         )) +
    geom_pointrange() +
    labs(y = "Specific richness") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() + 
    ggtitle(titre) + 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA))
}



plot_risks.overall_shannon <-function(risks.overall_shannon)
{   ggplot(risks.overall_shannon,
           aes(
             quantile,
             est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd
           )) +
    geom_pointrange() + 
    labs(y = "Shannon diversity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+ 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA))
  }

plot_risks.overall_faith <- function(risks.overall_faith)
{    ggplot(risks.overall_faith,
            aes(
              quantile,
              est,
              ymin = est - 1.96 * sd,
              ymax = est + 1.96 * sd
            )) +
    geom_pointrange() +
    labs(y = "Faith phylogenetic diversity")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+ 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA))
  } 



plot_risks.singvar <- function(risks.singvar_sperich, 
                               risks.singvar_shannon, 
                               risks.singvar_faith
                               #,plot_title
                               ) {
  
  
  plot_risks.singvar_sperich <- risks.singvar_sperich %>%
    mutate(
      variable = str_replace_all(variable,
                                 c("mo_" = "",
                                   "ch_" = "",
                                   "_total_cont_log_" = " ",
                                   "_cont_log" = " t2",
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
                                   "BPS" = "Bisphenol S",
                                   "OXBE" = "Benzophenone 3", 
                                   "TRCS" = "Triclosan")), 
      
      variable = fct_recode(variable,
                            "PFDA t2" = "PFDA ",
                            "PFHpS t2" = "PFHpS ",
                            "PFHxS t2" = "PFHxS",
                            "PFNA t2" = "PFNA",
                            "PFOA t2" = "PFOA",
                            "PFOS t2" = "PFOS",
                            "PFUnDA t2" = "PFUnDA "), 
      variable = fct_relevel(variable,
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
                             "Methylparaben Y1", "Methylparaben M2", "Methylparaben t3", "Methylparaben t2"))%>% 
    ggplot(
      aes(
        variable,
        est,
        ymin = est - 1.96 * sd,
        ymax = est + 1.96 * sd,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75), 
                    size = 0.15) +
    geom_hline(yintercept = 0, linetype="dashed") +
    labs(y = "Specific richness", 
         x ="Exposure") + 
    coord_flip()+
    theme_bw() +
    theme(legend.position = "none")
  
  plot_risks.singvar_shannon <- risks.singvar_shannon %>%
    mutate(
      variable = str_replace_all(variable,
                                 c("mo_" = "",
                                   "ch_" = "",
                                   "_total_cont_log_" = " ",
                                   "_cont_log" = " t2",
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
                                   "BPS" = "Bisphenol S",
                                   "OXBE" = "Benzophenone 3", 
                                   "TRCS" = "Triclosan")), 
      
      variable = fct_recode(variable,
                            "PFDA t2" = "PFDA ",
                            "PFHpS t2" = "PFHpS ",
                            "PFHxS t2" = "PFHxS",
                            "PFNA t2" = "PFNA",
                            "PFOA t2" = "PFOA",
                            "PFOS t2" = "PFOS",
                            "PFUnDA t2" = "PFUnDA "),
      variable = fct_relevel(variable,
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
                             "Methylparaben Y1", "Methylparaben M2", "Methylparaben t3", "Methylparaben t2"))%>% 
    ggplot(
      aes(
        variable,
        est,
        ymin = est - 1.96 * sd,
        ymax = est + 1.96 * sd,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75), 
                    size = 0.15) +
    geom_hline(yintercept = 0, linetype="dashed") + 
    labs(y = "Shannon diversity", 
         x ="Exposure")+
    coord_flip() +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
    theme(legend.position = "none")
  
  plot_risks.singvar_faith <- risks.singvar_faith %>%
    mutate(
      variable = str_replace_all(variable,
                                 c("mo_" = "",
                                   "ch_" = "",
                                   "_total_cont_log_" = " ",
                                   "_cont_log" = " t2",
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
                                   "BPS" = "Bisphenol S",
                                   "OXBE" = "Benzophenone 3", 
                                   "TRCS" = "Triclosan")), 
      
      variable = fct_recode(variable,
                            "PFDA t2" = "PFDA ",
                            "PFHpS t2" = "PFHpS ",
                            "PFHxS t2" = "PFHxS",
                            "PFNA t2" = "PFNA",
                            "PFOA t2" = "PFOA",
                            "PFOS t2" = "PFOS",
                            "PFUnDA t2" = "PFUnDA "), 
      variable = fct_relevel(variable,
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
                             "Methylparaben Y1", "Methylparaben M2", "Methylparaben t3", "Methylparaben t2"))%>% 
    ggplot(
      aes(
        variable,
        est,
        ymin = est - 1.96 * sd,
        ymax = est + 1.96 * sd,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75), 
                    size = 0.15) +
    geom_hline(yintercept = 0, linetype="dashed") +
    labs(y = "Faith phylogenetic diversity", 
         x ="Exposure") + 
    coord_flip()+
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  
  
  plot_risks.singvar <- 
    plot_risks.singvar_sperich + 
    plot_risks.singvar_shannon + 
    plot_risks.singvar_faith + 
    plot_layout(ncol = 3)
  #+ plot_annotation(title = plot_title)
  
  return(plot_risks.singvar)
  
}

# PARTIE 1: GENERER LES DONNEES ----
# Nettoyage des données ----
## selon méthode Claire 
## possibilité 3 : on prend toutes les expositions, séparée par fenetre d'exposition ----
bdd_exposure_num <- metadata %>%          
  select(
    ident,
    all_of(phenols_vec_2), 
    all_of(pfas_vec))%>%           # variables d'exposition (déjà logtransformé)
  select(
    ident, 
    where(is.numeric))                    # conserver que les lignes avec toutes les observations completes sans données manquantes 
keep_expo_1 <- colnames(bdd_exposure_1[, 2:30])
bdd_exposure_1[, keep_expo_1] <- lapply(bdd_exposure_1[, keep_expo_1], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_exposure_1) <- c("ident", keep_expo_1)



bdd_exposure_cat <- metadata %>%
  select(ident, 
         contains("ETPA"),
         contains("PRPA"),
         contains("BUPA"), 
         contains("BPS"), 
         contains("PFDoDA"),
         contains("PFHxPA"), 
         contains("PFHpA"),
         contains("PFTrDA"), 
         contains("PFBS"), 
         contains("PFOSA"), 
         contains("diPAP"), 
         phenols_vec_2, 
         pfas_vec) %>%
  rename(
    mo_8_2diPAP_LOD = "_8_2diPAP_LOD", 
    mo_6_2diPAP_LOD = "_6_2diPAP_LOD"
  ) %>%
  mutate(
    ch_ETPA_total_cat_M2 = as.numeric(fct_recode(ch_ETPA_total_cat_M2,"1" = "<LOD",  "2" = "LOD-LOQ", "3" = ">LOQ")), 
    ch_ETPA_total_string_M2 = as.numeric(ch_ETPA_total_string_M2),
    ch_ETPA_total_LOD_M2 = as.numeric(ch_ETPA_total_LOD_M2),
    ch_ETPA_total_cont_M2 = ifelse(ch_ETPA_total_cat_M2 > 1, ch_ETPA_total_string_M2, ch_ETPA_total_LOD_M2/sqrt(2)), 
    ch_ETPA_total_cont_M2 = as.numeric(ch_ETPA_total_cont_M2), 
    ch_ETPA_total_cont_log_M2 = scale(log(ch_ETPA_total_cont_M2)), 
    
    ch_PRPA_total_cat_M2 = as.numeric(fct_recode(ch_PRPA_total_cat_M2,"1" = "<LOD",  "2" = "LOD-LOQ", "3" = ">LOQ")), 
    ch_PRPA_total_string_M2 = as.numeric(ch_PRPA_total_string_M2),
    ch_PRPA_total_LOD_M2 = as.numeric(ch_PRPA_total_LOD_M2),
    ch_PRPA_total_cont_M2 = ifelse(ch_PRPA_total_cat_M2 > 1, ch_PRPA_total_string_M2, ch_PRPA_total_LOD_M2/sqrt(2)), 
    ch_PRPA_total_cont_M2 = as.numeric(ch_PRPA_total_cont_M2), 
    ch_PRPA_total_cont_log_M2 = scale(log(ch_PRPA_total_cont_M2)),
    
    mo_BUPA_total_cat_t2 = as.numeric(mo_BUPA_total_cat_t2),
    mo_BUPA_total_cont_t2 = ifelse(mo_BUPA_total_cat_t2 > 1, mo_BUPA_total_string_t2, mo_BUPA_LOD/sqrt(2)), 
    mo_BUPA_total_cont_t2 = as.numeric(mo_BUPA_total_cont_t2), 
    mo_BUPA_total_cont_log_t2 = scale(log(mo_BUPA_total_cont_t2)),
    
    mo_BUPA_total_cat_t3 = as.numeric(mo_BUPA_total_cat_t3), 
    mo_BUPA_total_cont_t3 = ifelse(mo_BUPA_total_cat_t3 > 1, mo_BUPA_total_string_t3, mo_BUPA_LOD/sqrt(2)), 
    mo_BUPA_total_cont_t3 = as.numeric(mo_BUPA_total_cont_t3), 
    mo_BUPA_total_cont_log_t3 = scale(log(mo_BUPA_total_cont_t3)), 
    
    ch_BUPA_total_cat_M2 = as.numeric(fct_recode(ch_BUPA_total_cat_M2,"1" = "<LOD",  "2" = "LOD-LOQ", "3" = ">LOQ")), 
    ch_BUPA_total_string_M2 = as.numeric(ch_BUPA_total_string_M2),
    ch_BUPA_total_LOD_M2 = as.numeric(ch_BUPA_total_LOD_M2),
    ch_BUPA_total_cont_M2 = ifelse(ch_BUPA_total_cat_M2 > 1, ch_BUPA_total_string_M2, ch_BUPA_total_LOD_M2/sqrt(2)), 
    ch_BUPA_total_cont_M2 = as.numeric(ch_BUPA_total_cont_M2), 
    ch_BUPA_total_cont_log_M2 = scale(log(ch_BUPA_total_cont_M2)), 
    
    ch_BUPA_total_cat_Y1 = as.numeric(ch_BUPA_total_cat_Y1), 
    ch_BUPA_total_LOD_Y1 = as.numeric(ch_BUPA_total_LOD_Y1),
    ch_BUPA_total_cont_Y1 = ifelse(ch_BUPA_total_cat_Y1 > 1, ch_BUPA_total_string_Y1, ch_BUPA_total_LOD_Y1/sqrt(2)), 
    ch_BUPA_total_cont_Y1 = as.numeric(ch_BUPA_total_cont_Y1), 
    ch_BUPA_total_cont_log_Y1 = scale(log(ch_BUPA_total_cont_Y1)), 
    
    mo_BPS_total_cat_t2 = as.numeric(fct_recode(mo_BPS_total_cat_t2,"1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")), 
    mo_BPS_total_string_t2 = as.numeric(mo_BPS_total_string_t2),
    mo_BPS_LOD = as.numeric(mo_BPS_LOD),
    mo_BPS_total_cont_t2 = ifelse(mo_BPS_total_cat_t2 > 1, mo_BPS_total_string_t2, mo_BPS_LOD/sqrt(2)), 
    mo_BPS_total_cont_t2 = as.numeric(mo_BPS_total_cont_t2), 
    mo_BPS_total_cont_log_t2 = scale(log(mo_BPS_total_cont_t2)), 
    
    mo_BPS_total_cat_t3 = as.numeric(fct_recode(mo_BPS_total_cat_t3,"1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")), 
    mo_BPS_total_string_t3 = as.numeric(mo_BPS_total_string_t3),
    mo_BPS_LOD= as.numeric(mo_BPS_LOD),
    mo_BPS_total_cont_t3 = ifelse(mo_BPS_total_cat_t3 > 1, mo_BPS_total_string_t3, mo_BPS_LOD/sqrt(2)), 
    mo_BPS_total_cont_t3 = as.numeric(mo_BPS_total_cont_t3), 
    mo_BPS_total_cont_log_t3 = scale(log(mo_BPS_total_cont_t3)), 
    
    ch_BPS_total_cat_M2 = as.numeric(fct_recode(ch_BPS_total_cat_M2,"1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")), 
    ch_BPS_total_string_M2 = as.numeric(ch_BPS_total_string_M2),
    ch_BPS_total_LOD_M2= as.numeric(ch_BPS_total_LOD_M2),
    ch_BPS_total_cont_M2 = ifelse(ch_BPS_total_cat_M2 > 1, ch_BPS_total_string_M2, ch_BPS_total_LOD_M2/sqrt(2)), 
    ch_BPS_total_cont_M2 = as.numeric(ch_BPS_total_cont_M2), 
    ch_BPS_total_cont_log_M2 = scale(log(ch_BPS_total_cont_M2)), 
    
    ch_BPS_total_cat_Y1 = as.numeric(ch_BPS_total_cat_Y1), 
    ch_BPS_total_LOD_Y1= as.numeric(ch_BPS_total_LOD_Y1),
    ch_BPS_total_cont_Y1 = ifelse(ch_BPS_total_cat_Y1 > 1, ch_BPS_total_string_Y1, ch_BPS_total_LOD_Y1/sqrt(2)), 
    ch_BPS_total_cont_Y1 = as.numeric(ch_BPS_total_cont_Y1), 
    ch_BPS_total_cont_log_Y1 = scale(log(ch_BPS_total_cont_Y1)), 
    
    mo_PFDoDA = as.numeric(mo_PFDoDA), 
    mo_PFDoDa_cat = as.numeric(mo_PFDoDa_cat),
    mo_PFDoDa_cont = ifelse(mo_PFDoDa_cat > 1, mo_PFDoDA, PFDoDA_LOD/sqrt(2)), 
    mo_PFDoDa_cont = as.numeric(mo_PFDoDa_cont), 
    mo_PFDoDa_cont_log = scale(log(mo_PFDoDa_cont)), 
    
    mo_PFHxPA = as.numeric(mo_PFHxPA), 
    mo_PFHxPA_cat  = as.numeric(mo_PFHxPA_cat),
    mo_PFHxPA_cont = ifelse(mo_PFHxPA_cat > 1, mo_PFHxPA, PFHxPA_LOD/sqrt(2)), 
    mo_PFHxPA_cont = as.numeric(mo_PFHxPA_cont), 
    mo_PFHxPA_cont_log = scale(log(mo_PFHxPA_cont)), 
    
    mo_PFHpA = as.numeric(mo_PFHpA), 
    mo_PFHpA_cat  = as.numeric(fct_recode(mo_PFHpA_cat, "1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")),
    mo_PFHpA_cont = ifelse(mo_PFHpA_cat > 1, mo_PFHpA, PFHpA_LOD/sqrt(2)), 
    mo_PFHpA_cont = as.numeric(mo_PFHpA_cont), 
    mo_PFHpA_cont_log = scale(log(mo_PFHpA_cont)), 
    
    mo_PFTrDA =as.numeric(mo_PFTrDA),
    mo_PFTrDa_cat = as.numeric(mo_PFTrDa_cat), 
    mo_PFTrDa_cont = ifelse(mo_PFTrDa_cat > 1, mo_PFTrDA, PFTrDA_LOD/sqrt(2)), 
    mo_PFTrDa_cont = as.numeric(mo_PFTrDa_cont), 
    mo_PFTrDa_cont_log = scale(log(mo_PFTrDa_cont)), 
    
    mo_PFBS = as.numeric(mo_PFBS), 
    mo_PFBS_cat  = as.numeric(fct_recode(mo_PFBS_cat, "1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")),
    mo_PFBS_cont = ifelse(mo_PFBS_cat > 1, mo_PFBS, PFBS_LOD/sqrt(2)), 
    mo_PFBS_cont = as.numeric(mo_PFBS_cont), 
    mo_PFBS_cont_log = scale(log(mo_PFBS_cont)), 
    
    mo_PFOSA = as.numeric(mo_PFOSA), 
    mo_PFOSA_cat  = as.numeric(fct_recode(mo_PFOSA_cat, "1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")),
    mo_PFOSA_cont = ifelse(mo_PFOSA_cat > 1, mo_PFOSA, PFOSA_LOD/sqrt(2)), 
    mo_PFOSA_cont = as.numeric(mo_PFOSA_cont), 
    mo_PFOSA_cont_log = scale(log(mo_PFOSA_cont)), 
    
    mo_6_2diPAP_cat  = as.numeric(fct_recode(mo_6_2diPAP_cat, "1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")),
    mo__6_2diPAP = as.numeric(mo__6_2diPAP), 
    mo_6_2diPAP_LOD = as.numeric(mo_6_2diPAP_LOD),
    mo_6_2diPAP_cont = ifelse(mo_6_2diPAP_cat > 1, mo__6_2diPAP, mo_6_2diPAP_LOD/sqrt(2)), 
    mo_6_2diPAP_cont = as.numeric(mo_6_2diPAP_cont), 
    mo_6_2diPAP_cont_log = scale(log(mo_6_2diPAP_cont)), 
    
    mo_8_2diPAP_cat  = as.numeric(fct_recode(mo_8_2diPAP_cat, "1" = "<LOD", "2" = "LOD-LOQ", "3" = ">LOQ")),
    mo__8_2diPAP = as.numeric(mo__8_2diPAP), 
    mo_8_2diPAP_LOD  = as.numeric(mo_8_2diPAP_LOD),
    mo_8_2diPAP_cont = ifelse(mo_8_2diPAP_cat > 1, mo__8_2diPAP, mo_8_2diPAP_LOD/sqrt(2)), 
    mo_8_2diPAP_cont = as.numeric(mo_8_2diPAP_cont), 
    mo_8_2diPAP_cont_log = scale(log(mo_8_2diPAP_cont))) %>%
  select(ident, 
         "ch_ETPA_total_cont_log_M2",
         
         "ch_PRPA_total_cont_log_M2",
         
         "mo_BUPA_total_cont_log_t2",
         "mo_BUPA_total_cont_log_t3",
         "ch_BUPA_total_cont_log_M2", 
         "ch_BUPA_total_cont_log_Y1", 
         
         "mo_BPS_total_cont_log_t2",
         "mo_BPS_total_cont_log_t3",
         "ch_BPS_total_cont_log_M2", 
         "ch_BPS_total_cont_log_Y1", 
         
         "mo_PFDoDa_cont_log", 
         "mo_PFHxPA_cont_log",
         "mo_PFHpA_cont_log",
         "mo_PFTrDa_cont_log", 
         "mo_PFBS_cont_log",
         "mo_PFOSA_cont_log",
         "mo_6_2diPAP_cont_log",
         "mo_8_2diPAP_cont_log")


bdd_exposure <- left_join(bdd_exposure_num, 
                          bdd_exposure_cat, 
                          by = "ident")%>%                   
  select(ident, 
         "mo_MEPA_total_i_cor_t2_ln", "mo_MEPA_total_i_cor_t3_ln", "ch_MEPA_total_i_cor_M2_ln", "ch_MEPA_conj_i_cor_Y1_ln",
         "mo_ETPA_total_i_cor_t2_ln", "mo_ETPA_total_i_cor_t3_ln", "ch_ETPA_total_cont_log_M2", "ch_ETPA_total_i_cor_Y1_ln",
         "mo_PRPA_total_i_cor_t2_ln", "mo_PRPA_total_i_cor_t3_ln", "ch_PRPA_total_cont_log_M2", "ch_PRPA_total_i_cor_Y1_ln",
         "mo_BUPA_total_cont_log_t2", "mo_BUPA_total_cont_log_t3", "ch_BUPA_total_cont_log_M2", "ch_BUPA_total_cont_log_Y1", 
         "mo_BPA_total_i_cor_t2_ln","mo_BPA_total_i_cor_t3_ln", "ch_BPA_conj_i_cor_M2_ln", "ch_BPA_conj_i_cor_Y1_ln",  
         "mo_BPS_total_cont_log_t2", "mo_BPS_total_cont_log_t3", "ch_BPS_total_cont_log_M2", "ch_BPS_total_cont_log_Y1",
         "mo_OXBE_total_i_cor_t2_ln", "mo_OXBE_total_i_cor_t3_ln", "ch_OXBE_total_i_cor_M2_ln", "ch_OXBE_total_i_cor_Y1_ln",
         "mo_TRCS_total_i_cor_t2_ln", "mo_TRCS_total_i_cor_t3_ln","ch_TRCS_total_i_cor_M2_ln", "ch_TRCS_total_i_cor_Y1_ln", 
         "mo_PFOA_cor_ln", 
         "mo_PFNA_ln",
         "mo_PFHxS_cor_ln",
         "mo_PFOS_cor_ln",
         "mo_PFDA_i_cor_ln",
         "mo_PFUnDA_i_cor_ln",
         "mo_PFHpS_i_cor_ln",
         "mo_PFDoDa_cont_log",
         "mo_PFHxPA_cont_log",
         "mo_PFHpA_cont_log",
         "mo_PFTrDa_cont_log",
         "mo_PFBS_cont_log",         
         "mo_PFOSA_cont_log",
         "mo_6_2diPAP_cont_log",
         "mo_8_2diPAP_cont_log") 

## possibilité 4 : on ne prend que les expositions qui sont numériques, séparée par fenetre d'exposition ----
bdd_exposure_t2_num <- bdd_exposure %>%          
  select(
    ident,
    !contains(c("t3", "M2", "Y1")))%>%           # variables d'exposition (déjà logtransformé)
  select(
    ident, 
    !contains("cont"))            # conserver que les expositions codées en numérique)           


bdd_exposure_t2_all <- bdd_exposure %>%          
  select(
    ident,
    !contains(c("t3", "M2", "Y1")))         # variables d'exposition (déjà logtransformé)


bdd_exposure_t3_phenols_num <- bdd_exposure %>%          
  select(
    ident,
    contains("t3"))%>%           # variables d'exposition (déjà logtransformé)
  select(
    ident, 
    !contains("cont"))          # conserver que les expositions codées en numérique)           


bdd_exposure_t3_phenols_all <- bdd_exposure %>%          
  select(
    ident,
    contains("t3"))          # variables d'exposition (déjà logtransformé)



bdd_exposure_M2_phenols_num <- bdd_exposure %>%          
  select(
    ident,
    contains("M2"))%>%           # variables d'exposition (déjà logtransformé)
  select(
    ident, 
    !contains("cont"))           # conserver que les expositions codées en numérique)           
 

bdd_exposure_M2_phenols_all <- bdd_exposure %>%          
  select(
    ident,
    contains("M2"))         # variables d'exposition (déjà logtransformé)



bdd_exposure_Y1_phenols_num <- bdd_exposure %>%          
  select(
    ident,
    contains("Y1"))%>%           # variables d'exposition (déjà logtransformé)
  select(
    ident, 
    !contains("cont"))          # conserver que les expositions codées en numérique)           


bdd_exposure_Y1_phenols_all <- bdd_exposure %>%          
  select(
    ident,
    contains("Y1"))


# Préparation des matrices ----
keep_expo_all <- colnames(bdd_exposure[, 2:48])
keep_expo_num <- colnames(bdd_exposure_num[, 2:30])
keep_expo_t2_all <- colnames(bdd_exposure_t2_all[, 2:24])
keep_expo_t2_num <- colnames(bdd_exposure_t2_num[, 2:14])

keep_expo_t3_all <- colnames(bdd_exposure_t3_phenols_all[, 2:9])
keep_expo_t3_num <- colnames(bdd_exposure_t3_phenols_num[, 2:7])

keep_expo_M2_all <- colnames(bdd_exposure_M2_phenols_all[, 2:9])
keep_expo_M2_num <- colnames(bdd_exposure_M2_phenols_num[, 2:5])

keep_expo_Y1_all <- colnames(bdd_exposure_Y1_phenols_all[, 2:9])
keep_expo_Y1_num <- colnames(bdd_exposure_Y1_phenols_num[, 2:7])

keep_covar <- colnames(bdd_covariates_i[, 2:25])

bdd_bkmr_alpha_all <- 
  list(bdd_outcomes_alpha, bdd_covariates_i, bdd_exposure) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

outcome_specrich <- bdd_bkmr_alpha_all %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon <- bdd_bkmr_alpha_all %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()
outcome_faith <- bdd_bkmr_alpha_all %>% select(ch_feces_Faith_5000_ASV_Y1) %>% as.matrix()
outcome_vec <- c("outcome_specrich", "outcome_shannon", "outcome_faith")

covariates <- bdd_bkmr_alpha_all %>% select(all_of(keep_covar)) %>% as.matrix()

mixture_all <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_all)) %>% as.matrix()
mixture_num <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_num)) %>% as.matrix()

mixture_t2_all <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_t2_all)) %>% as.matrix()
mixture_t2_num <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_t2_num)) %>% as.matrix()
mixture_t3_all <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_t3_all)) %>% as.matrix()
mixture_t3_num <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_t3_num)) %>% as.matrix()
mixture_M2_all <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_M2_all)) %>% as.matrix()
mixture_M2_num <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_M2_num)) %>% as.matrix()
mixture_Y1_all <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_Y1_all)) %>% as.matrix()
mixture_Y1_num <- bdd_bkmr_alpha_all %>% select(all_of(keep_expo_Y1_num)) %>% as.matrix()
mixture_window_vec <- c("mixture_all", "mixture_num", 
                        "mixture_t2_all", "mixture_t2_num", 
                        "mixture_t3_all", "mixture_t3_num", 
                        "mixture_M2_all", "mixture_M2_num",
                        "mixture_Y1_all", "mixture_Y1_num")


# Fit BKMR ----
## BKMR not hierarchical ----
plan(multisession)                    # Spécifier le planificateur future et le nombre de noyaux à utiliser
ncores <- availableCores()            # ou spécifier le nombre de noyaux à utiliser
futures_not_hierar <-                 # future_lapply() pour exécuter bkmr en parallèle pour les 3 outcomes
  future_lapply(list(outcome_specrich, 
                     outcome_shannon, 
                     outcome_faith), 
                bkmr_not_hierar, 
                future.seed = TRUE)

bkmr_sperich_not_hierar_num <- futures_not_hierar[[1]][[1]]   # extraire les résultats
bkmr_sperich_not_hierar_all <- futures_not_hierar[[1]][[2]]

bkmr_shannon_not_hierar_num <- futures_not_hierar[[2]][[1]]
bkmr_shannon_not_hierar_all <- futures_not_hierar[[2]][[2]]

bkmr_faith_not_hierar_num <- futures_not_hierar[[3]][[1]]
bkmr_faith_not_hierar_all <- futures_not_hierar[[3]][[2]]


## BKMR hierarchical ----
plan(multisession)         # Spécifier le planificateur future et le nombre de noyaux à utiliser
ncores <- availableCores() # ou spécifier le nombre de noyaux à utiliser
futures_4 <- 
  future_lapply(list(outcome_specrich,  # future_lapply() pour exécuter bkmr en parallèle pour les 3 outcomes
                     outcome_shannon, 
                     outcome_faith),   
                bkmr_hierar,
                future.seed = TRUE)

bkmr_sperich_hierar_num_family <- futures_4[[1]][[1]]   # extraire les résultats
bkmr_sperich_hierar_all_family <- futures_4[[1]][[2]]
bkmr_sperich_hierar_num_window <- futures_4[[1]][[3]]
bkmr_sperich_hierar_all_window <- futures_4[[1]][[4]]

bkmr_shannon_hierar_num_family <- futures_4[[2]][[1]]
bkmr_shannon_hierar_all_family <- futures_4[[2]][[2]]
bkmr_shannon_hierar_num_window <- futures_4[[2]][[3]]
bkmr_shannon_hierar_all_window <- futures_4[[2]][[4]]

bkmr_faith_hierar_num_family <- futures_4[[3]][[1]]
bkmr_faith_hierar_all_family <- futures_4[[3]][[2]]
bkmr_faith_hierar_num_window <- futures_4[[3]][[3]]
bkmr_faith_hierar_all_window <- futures_4[[3]][[4]]

## BKMR hierarchical bis ----
plan(multisession)         # Spécifier le planificateur future et le nombre de noyaux à utiliser
ncores <- availableCores() # ou spécifier le nombre de noyaux à utiliser
futures_hierar_bis <- 
  future_lapply(list(outcome_specrich,  # future_lapply() pour exécuter bkmr en parallèle pour les 3 outcomes
                     outcome_shannon, 
                     outcome_faith),   
                bkmr_hierar_bis,
                future.seed = TRUE)

bkmr_sperich_hierar_num_window_5 <- futures_hierar_bis[[1]][[1]]
bkmr_sperich_hierar_all_window_5 <- futures_hierar_bis[[1]][[2]]
bkmr_sperich_sep_window_hierar_num_t2 <- futures_hierar_bis[[1]][[3]]
bkmr_sperich_sep_window_hierar_all_t2 <- futures_hierar_bis[[1]][[4]]

bkmr_shannon_hierar_num_window_5 <- futures_hierar_bis[[2]][[1]]
bkmr_shannon_hierar_all_window_5 <- futures_hierar_bis[[2]][[2]]
bkmr_shannon_sep_window_hierar_num_t2 <- futures_hierar_bis[[2]][[3]]
bkmr_shannon_sep_window_hierar_all_t2 <- futures_hierar_bis[[2]][[4]]

bkmr_faith_hierar_num_window_5 <- futures_hierar_bis[[3]][[1]]
bkmr_faith_hierar_all_window_5 <- futures_hierar_bis[[3]][[2]]
bkmr_faith_sep_window_hierar_num_t2 <- futures_hierar_bis[[3]][[3]]
bkmr_faith_sep_window_hierar_all_t2 <- futures_hierar_bis[[3]][[4]]

## BKMR fenetres séparées ----                    
plan(multisession)                    # Spécifier le planificateur future et le nombre de noyaux à utiliser
ncores <- availableCores()            # ou spécifier le nombre de noyaux à utiliser
futures_sep_window <-                 # future_lapply() pour exécuter bkmr en parallèle pour toutes les fenetres 
  future_lapply(list(mixture_t2_all, mixture_t2_num,  
                     mixture_t3_all, mixture_t3_num,
                     mixture_M2_all, mixture_M2_num,
                     mixture_Y1_all, mixture_Y1_num), 
                bkmr_sep_window, 
                future.seed = TRUE)

bkmr_sperich_sep_window_t2_all <- futures_sep_window[[1]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_t2_all <- futures_sep_window[[1]][[2]]
bkmr_faith_sep_window_t2_all <- futures_sep_window[[1]][[3]]

bkmr_sperich_sep_window_t2_num <- futures_sep_window[[2]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_t2_num <- futures_sep_window[[2]][[2]]
bkmr_faith_sep_window_t2_num <- futures_sep_window[[2]][[3]]

bkmr_sperich_sep_window_t3_all <- futures_sep_window[[3]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_t3_all <- futures_sep_window[[3]][[2]]
bkmr_faith_sep_window_t3_all <- futures_sep_window[[3]][[3]]

bkmr_sperich_sep_window_t3_num <- futures_sep_window[[4]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_t3_num <- futures_sep_window[[4]][[2]]
bkmr_faith_sep_window_t3_num <- futures_sep_window[[4]][[3]]

bkmr_sperich_sep_window_M2_all <- futures_sep_window[[5]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_M2_all <- futures_sep_window[[5]][[2]]
bkmr_faith_sep_window_M2_all <- futures_sep_window[[5]][[3]]

bkmr_sperich_sep_window_M2_num <- futures_sep_window[[6]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_M2_num <- futures_sep_window[[6]][[2]]
bkmr_faith_sep_window_M2_num <- futures_sep_window[[6]][[3]]

bkmr_sperich_sep_window_Y1_all <- futures_sep_window[[7]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_Y1_all <- futures_sep_window[[7]][[2]]
bkmr_faith_sep_window_Y1_all <- futures_sep_window[[7]][[3]]

bkmr_sperich_sep_window_Y1_num <- futures_sep_window[[8]][[1]]   # extraire les résultats
bkmr_shannon_sep_window_Y1_num <- futures_sep_window[[8]][[2]]
bkmr_faith_sep_window_Y1_num <- futures_sep_window[[8]][[3]]

save(futures_4, futures_hierar_bis, futures_sep_window, file = "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/pollutants_gut_microbiota_Y1/sauvegarde.bkmr.RData")

# PARTIE 2: PRESENTER LES RESULTATS ----
load("2_final_data/bdd_bkmr_phenols_pfas.RData")

# Model convergence ----
## BKMR hierachical ----
TracePlot_group(bkmr_sperich_hierar_num_window_5, bkmr_shannon_hierar_num_window_5, bkmr_faith_hierar_num_window_5, 
                titre = "Modèles BKMR hierarchical 5 groupes (phenols t2, t3, M2, Y1 et pfas), expo numériques uniquement (Specific richness, Shannon et Faith de haut en bas)")

## BKMR fenetres séparées ----
TracePlot_group(bkmr_sperich_sep_window_t2_all, bkmr_shannon_sep_window_t2_all, bkmr_faith_sep_window_t2_all, 
                titre = "Modèles BKMR trim.2, toutes expositions (Specific richness, Shannon et Faith de haut en bas)" ) 
TracePlot_group(bkmr_sperich_sep_window_t3_all, bkmr_shannon_sep_window_t3_all, bkmr_faith_sep_window_t3_all, 
                titre = "Modèles BKMR trim.3, toutes expositions (Specific richness, Shannon et Faith de haut en bas)" ) 
TracePlot_group(bkmr_sperich_sep_window_M2_all, bkmr_shannon_sep_window_M2_all, bkmr_faith_sep_window_M2_all, 
                titre = "Modèles BKMR at 2 months, toutes expositions (Specific richness, Shannon et Faith de haut en bas)" ) 
TracePlot_group(bkmr_sperich_sep_window_Y1_all, bkmr_shannon_sep_window_Y1_all, bkmr_faith_sep_window_Y1_all, 
                titre = "Modèles BKMR at 12 months, toutes expositions (Specific richness, Shannon et Faith de haut en bas)" ) 

TracePlot_group(bkmr_sperich_sep_window_t2_num, bkmr_shannon_sep_window_t2_num, bkmr_faith_sep_window_t2_num, 
                titre = "Modèles BKMR trim.2, expo numériques uniquement (Specific richness, Shannon et Faith de haut en bas)" ) 
TracePlot_group(bkmr_sperich_sep_window_t3_num, bkmr_shannon_sep_window_t3_num, bkmr_faith_sep_window_t3_num, 
                titre = "Modèles BKMR trim.3, expo numériques uniquement (Specific richness, Shannon et Faith de haut en bas)" ) 
TracePlot_group(bkmr_sperich_sep_window_M2_num, bkmr_shannon_sep_window_M2_num, bkmr_faith_sep_window_M2_num, 
                titre = "Modèles BKMR at 2 months, expo numériques uniquement (Specific richness, Shannon et Faith de haut en bas)" ) 
TracePlot_group(bkmr_sperich_sep_window_Y1_num, bkmr_shannon_sep_window_Y1_num, bkmr_faith_sep_window_Y1_num, 
                titre = "Modèles BKMR at 12 months, expo numériques uniquement (Specific richness, Shannon et Faith de haut en bas)" ) 

TracePlot_group(bkmr_sperich_sep_window_hierar_all_t2, bkmr_shannon_sep_window_hierar_all_t2, bkmr_faith_sep_window_hierar_all_t2, 
                titre = "Modèles BKMR trim.2 (2groupes phenols / pfas), toutes expositions (Specific richness, Shannon et Faith de haut en bas)" ) 
TracePlot_group(bkmr_sperich_sep_window_hierar_num_t2, bkmr_shannon_sep_window_hierar_num_t2, bkmr_faith_sep_window_hierar_num_t2, 
                titre = "Modèles BKMR trim.2 (2groupes phenols / pfas), expo numériques uniquement (Specific richness, Shannon et Faith de haut en bas)" ) 

par(mfrow=c(1,1)) # Revenir à une fenêtre graphique normale

# BKMR stats ----
## PIP ----
pip_results_hierar <- list(
  
  pip_hierar_num_family = pip_results_hierar(bkmr_sperich_hierar_num_family, bkmr_shannon_hierar_num_family, bkmr_faith_hierar_num_family),
  pip_hierar_all_family = pip_results_hierar(bkmr_sperich_hierar_all_family, bkmr_shannon_hierar_all_family, bkmr_faith_hierar_all_family),
  pip_hierar_num_window = pip_results_hierar(bkmr_sperich_hierar_num_window, bkmr_shannon_hierar_num_window, bkmr_faith_hierar_num_window),
  pip_hierar_all_window = pip_results_hierar(bkmr_sperich_hierar_all_window, bkmr_shannon_hierar_all_window, bkmr_faith_hierar_all_window),
  pip_hierar_num_window_5 = pip_results_hierar(bkmr_sperich_hierar_num_window_5, bkmr_shannon_hierar_num_window_5, bkmr_faith_hierar_num_window_5),
  pip_hierar_all_window_5 = pip_results_hierar(bkmr_sperich_hierar_all_window_5, bkmr_shannon_hierar_all_window_5, bkmr_faith_hierar_all_window_5), 
  
  pip_sep_window_hierar_num_t2 = pip_results_hierar(bkmr_sperich_sep_window_hierar_num_t2, bkmr_sperich_sep_window_hierar_num_t2, bkmr_sperich_sep_window_hierar_num_t2), 
  pip_sep_window_hierar_all_t2 = pip_results_hierar(bkmr_sperich_sep_window_hierar_all_t2, bkmr_sperich_sep_window_hierar_all_t2, bkmr_sperich_sep_window_hierar_all_t2))

pip_results_not_hierar <- list(
  
  pip_sep_window_num_t2 = pip_results_not_hierar(bkmr_sperich_sep_window_t2_num, bkmr_shannon_sep_window_t2_num, bkmr_faith_sep_window_t2_num), 
  pip_sep_window_num_t3 = pip_results_not_hierar(bkmr_sperich_sep_window_t3_num, bkmr_shannon_sep_window_t3_num, bkmr_faith_sep_window_t3_num), 
  pip_sep_window_num_M2 = pip_results_not_hierar(bkmr_sperich_sep_window_M2_num, bkmr_shannon_sep_window_M2_num, bkmr_faith_sep_window_M2_num), 
  pip_sep_window_num_Y1 = pip_results_not_hierar(bkmr_sperich_sep_window_Y1_num, bkmr_shannon_sep_window_Y1_num, bkmr_faith_sep_window_Y1_num),
  
  pip_sep_window_all_t2 = pip_results_not_hierar(bkmr_sperich_sep_window_t2_all, bkmr_shannon_sep_window_t2_all, bkmr_faith_sep_window_t2_all), 
  pip_sep_window_all_t3 = pip_results_not_hierar(bkmr_sperich_sep_window_t3_all, bkmr_shannon_sep_window_t3_all, bkmr_faith_sep_window_t3_all), 
  pip_sep_window_all_M2 = pip_results_not_hierar(bkmr_sperich_sep_window_M2_all, bkmr_shannon_sep_window_M2_all, bkmr_faith_sep_window_M2_all), 
  pip_sep_window_all_Y1 = pip_results_not_hierar(bkmr_sperich_sep_window_Y1_all, bkmr_shannon_sep_window_Y1_all, bkmr_faith_sep_window_Y1_all)) 

write_xlsx(pip_results_hierar, path = "4_output/pip_results_hierar.xlsx")
write_xlsx(pip_results_not_hierar, path = "4_output/pip_results_not_hierar.xlsx")
  


## Overall ----
### Tables ----
#bkmr_risks_overall_results_specrich_all_window <- risks_overall(bkmr_sperich_hierar_all_window, outcome_specrich, mixture_all)
#bkmr_risks_overall_results_shannon_all_window <- risks_overall(bkmr_shannon_hierar_all_window, outcome_shannon, mixture_all)
#bkmr_risks_overall_results_faith_all_window <- risks_overall(bkmr_faith_hierar_all_window, outcome_faith, mixture_all)

#bkmr_risks_overall_results_specrich_num_window <- risks_overall(bkmr_sperich_hierar_num_window, outcome_specrich, mixture_num)
#bkmr_risks_overall_results_shannon_num_window <- risks_overall(bkmr_shannon_hierar_num_window, outcome_shannon, mixture_num)
#bkmr_risks_overall_results_faith_num_window <- risks_overall(bkmr_faith_hierar_num_window, outcome_faith, mixture_num)

#bkmr_risks_overall_results_specrich_all_window_5 <- risks_overall(bkmr_sperich_hierar_all_window_5, outcome_specrich, mixture_all)
#bkmr_risks_overall_results_shannon_all_window_5 <- risks_overall(bkmr_shannon_hierar_all_window_5, outcome_shannon, mixture_all)
#bkmr_risks_overall_results_faith_all_window_5 <- risks_overall(bkmr_faith_hierar_all_window_5, outcome_faith, mixture_all)

bkmr_risks_overall_results_specrich_num_window_5 <- risks_overall(bkmr_sperich_hierar_num_window_5, outcome_specrich, mixture_num)
bkmr_risks_overall_results_shannon_num_window_5 <- risks_overall(bkmr_shannon_hierar_num_window_5, outcome_shannon, mixture_num)
bkmr_risks_overall_results_faith_num_window_5 <- risks_overall(bkmr_faith_hierar_num_window_5, outcome_faith, mixture_num)

#bkmr_risks_overall_results_specrich_all_family <- risks_overall(bkmr_sperich_hierar_all_family, outcome_specrich, mixture_all)
#bkmr_risks_overall_results_shannon_all_family <- risks_overall(bkmr_shannon_hierar_all_family, outcome_shannon, mixture_all)
#bkmr_risks_overall_results_faith_all_family <- risks_overall(bkmr_faith_hierar_all_family, outcome_faith, mixture_all)

#bkmr_risks_overall_results_specrich_num_family <- risks_overall(bkmr_sperich_hierar_num_family, outcome_specrich, mixture_num)
#bkmr_risks_overall_results_shannon_num_family <- risks_overall(bkmr_shannon_hierar_num_family, outcome_shannon, mixture_num)
#bkmr_risks_overall_results_faith_num_family <- risks_overall(bkmr_faith_hierar_num_family, outcome_faith, mixture_num)

  
#bkmr_risks_overall_specrich_sep_window_all_t2 <- risks_overall(bkmr_sperich_sep_window_t2_all, outcome_specrich, mixture_t2_all)
#bkmr_risks_overall_shannon_sep_window_all_t2 <- risks_overall(bkmr_shannon_sep_window_t2_all, outcome_shannon, mixture_t2_all)
#bkmr_risks_overall_faith_sep_window_all_t2 <- risks_overall(bkmr_faith_sep_window_t2_all, outcome_faith, mixture_t2_all)
#bkmr_risks_overall_specrich_sep_window_num_t2 <- risks_overall(bkmr_sperich_sep_window_t2_num, outcome_specrich, mixture_t2_num)
#bkmr_risks_overall_shannon_sep_window_num_t2 <- risks_overall(bkmr_shannon_sep_window_t2_num, outcome_shannon, mixture_t2_num)
#bkmr_risks_overall_faith_sep_window_num_t2 <- risks_overall(bkmr_faith_sep_window_t2_num, outcome_faith, mixture_t2_num)


#bkmr_risks_overall_specrich_sep_window_hierar_all_t2 <- risks_overall(bkmr_sperich_sep_window_hierar_all_t2, outcome_specrich, mixture_t2_all)
#bkmr_risks_overall_shannon_sep_window_hierar_all_t2 <- risks_overall(bkmr_shannon_sep_window_hierar_all_t2, outcome_shannon, mixture_t2_all)
#bkmr_risks_overall_faith_sep_window_hierar_all_t2 <- risks_overall(bkmr_faith_sep_window_hierar_all_t2, outcome_faith, mixture_t2_all)
bkmr_risks_overall_specrich_sep_window_hierar_num_t2 <- risks_overall(bkmr_sperich_sep_window_hierar_num_t2, outcome_specrich, mixture_t2_num)
bkmr_risks_overall_shannon_sep_window_hierar_num_t2 <- risks_overall(bkmr_shannon_sep_window_hierar_num_t2, outcome_shannon, mixture_t2_num)
bkmr_risks_overall_faith_sep_window_hierar_num_t2 <- risks_overall(bkmr_faith_sep_window_hierar_num_t2, outcome_faith, mixture_t2_num)


#bkmr_risks_overall_specrich_sep_window_all_t3 <- risks_overall(bkmr_sperich_sep_window_t3_all, outcome_specrich, mixture_t3_all)
#bkmr_risks_overall_shannon_sep_window_all_t3 <- risks_overall(bkmr_shannon_sep_window_t3_all, outcome_shannon, mixture_t3_all)
#bkmr_risks_overall_faith_sep_window_all_t3 <- risks_overall(bkmr_faith_sep_window_t3_all, outcome_faith, mixture_t3_all)
bkmr_risks_overall_specrich_sep_window_num_t3 <- risks_overall(bkmr_sperich_sep_window_t3_num, outcome_specrich, mixture_t3_num)
bkmr_risks_overall_shannon_sep_window_num_t3 <- risks_overall(bkmr_shannon_sep_window_t3_num, outcome_shannon, mixture_t3_num)
bkmr_risks_overall_faith_sep_window_num_t3 <- risks_overall(bkmr_faith_sep_window_t3_num, outcome_faith, mixture_t3_num)

#bkmr_risks_overall_specrich_sep_window_all_M2 <- risks_overall(bkmr_sperich_sep_window_M2_all, outcome_specrich, mixture_M2_all)
#bkmr_risks_overall_shannon_sep_window_all_M2 <- risks_overall(bkmr_shannon_sep_window_M2_all, outcome_shannon, mixture_M2_all)
#bkmr_risks_overall_faith_sep_window_all_M2 <- risks_overall(bkmr_faith_sep_window_M2_all, outcome_faith, mixture_M2_all)
bkmr_risks_overall_specrich_sep_window_num_M2 <- risks_overall(bkmr_sperich_sep_window_M2_num, outcome_specrich, mixture_M2_num)
bkmr_risks_overall_shannon_sep_window_num_M2 <- risks_overall(bkmr_shannon_sep_window_M2_num, outcome_shannon, mixture_M2_num)
bkmr_risks_overall_faith_sep_window_num_M2 <- risks_overall(bkmr_faith_sep_window_M2_num, outcome_faith, mixture_M2_num)

#bkmr_risks_overall_specrich_sep_window_all_Y1 <- risks_overall(bkmr_sperich_sep_window_Y1_all, outcome_specrich, mixture_Y1_all)
#bkmr_risks_overall_shannon_sep_window_all_Y1 <- risks_overall(bkmr_shannon_sep_window_Y1_all, outcome_shannon, mixture_Y1_all)
#bkmr_risks_overall_faith_sep_window_all_Y1 <- risks_overall(bkmr_faith_sep_window_Y1_all, outcome_faith, mixture_Y1_all)
bkmr_risks_overall_specrich_sep_window_num_Y1 <- risks_overall(bkmr_sperich_sep_window_Y1_num, outcome_specrich, mixture_Y1_num)
bkmr_risks_overall_shannon_sep_window_num_Y1 <- risks_overall(bkmr_shannon_sep_window_Y1_num, outcome_shannon, mixture_Y1_num)
bkmr_risks_overall_faith_sep_window_num_Y1 <- risks_overall(bkmr_faith_sep_window_Y1_num, outcome_faith, mixture_Y1_num)

### plots ----

#plot_risks.overall_all_family <-
#  plot_risks.overall(
#    risks.overall_sperich = bkmr_risks_overall_results_specrich_all_family,
#    risks.overall_shannon = bkmr_risks_overall_results_shannon_all_family,
#    risks.overall_faith = bkmr_risks_overall_results_faith_all_family,
#    plot_title = "All exposure, hierarchical family (2 groups)")

#plot_risks.overall_all_window <-
#    plot_risks.overall(
#      risks.overall_sperich = bkmr_risks_overall_results_specrich_all_window,
#      risks.overall_shannon = bkmr_risks_overall_results_shannon_all_window,
#      risks.overall_faith = bkmr_risks_overall_results_faith_all_window,
#      plot_title = "All exposure, hierarchical window (3 groups)")

#plot_risks.overall_all_window_5 <-
#  plot_risks.overall(
#    risks.overall_sperich = bkmr_risks_overall_results_specrich_all_window_5,
#    risks.overall_shannon = bkmr_risks_overall_results_shannon_all_window_5,
#    risks.overall_faith = bkmr_risks_overall_results_faith_all_window_5,
#    plot_title = "All exposure, hierarchical window (5 groups)")

#plot_risks.overall_num_family <-
#  plot_risks.overall(
#    risks.overall_sperich = bkmr_risks_overall_results_specrich_num_family,
#    risks.overall_shannon = bkmr_risks_overall_results_shannon_num_family,
#    risks.overall_faith = bkmr_risks_overall_results_faith_num_family,
#    plot_title = "Continuous exposures, hierarchical family (2 groups)")

#plot_risks.overall_num_window <-
#  plot_risks.overall(
#    risks.overall_sperich = bkmr_risks_overall_results_specrich_num_window,
#    risks.overall_shannon = bkmr_risks_overall_results_shannon_num_window,
#    risks.overall_faith = bkmr_risks_overall_results_faith_num_window,
#    plot_title = "Continuous exposures, hierarchical window (3 groups)")

plot_risks.overall_num_window_5 <-
  plot_risks.overall(
    risks.overall_sperich = bkmr_risks_overall_results_specrich_num_window_5,
    risks.overall_shannon = bkmr_risks_overall_results_shannon_num_window_5,
    risks.overall_faith = bkmr_risks_overall_results_faith_num_window_5
    #,plot_title = "Continuous exposures, hierarchical window (5 groups)"
    )

#plot_bkmr_risks_overall_sep_window_hierar_all_t2 <-
#  plot_risks.overall(
#    risks.overall_sperich = bkmr_risks_overall_specrich_sep_window_hierar_all_t2,
#    risks.overall_shannon = bkmr_risks_overall_shannon_sep_window_hierar_all_t2,
#    risks.overall_faith = bkmr_risks_overall_faith_sep_window_hierar_all_t2, 
#    plot_title = "All exposures, hierarchical t2 (phenols/pfas)")

plot_risks.overall <- 
  plot_risks.overall_sperich(bkmr_risks_overall_specrich_sep_window_hierar_num_t2, titre = "Second trimester exposure") + 
  plot_risks.overall_shannon(bkmr_risks_overall_shannon_sep_window_hierar_num_t2) + 
  plot_risks.overall_faith(bkmr_risks_overall_faith_sep_window_hierar_num_t2) + 
  
  plot_risks.overall_sperich(bkmr_risks_overall_specrich_sep_window_num_t3, titre = "Third trimester exposure") + 
  plot_risks.overall_shannon(bkmr_risks_overall_shannon_sep_window_num_t3) + 
  plot_risks.overall_faith(bkmr_risks_overall_faith_sep_window_num_t3) + 
  
  plot_risks.overall_sperich(bkmr_risks_overall_specrich_sep_window_num_M2, titre = "2 months exposure") + 
  plot_risks.overall_shannon(bkmr_risks_overall_shannon_sep_window_num_M2) + 
  plot_risks.overall_faith(bkmr_risks_overall_faith_sep_window_num_M2) + 
  
  plot_risks.overall_sperich(bkmr_risks_overall_specrich_sep_window_num_Y1, titre = "12 months exposure") + 
  plot_risks.overall_shannon(bkmr_risks_overall_shannon_sep_window_num_Y1) + 
  plot_risks.overall_faith(bkmr_risks_overall_faith_sep_window_num_Y1) + 
  
  plot_layout(ncol = 3, nrow = 4) 



  
ggsave("4_output/plot_risks.overall_all_family.png", 
       plot_risks.overall_all_family, 
       device = "png",
       units = "cm",
       width = 30,
       height = 15)

ggsave("4_output/plot_risks.overall_all_window.png", 
       plot_risks.overall_all_window, 
       device = "png",
       units = "cm",
       width = 30,
       height = 15)

ggsave("4_output/plot_risks.overall_all_window_5.png", 
       plot_risks.overall_all_window_5, 
       device = "png",
       units = "cm",
       width = 30,
       height = 15)

ggsave("4_output/plot_risks.overall_num_family.png", 
       plot_risks.overall_num_family, 
       device = "png",
       units = "cm",
       width = 30,
       height = 15)

ggsave("4_output/plot_risks.overall_num_window.png", 
       plot_risks.overall_num_window, 
       device = "png",
       units = "cm",
       width = 30,
       height = 15)

ggsave("4_output/plot_risks.overall_num_window_5.tiff", 
       plot_risks.overall_num_window_5, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 180,
       height = 80)

ggsave("4_output/plot_risks.overall.png", 
       plot_risks.overall, 
       device = "png",
       units = "cm",
       width = 30,
       height = 45)

ggsave("4_output/plot_bkmr_risks_overall_sep_window_hierar_num_t2.png", 
       plot_bkmr_risks_overall_sep_window_hierar_num_t2, 
       device = "png",
       units = "cm",
       width = 30,
       height = 15)




## Singvar ----
### Tables ----
bkmr_risks_singvar_specrich_num_window_5 <- risks_singvar(bkmr_sperich_hierar_num_window_5, outcome_specrich, mixture_num)
bkmr_risks_singvar_shannon_num_window_5 <- risks_singvar(bkmr_shannon_hierar_num_window_5, outcome_shannon, mixture_num)
bkmr_risks_singvar_faith_num_window_5 <- risks_singvar(bkmr_faith_hierar_num_window_5, outcome_faith, mixture_num)


### Plots ----
plot_risks.singvar_num_window_5 <-
    plot_risks.singvar(
      risks.singvar_sperich = bkmr_risks_singvar_specrich_num_window_5,
      risks.singvar_shannon = bkmr_risks_singvar_shannon_num_window_5,
      risks.singvar_faith = bkmr_risks_singvar_faith_num_window_5
      #,plot_title = "Num, 5 groups"
      )
  
ggsave(
  "4_output/plot_risks.singvar_num_window_5.tiff",
  plot = plot_risks.singvar_num_window_5,
  units = "mm",
  device = "tiff",
  dpi = 300,
  width = 180,
  height = 200,
  limitsize = FALSE)
  
  
  plot_risks.singvar_sperich_all_window <- bkmr_risks_singvar_results_specrich_all_window %>%
    ggplot(
      aes(
        variable,
        est,
        ymin = est - 1.96 * sd,
        ymax = est + 1.96 * sd,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    labs(y = "Specific richness", 
         x ="Exposure") + 
    coord_flip()+
    theme_bw() +
    theme(legend.position = "none")
  
  
  
  plot_risks.singvar_sperich_all_window <- bkmr_risks_singvar_results_specrich_all_window %>%
    ggplot(
      aes(
        variable,
        est,
        ymin = est - 1.96 * sd,
        ymax = est + 1.96 * sd,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    labs(y = "Specific richness", 
         x ="Exposure") + 
    coord_flip()+
    theme_bw() +
    theme(legend.position = "none")
  
  
# Sauvegarde des résultats bkmr ----
list_bkmr <- list(futures,    # old
                  futures_2,  # old
                  futures_3,  # old
                  futures_4,  # future hierar
                  futures_hierar_bis, # future hierar 5 groups + future hierar t2
                  futures_sep_window,   # future par fenetre d'expo
                  
                  pip_results_hierar, 
                  pip_results_not_hierar, 
                  
                  bkmr_risks_overall_results_specrich_num_window_5,
                  bkmr_risks_overall_results_shannon_num_window_5,
                  bkmr_risks_overall_results_faith_num_window_5, 
                  
                  bkmr_risks_overall_specrich_sep_window_hierar_num_t2,
                  bkmr_risks_overall_shannon_sep_window_hierar_num_t2,
                  bkmr_risks_overall_faith_sep_window_hierar_num_t2,
                  
                  bkmr_risks_overall_specrich_sep_window_num_t3, 
                  bkmr_risks_overall_shannon_sep_window_num_t3, 
                  bkmr_risks_overall_faith_sep_window_num_t3, 
                  
                  bkmr_risks_overall_specrich_sep_window_num_M2, 
                  bkmr_risks_overall_shannon_sep_window_num_M2, 
                  bkmr_risks_overall_faith_sep_window_num_M2, 
                  
                  bkmr_risks_overall_specrich_sep_window_num_Y1, 
                  bkmr_risks_overall_shannon_sep_window_num_Y1, 
                  bkmr_risks_overall_faith_sep_window_num_Y1, 
                  
                  bkmr_risks_singvar_results_specrich_all_window,
                  bkmr_risks_singvar_results_shannon_all_window,
                  bkmr_risks_singvar_results_faith_all_window,
                  
                  bkmr_risks_singvar_specrich_num_window_5,
                  bkmr_risks_singvar_shannon_num_window_5,
                  bkmr_risks_singvar_faith_num_window_5 )
  
  
save(list_bkmr, 
     file = "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/pollutants_gut_microbiota_Y1/4_output/data_bkmr.RData")


# Figure pour poster SF DohaD ----
plot_risks.overall_sperich <- function(risks.overall_sperich){
  ggplot(risks.overall_sperich,
         aes(
           quantile,
           est,
           ymin = est - 1.96 * sd,
           ymax = est + 1.96 * sd
         )) +
    geom_pointrange(size = 1) +
    labs(y = "Specific richness") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() + 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA), 
      axis.title = element_text(size = 20), 
      axis.text = element_text(size = 20))
}

plot_risks.overall_shannon <-function(risks.overall_shannon)
{   ggplot(risks.overall_shannon,
           aes(
             quantile,
             est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd
           )) +
    geom_pointrange(size = 1) + 
    labs(y = "Shannon diversity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+ 
    theme_bw() + 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA), 
      axis.title = element_text(size = 20), 
      axis.text = element_text(size = 20))
}

plot_risks.overall_faith <- function(risks.overall_faith)
{    ggplot(risks.overall_faith,
            aes(
              quantile,
              est,
              ymin = est - 1.96 * sd,
              ymax = est + 1.96 * sd
            )) +
    geom_pointrange(size = 1) +
    labs(y = "Faith phylogenetic diversity")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+ 
    theme_bw() + 
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),  # enleve tout le fond sauf legend 
      legend.background = element_rect(fill = "transparent", colour = NA), 
      axis.title = element_text(size = 20), 
      axis.text = element_text(size = 20))
} 

plot_risks.overall_num_window_5_specrich <- plot_risks.overall_sperich(bkmr_risks_overall_results_specrich_num_window_5)
plot_risks.overall_num_window_5_shannon <- plot_risks.overall_shannon(bkmr_risks_overall_results_shannon_num_window_5)
plot_risks.overall_num_window_5_faith <- plot_risks.overall_faith(bkmr_risks_overall_results_faith_num_window_5)


ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/plot_risks.overall_num_window_5_specrich.png", 
       plot_risks.overall_num_window_5_specrich, 
       device = "png",
       units = "mm",
       width = 120,
       height = 140)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/plot_risks.overall_num_window_5_shannon.png", 
       plot_risks.overall_num_window_5_shannon, 
       device = "png",
       units = "mm",
       width = 120,
       height = 140)

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/13. SF-DohaD 16.11.2023/plot_risks.overall_num_window_5_faith.png", 
       plot_risks.overall_num_window_5_faith, 
       device = "png",
       units = "mm",
       width = 120,
       height = 140)
