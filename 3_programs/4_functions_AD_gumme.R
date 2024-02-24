## Fonctions Aline
## GUMME project

# Chargement des packages ----
library(tidyverse)
library(haven)
library(reshape2)
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
library(parameters)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(rstatix)
library(grDevices)
library(lazyeval)
library(mice)
library(car)
library(expss)

theme_gtsummary_language("en", decimal.mark = ".", big.mark = " ")
theme_gtsummary_compact(set_theme = TRUE)


# Descriptif -----
## Tableau descriptif variables numériques
descrip_num <- function(data, vars) {                                           
  data %>%
    select(all_of({{vars}})) %>%
    tbl_summary(
      missing = "no", 
      type = list(where(is.numeric)~ "continuous"), 
      statistic = all_continuous()~ "{min}/{p25}/{median}/{mean}/{p75}/{max}/{N_nonmiss}", 
      digits = all_continuous() ~ 3
    ) %>%
    bold_labels() %>% 
    as_gt()  %>% 
    as.data.frame()%>% 
    select(variable, label, stat_0) %>%
    separate(
      col = stat_0, 
      into = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N"), 
      sep = "/", 
      remove = TRUE) %>%
    rename(
      "Variable names" = variable, 
      "Variable labels" = label) 
  #%>% kable()
}

## Tableau comparaison d'effectifs
comp_effectifs <- function(data, vars_col1, vars_col2, name_col1, name_col2){ 
  table_col1 <- data %>% select(all_of({{vars_col1}})) 
  table_col2 <- data %>% select(all_of({{vars_col2}}))
  colnames(table_col2) <- colnames(table_col1)
  
  table_col1 <- table_col1 %>% tbl_summary()
  table_col2 <- table_col2 %>% tbl_summary()
  comp <- tbl_merge(
    tbls = list(table_col1, table_col2), 
    tab_spanner = c({{name_col1}}, {{name_col2}}))
  return(comp)
}


## Scatterplots
scatterplot <- function(data, outcome, vars) {
  data_long <-
    data %>%
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric)
  data_long <- data_long %>%
    pivot_longer(cols = -ident,
                 names_to = "Variable",
                 values_to = "Value")
  bdd_outcome <-
    data %>% 
    select(ident, {{outcome}})
  data_long <- 
    left_join(data_long,
              bdd_outcome,
              by = "ident") %>%
    rename(Outcome = {{outcome}})
  
  scatterplot <- data_long %>%
    ggplot() +
    aes(x = Outcome, y = Value) +
    geom_point(
      shape = "circle",
      size = 1.55,
      colour = "#112446") +
    labs(y ="") +
    theme_bw() +
    facet_wrap(~Variable, scales = "free", ncol = 4L)
  return(scatterplot)
    
}

## Boxplots
boxplot <- function(data, vars) {                                               # data = df large dans lequel récuperer les variables                                                                                         # vars = vecteur de noms de variables (toutes numériques)
  data_long <-                                                     
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  boxplot <- data_long %>%
    ggplot() +
    aes(x = "", y = Value) +
    geom_boxplot(shape = "circle", fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap(~Variable, scales = "free", ncol = 2L)
  
  return(boxplot)
}

## Histogrammes
histogram <- function(data, vars, order_var) {                                  # data = df large dans lequel récuperer les variables  
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes numériques)
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = fct_relevel(Variable, order_var))
  
  histogram <- 
    ggplot(data_long) +
    aes(x = Value) +
    geom_histogram(bins = 30L, fill = "#112446") +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap( ~ Variable, scales = "free", ncol = 2L)
  
  return(histogram)
} 

## Barplots 
barplot <- function(data, vars) {                                               # data = df large dans lequel récuperer les variables 
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes catégorielles)
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], factor) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  barplot <- data_long %>%
    filter(!is.na(Value)) %>%
    ggplot() +
    aes(x = Value) +
    geom_bar(fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap( ~ Variable, scales = "free", ncol = 3L) 
  
  return(barplot)
}

## Densityplots 
densityplot <- function(data, vars) {                                           # data = df large dans lequel récuperer les variables  
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes numériques)
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  densityplot <- 
    ggplot(data_long) +
    aes(x = Value) +
    geom_density(fill = "lightgray") +
    theme_bw() +
    facet_wrap( ~ Variable, scales = "free", ncol = 3L) +
    stat_overlay_normal_density(color = "red", linetype = "dashed") +
    theme(axis.title = element_blank())
  
  return(densityplot)
} 


## Heatmap de correlation 
heatmap_cor <- function(cormat, decimal) {                  # cormat = df avec seulement les variables à mettre dans l'heatmap + des noms raccourcis
  
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){              # Utiliser la corrélation entre les variables comme mesure de distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  cormat <- round(cor(cormat, 
                      use = "pairwise.complete.obs", 
                      method = "pearson"), decimal)
  cormat <- reorder_cormat(cormat)                 # réordonner les coef de cor
  upper_tri <- get_upper_tri(cormat)               # obtenir que le triangle sup
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)   # passer en df long rapidement 
  
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
      name = "Pearson\nCorrelation"
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
  
  return(heatmap)
}



heatmap_cor_pairwise <- function(data, vars_1, vars_2, decimal){
  
  
  bdd_cormat <- data %>% select(all_of({{vars_1}}), all_of({{vars_2}}))
  
  cormat <- round(cor(bdd_cormat, 
                      use = "pairwise.complete.obs", 
                      method = "spearman"), decimal)
  
  cormat <- cormat %>% 
    as.data.frame() %>% 
    select(all_of({{vars_1}})) %>% 
    t() %>%
    as.data.frame() %>%
    select(all_of({{vars_2}})) %>%
    as.matrix()
  
  # colnames(cormat) <- colnames(cormat) %>%
  #   str_replace_all(  
  #     c("mo_" = "",
  #       "ch_" = "",
  #       "_total_i_cor_sg_" = " ",
  #       "_i_cor_sg_" = " "))
  rownames(cormat) <- rownames(cormat) %>%
    str_replace_all(
      c("mo_" = "",
        "ch_" = "",
        "_total_i_cor_" = " ",
        "_i_cor_" = " ", 
        "_ms" = "", 
        "DINCH" = "ΣDINCH",
        "DiNP" = "ΣDiNP",
        "DEHP" = "ΣDEHP"))
  
  # cormat <- cormat %>%
  #   reorder_cormat()  %>%    # réordonner les coef de cor
  #   get_upper_tri() %>%      # obtenir que le triangle sup
  #   as.matrix()
  
  
  cormat_long <- 
    reshape2::melt(cormat, na.rm = TRUE) %>% # passer en df long rapidement 
    mutate(
      Var1 = fct_relevel(
        Var1, 
        "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", 
        "ohMPHP Y1", "ohMPHP t3", "ohMPHP t2", 
        "MEP Y1", "MEP M2", "MEP t2", "MEP t3", 
        "MBzP Y1", "MBzP M2", "MBzP t3", "MBzP t2", 
        "MiBP Y1", "MiBP M2", "MiBP t3", "MiBP t2", 
        "ΣDiNP Y1", "ΣDiNP M2", "ΣDiNP t3", "ΣDiNP t2", 
        "MnBP Y1", "MnBP M2", "MnBP t3", "MnBP t2", 
        "ΣDEHP Y1", "ΣDEHP M2", "ΣDEHP t3", "ΣDEHP t2"))
  
  heatmap <-                                       # faire la heatmap
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
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 12,
      hjust = 1
    )) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value),
              color = "black",
              size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(1.5, 0.5),
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 7,
      barheight = 1,
      title.position = "top",
      title.hjust = 0.5
    ))
  return(heatmap)
}



## Visualiser les outliers d'une variable + son ident 
outliers <- function(data, var) {
  data %>% 
    select(ident, {{var}}) %>% 
    filter({{var}} %in% boxplot.stats({{var}})$out) %>% 
    arrange({{var}}) %>% 
    kable()
}


## Visualiser des corrélations entre deux groupes de variables identiques (ex: polluant i_cor VS polluants i_cor_sg)
table_cor_sg <- function(data, exposure_vec, exposure_vec_sg) {
  
  bdd_cor_sg <- data %>% 
    select(all_of({{exposure_vec}}), 
           all_of({{exposure_vec_sg}}))
  
  table_cor_sg <- 
    round(cor(bdd_cor_sg, 
              use = "pairwise.complete.obs", 
              method = "spearman"), 
          2) %>%
    as.data.frame() %>% 
    select(all_of({{exposure_vec}})) %>%
    t() %>%
    as.data.frame() %>%
    select(all_of({{exposure_vec_sg}})) %>%
    as.matrix()
  
  colnames(table_cor_sg) <- colnames(table_cor_sg) %>%
    str_replace_all(  
      c("mo_" = "",
        "ch_" = "",
        "_total_i_cor_sg_" = " ",
        "_i_cor_sg_" = " "))
  rownames(table_cor_sg) <- rownames(table_cor_sg) %>%
    str_replace_all(
      c("mo_" = "",
        "ch_" = "",
        "_total_i_cor_" = " ",
        "_i_cor_" = " "))
  table_cor_sg <- table_cor_sg %>%
    melt(na.rm = TRUE) %>%  # passer en df long rapidement 
    filter(Var1 == Var2)
  
  return(table_cor_sg)
}


## Visualiser des corrélations entre une variable et un groupe de variables (ex: xx_pool_sg_xx VS polluants)
table_cor <- function(data, var, vars){
  bdd_cor <- data %>% 
    select(all_of({{vars}}), 
           {{var}})
  table_cor <- 
    round(cor(bdd_cor, 
              method = "spearman",
              use = "pairwise.complete.obs"), 2) %>%
    as.data.frame() %>% 
    select(-{{var}}) %>%
    t()%>%
    as.data.frame()%>%
    select({{var}}) 
  return(table_cor)
}



# Statistiques (main analysis) ----
## Modeles (adjusted on MiSeq batch effect and fully adjusted) ----
model_summary <- function(data, outcome, exposure){
  
  model_univ <-            # modèle seulement ajusté sur MiSeq batch effect 
    lm(outcome ~ 
         exposure + 
         ch_feces_RUN_Y1, 
       data = data)
  
  ks.test_model_univ <- 
    ks.test(rstudent(model_univ), 
            pnorm, 
            mean = mean(rstudent(model_univ)), 
            sd = sd(rstudent(model_univ)))  
  
  model_multi <-           # modèle complet
    lm(
      outcome ~
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
      data = data
    )
  ks.test_model_multi <- 
    ks.test(rstudent(model_multi), pnorm, mean = mean(rstudent(model_multi)), sd = sd(rstudent(model_multi)))  
  
  qqnorm <- qqnorm(rstudent(model_multi)) ; abline(0, 1) 
  
  return(list("Adjusted on the MiSeq run model" = summary(model_univ),
              "Residuals normality test, adjusted on the MiSeq run model" = ks.test_model_univ,
              "Fully adjusted model" = summary(model_multi),
              "Residuals normality test, fully adjusted model" = ks.test_model_multi,
              view(qqnorm)))
}



## Tableau de resultats (gt_summary) ----
### Tableaux de resultats effets des coavriables seules (effectif + beta + IC95% + p fully adjusted)
model_covar <- function(data, outcome){
  model_multi <-           # modèle complet
    lm(
      outcome ~
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
  result <- model_multi %>% 
    tbl_regression(
      hide_n = FALSE,
      pvalue_fun = ~ style_pvalue(.x, digits = 3),
      estimate_fun = ~ style_sigfig(.x, digits = 2)) %>%
    add_global_p(keep = TRUE, singular.ok = TRUE) %>%
    bold_labels() 
  return(result)
} 

### Tableaux de résultats (effectif + beta + IC95% + p avant et après ajustement)
model_univ_multi <- function(data, outcome, exposure_vec){
  
  effectif <- data %>%                   # Création d'une colonne effectif
    filter(!is.na({{outcome}}))%>%
    select(all_of({{exposure_vec}})) %>%
    tbl_summary(missing = "no", 
                statistic = all_continuous() ~ "{N_nonmiss}") %>%
    bold_labels()
  
  univ <- data %>%                      # Création des colonnes modèle simple ajusté sur le run
    select(
      {{outcome}}, 
      all_of({{exposure_vec}}),
      ch_feces_RUN_Y1                         
    ) %>%
    
    tbl_uvregression(
      method = lm ,
      y = {{outcome}},
      formula = "{y} ~ {x} + ch_feces_RUN_Y1", 
      hide_n = TRUE, 
      pvalue_fun = ~style_pvalue(.x, digits = 3), 
      estimate_fun = ~style_sigfig(.x,digits = 2)) %>%
    add_global_p() %>%
    bold_labels() 
  
  multi <- data %>%                  # Création des colonnes modèle complet
    select(
      {{outcome}},
      all_of({{exposure_vec}}),
      all_of(covar_vec_i)) %>% 
    
    tbl_uvregression(
      method = lm ,
      y = {{outcome}},
      formula = "{y} ~ {x} +
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
      estimate_fun = ~ style_sigfig(.x, digits = 1)
    ) %>%
    add_global_p(keep = TRUE, singular.ok = TRUE) %>%
    bold_labels()
  
  result <- tbl_merge(                           # merge les colonnes 
    tbls = list(effectif, univ, multi),
    tab_spanner = c("**N**", 
                    "**Adjusted on MiSeq batch effect**",
                    "**Fully adjusted**"))
  return(result)
}

### Tableaux de résultats (effectif + beta + IC95% + p)
model_multi <- function(data, outcome, exposure_vec){
  
  effectif <- data %>%                   # Création d'une colonne effectif
    filter(!is.na({{outcome}}))%>%
    select(all_of({{exposure_vec}})) %>%
    tbl_summary(missing = "no", 
                statistic = all_continuous() ~ "{N_nonmiss}") %>%
    bold_labels()
  
  multi <- data %>%                  
    select(
      {{outcome}},
      all_of({{exposure_vec}}),
      all_of(covar_vec_i)) %>% 
    
    tbl_uvregression(
      method = lm ,
      y = {{outcome}},
      formula = "{y} ~ {x} +
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
      estimate_fun = ~ style_sigfig(.x, digits = 1)
    ) %>%
    add_global_p(keep = TRUE, singular.ok = TRUE) %>%
    bold_labels() 
  
  result <- tbl_merge(
    tbls = list(effectif, multi), 
    tab_spanner = c("**N**", 
                    "**Fully adjusted**"))
  
  return(result)
}


# Statistiques (sensitivity analysis) ----
test_sensi_sg <- function(data, outcome, exposure_vec) {
  
  exposure_vec_t2 <- data %>% select(all_of({{exposure_vec}})) %>% select(contains("t2")) %>% colnames()
  exposure_vec_t3 <- data %>% select(all_of({{exposure_vec}})) %>% select(contains("t3")) %>% colnames()
  exposure_vec_M2 <- data %>% select(all_of({{exposure_vec}})) %>% select(contains("M2")) %>% colnames()
  exposure_vec_Y1 <- data %>% select(all_of({{exposure_vec}})) %>% select(contains("Y1")) %>% colnames()
  
i_cor_pool_sg_continu_t2 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_t2}}),
    all_of(covar_vec_i), 
    mo_pool_sg_T1_2) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         mo_pool_sg_T1_2 +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 


i_cor_pool_sg_ter_t2 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_t2}}),
    all_of(covar_vec_i), 
    mo_pool_sg_T1_2_ter) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         mo_pool_sg_T1_2_ter +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 


i_cor_pool_sg_continu_t3 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_t3}}),
    all_of(covar_vec_i), 
    mo_pool_sg_T3_2) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         mo_pool_sg_T3_2 +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 


i_cor_pool_sg_ter_t3 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_t3}}),
    all_of(covar_vec_i), 
    mo_pool_sg_T3_2_ter) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         mo_pool_sg_T3_2_ter +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 

i_cor_pool_sg_continu_M2 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_M2}}),
    all_of(covar_vec_i), 
    ch_pool_sg_M2_2) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         ch_pool_sg_M2_2 +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 


i_cor_pool_sg_ter_M2 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_M2}}),
    all_of(covar_vec_i), 
    ch_pool_sg_M2_2_ter) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         ch_pool_sg_M2_2_ter +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 

i_cor_pool_sg_continu_Y1 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_Y1}}),
    all_of(covar_vec_i), 
    ch_pool_sg_Y1_2) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         ch_pool_sg_Y1_2 +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 


i_cor_pool_sg_ter_Y1 <- 
  
  data %>%                  # Création des colonnes modèle complet
  select(
    {{outcome}},
    all_of({{exposure_vec_Y1}}),
    all_of(covar_vec_i), 
    ch_pool_sg_Y1_2_ter) %>% 
  
  tbl_uvregression(
    method = lm ,
    y = {{outcome}},
    formula = "{y} ~ {x} +
         ch_feces_RUN_Y1  +
         ch_pool_sg_Y1_2_ter +
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
    pvalue_fun = ~ style_pvalue(.x, digits = 3),
    estimate_fun = ~ style_sigfig(.x, digits = 2)
  ) %>%
  add_global_p(keep = TRUE, singular.ok = TRUE) %>%
  bold_labels() 

result <- tbl_merge(
  tbls = list(i_cor_pool_sg_continu_t2, i_cor_pool_sg_ter_t2, 
              i_cor_pool_sg_continu_t3, i_cor_pool_sg_ter_t3, 
              i_cor_pool_sg_continu_M2, i_cor_pool_sg_ter_M2,
              i_cor_pool_sg_continu_Y1, i_cor_pool_sg_ter_Y1), 
  tab_spanner = c("**Fully adjusted i_cor + pool_sg continu t2**", 
                  "**Fully adjusted i_cor + pool_sg tertiles t2**", 
                  "**Fully adjusted i_cor + pool_sg continu t3**", 
                  "**Fully adjusted i_cor + pool_sg tertiles t3**", 
                  "**Fully adjusted i_cor + pool_sg continu M2**", 
                  "**Fully adjusted i_cor + pool_sg tertiles M2**", 
                  "**Fully adjusted i_cor + pool_sg continu Y1**", 
                  "**Fully adjusted i_cor + pool_sg tertiles Y1**")
)
return(result)
}
