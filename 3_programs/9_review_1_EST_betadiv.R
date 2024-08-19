# 10. review_1_EST

# Packages, functions and data loading ----
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
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')
rm(model_covar, model_multi, model_summary, model_univ_multi, table_cor, table_cor_sg, test_sensi_sg)
source("3_programs/4_vectors_AD_gumme.R")
rm(phenols_vec_2, 
   pfas_vec_exclu17673, pfas_vec_2, pfas_vec_3, 
   pollutant_vec_t2, pollutant_vec_t3, pollutant_vec_M2, pollutant_vec_Y1, 
   covar_vec, covar_vec_cat, covar_vec_cat_i,
   list = ls()[grep("sg", ls())])
rm(list = ls()[grep("phthalates", ls())])
rm(list = ls()[grep("num", ls())])
asv_raw_not_rarefied <- read_labelled_csv("0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") 

# Beta diversity ----
## Data cleaning ----
### Rarefaction ----
#### Data preparation 
asv_raw_not_rarefied <- 
  asv_raw_not_rarefied %>%
  select(ident, 
         starts_with("ch_feces_raw_asv")) %>%
  na.omit()
row.names(asv_raw_not_rarefied) <- NULL 
asv_raw_not_rarefied <- asv_raw_not_rarefied %>%
  column_to_rownames("ident") %>%
  t() %>%
  as.data.frame()%>% 
  rownames_to_column("ch_feces_ASV_ID_Y1") %>%
  mutate(
    ch_feces_ASV_ID_Y1 = str_replace_all(ch_feces_ASV_ID_Y1, c("ch_feces_raw_" = "", "_Y1"=""))) %>%
  column_to_rownames("ch_feces_ASV_ID_Y1") %>%
  otu_table(taxa_are_rows = TRUE)     

#### from the raw ASV table, choose a threshold for the sequencing depth 
#### define a subset of samples to keep = samples with a sequencing depth > the chosen threshold  
#### samples with a sequencing depth < the chosen threshold become missing data 
keep_5000 <- 
  names(which(sample_sums(asv_raw_not_rarefied)>= 5000)) %>%
  prune_samples(asv_raw_not_rarefied) %>% 
  as.data.frame()                # Loss of 6 samples 

#### reduce the sequencing depth of the samples to the chosen threshold
#### the sequences kept within each sample are randomly selected 
ASV_rarefied_5000_Y1 <- keep_5000 %>%
  SRS(5000, set_seed = TRUE, seed = 1)
rownames(ASV_rarefied_5000_Y1)<- rownames(keep_5000)

#### put the dataframe with rarefied ASVs in columns and samples in rows
ASV_rarefied_5000_Y1 <- 
  ASV_rarefied_5000_Y1 %>%
  t() %>%
  as.data.frame()
rm(keep_5000)

#### check if the rarefaction worked properly 
#### ok the rowsums are equal to 5000, the threshold we chose
rowSums(ASV_rarefied_5000_Y1)   

### Metric calculation (Bray Curtis) ----
### Calcul de la metric Bray Curtis
set.seed(1996)
metric_bray_curtis <- vegan::vegdist(ASV_rarefied_5000_Y1, method = "bray")
### warning message: Plus d’une classe "dist" est trouvée en cache : Utilisation de la première, depuis l’espace de noms 'BiocGenerics'. Aussi défini par ‘spam’
### ok, phyloseq package is supposed to use BiocGeneric 
nmds <- metaMDS(metric_bray_curtis)
scores(nmds) %>%
  as_tibble(rownames = "ident") %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point()
# nmds <- metaMDS(ASV_rarefied_5000_Y1, autotransform = FALSE) # ?
# scores(nmds) %>%
#   as_tibble(rownames = "ident") %>%
#   ggplot(aes(x=NMDS1, y=NMDS2)) +
#   geom_point()

metric_bray_curtis <- 
  metric_bray_curtis %>% 
  as.matrix %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ident")

### Metadata ----
metadata <- metadata %>%
  select(ident, 
         all_of(covar_vec_i), 
         all_of(phenols_vec_ter), 
         all_of(phenols_vec_cat), 
         all_of(pfas_vec_ter), 
         all_of(pfas_vec_cat)) %>%
  mutate(ident = as.character(ident))

#### merge betadiversity data and metadata
bdd_final <- inner_join(metric_bray_curtis, metadata, by = "ident")
bdd_final[c(phenols_vec_ter, pfas_vec_ter)] <- lapply(bdd_final[c(phenols_vec_ter, pfas_vec_ter)], as.factor)
covariables <- bdd_final %>% select(all_of(covar_vec_i)) %>% select(-mo_interpreg_3cat) %>% colnames()

#### filter t2 phenols (because of the NA on the pollutants)
bdd_final %>% select(ident, contains("t2")) %>% filter_all(any_vars(is.na(.)))
bdd_final_t2 <- bdd_final %>% 
  select(ident, contains("t2"), everything()) %>%
  filter(ident != 15804) %>%
  select(-"15804")
all_dist_t2 <- bdd_final_t2 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter t2 pfas (because of the NA on the pollutants)
bdd_final %>% select(ident, 
                     "mo_PFOA_cor_ter", 
                     "mo_PFNA_ter", 
                     "mo_PFHxS_cor_ter", 
                     "mo_PFOS_cor_ter", 
                     "mo_PFDA_i_cor_ter", 
                     "mo_PFUnDA_i_cor_ter", 
                     "mo_PFHpS_i_cor_ter",
                     "mo_PFDoDa_cat", 
                     "mo_PFHxPA_cat",
                     "mo_PFHpA_cat_2", 
                     "mo_PFTrDa_cat", 
                     "mo_PFBS_cat_2", 
                     "mo_PFOSA_cat_2", 
                     "mo_6_2diPAP_cat_2", 
                     "mo_8_2diPAP_cat_2") %>% filter_all(any_vars(is.na(.)))
bdd_final_pfas_t2 <- bdd_final %>% 
  select(ident, 
         "mo_PFOA_cor_ter", 
         "mo_PFNA_ter", 
         "mo_PFHxS_cor_ter", 
         "mo_PFOS_cor_ter", 
         "mo_PFDA_i_cor_ter", 
         "mo_PFUnDA_i_cor_ter", 
         "mo_PFHpS_i_cor_ter",
         "mo_PFDoDa_cat", 
         "mo_PFHxPA_cat",
         "mo_PFHpA_cat_2", 
         "mo_PFTrDa_cat", 
         "mo_PFBS_cat_2", 
         "mo_PFOSA_cat_2", 
         "mo_6_2diPAP_cat_2", 
         "mo_8_2diPAP_cat_2", everything()) %>%
  filter(!ident %in% c("15463" , "15736", "16318", "16958", "17154", "17958", 
                       "18437", "19378", "19971", "23311", "23371", "23901", 
                       "24862", "26891", "26960", "27015")) %>%
  select(-"15463" , -"15736", -"16318", -"16958", -"17154", -"17958", -"18437", 
         -"19378", -"19971", -"23311", -"23371", -"23901", -"24862", -"26891", 
         -"26960", -"27015")
all_dist_pfas_t2 <- bdd_final_pfas_t2 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter t3 (because of the NA on the pollutants)
bdd_final %>% select(ident, contains("t3")) %>% filter_all(any_vars(is.na(.)))
bdd_final_t3 <- bdd_final %>% 
  select(ident, contains("t3"), everything()) %>%
  filter(!ident %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  select(-"17827", -"15929", -"26891", -"25668", -"23330", -"28199")
all_dist_t3 <- bdd_final_t3 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter M2 (because of the NA on the pollutants)
bdd_final %>% select(ident, contains("M2")) %>% filter_all(any_vars(is.na(.)))
bdd_final_M2 <- bdd_final %>% 
  select(ident, contains("M2"), everything()) %>%
  filter(!ident %in% c(14461, 27630, 25673, 25613, 22869, 17717, 15741, 15463, 23330, 22620, 18111, 22126)) %>%
  select(-"14461", -"27630", -"25673", -"25613", -"22869", -"17717", -"15741", -"15463", -"23330", -"22620", -"18111", -"22126")
all_dist_M2 <- bdd_final_M2 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter Y1 (because of the NA on the pollutants)
bdd_final %>% select(ident, contains("Y1")) %>% filter_all(any_vars(is.na(.)))
bdd_final_Y1 <- bdd_final %>% 
  select(ident, contains("Y1"), everything()) %>%
  filter(!ident %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  select(-"23994", -"25166", -"26766", -"14668", -"26923")
all_dist_Y1 <- bdd_final_Y1 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()



## Statistical analysis Adonis2 ----
### Trim 2. phenols ----
explanatory_vars_t2 <- 
  bdd_final_t2 %>% 
  select("mo_MEPA_total_i_cor_t2_ter", 
         "mo_ETPA_total_i_cor_t2_ter", 
         "mo_PRPA_total_i_cor_t2_ter", 
         "mo_BUPA_total_cat_t2", 
         "mo_TRCS_total_i_cor_t2_ter", 
         "mo_BPA_total_i_cor_t2_ter", 
         "mo_BPS_total_cat_t2_2", 
         "mo_OXBE_total_i_cor_t2_ter") %>%
  colnames()

results_betadiv_univar_bray_curtis_t2 <- 
  lapply(explanatory_vars_t2, function(x) {
    formula <- reformulate(x, response = "all_dist_t2")
    adonis2(formula, data = bdd_final_t2, permutations = 999)
  })
results_betadiv_univar_bray_curtis_t2 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_t2) %>%
  rownames_to_column(var = "Explanatory variables") %>%
    mutate(
      Pollutants = c(rep("mo_MEPA_total_i_cor_t2_ter", times = 3), 
                     rep("mo_ETPA_total_i_cor_t2_ter", times = 3), 
                     rep("mo_PRPA_total_i_cor_t2_ter", times = 3),
                     rep("mo_BUPA_total_cat_t2", times = 3),
                     rep("mo_TRCS_total_i_cor_t2_ter", times = 3),
                     rep("mo_BPA_total_i_cor_t2_ter", times = 3),
                     rep("mo_BPS_total_cat_t2_2", times = 3),
                     rep("mo_OXBE_total_i_cor_t2_ter", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_t2 <-
  lapply(explanatory_vars_t2, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_t2")
    adonis2(formula, data = bdd_final_t2, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_t2 <- 
  do.call(rbind, results_betadiv_multivar_bray_curtis_t2) %>%
  rownames_to_column(var = "Explanatory variables") %>% 
  mutate(
    Pollutants = c(rep("mo_MEPA_total_i_cor_t2_ter", times = 23), 
                   rep("mo_ETPA_total_i_cor_t2_ter", times = 23), 
                   rep("mo_PRPA_total_i_cor_t2_ter", times = 23),
                   rep("mo_BUPA_total_cat_t2", times = 23),
                   rep("mo_TRCS_total_i_cor_t2_ter", times = 23),
                   rep("mo_BPA_total_i_cor_t2_ter", times = 23),
                   rep("mo_BPS_total_cat_t2_2", times = 23),
                   rep("mo_OXBE_total_i_cor_t2_ter", times = 23))) %>%
  select(Pollutants, everything())

### Trim 2. pfas ----
explanatory_vars_pfas_t2 <- 
  bdd_final_pfas_t2 %>% 
  select( "mo_PFOA_cor_ter", 
          "mo_PFNA_ter", 
          "mo_PFHxS_cor_ter", 
          "mo_PFOS_cor_ter", 
          "mo_PFDA_i_cor_ter", 
          "mo_PFUnDA_i_cor_ter", 
          "mo_PFHpS_i_cor_ter",
          "mo_PFDoDa_cat", 
          "mo_PFHxPA_cat",
          "mo_PFHpA_cat_2", 
          "mo_PFTrDa_cat", 
          "mo_PFBS_cat_2", 
          "mo_PFOSA_cat_2", 
          "mo_6_2diPAP_cat_2", 
          "mo_8_2diPAP_cat_2") %>%
  colnames()

results_betadiv_univar_bray_curtis_pfas_t2 <- 
  lapply(explanatory_vars_pfas_t2, function(x) {
    formula <- reformulate(x, response = "all_dist_pfas_t2")
    adonis2(formula, data = bdd_final_pfas_t2, permutations = 999)
  })
results_betadiv_univar_bray_curtis_pfas_t2 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_pfas_t2) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("mo_PFOA_cor_ter", times = 3), 
                   rep("mo_PFNA_ter", times = 3), 
                   rep("mo_PFHxS_cor_ter", times = 3),
                   rep("mo_PFOS_cor_ter", times = 3),
                   rep("mo_PFDA_i_cor_ter", times = 3),
                   rep("mo_PFUnDA_i_cor_ter", times = 3),
                   rep("mo_PFHpS_i_cor_ter", times = 3),
                   rep("mo_PFDoDa_cat", times = 3), 
                   rep("mo_PFHxPA_cat", times = 3),
                   rep("mo_PFHpA_cat_2", times = 3),
                   rep("mo_PFTrDa_cat", times = 3),
                   rep("mo_PFBS_cat_2", times = 3),
                   rep("mo_PFOSA_cat_2", times = 3),
                   rep("mo_6_2diPAP_cat_2", times = 3),
                   rep("mo_8_2diPAP_cat_2", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_pfas_t2 <-
  lapply(explanatory_vars_pfas_t2, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_pfas_t2")
    adonis2(formula, data = bdd_final_pfas_t2, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_pfas_t2 <- 
  do.call(rbind, results_betadiv_multivar_bray_curtis_pfas_t2) %>%
  rownames_to_column(var = "Explanatory variables") %>% 
  mutate(
    Pollutants = c(rep("mo_PFOA_cor_ter", times = 23), 
                   rep("mo_PFNA_ter", times = 23), 
                   rep("mo_PFHxS_cor_ter", times = 23),
                   rep("mo_PFOS_cor_ter", times = 23),
                   rep("mo_PFDA_i_cor_ter", times = 23),
                   rep("mo_PFUnDA_i_cor_ter", times = 23),
                   rep("mo_PFHpS_i_cor_ter", times = 23),
                   rep("mo_PFDoDa_cat", times = 23), 
                   rep("mo_PFHxPA_cat", times = 23),
                   rep("mo_PFHpA_cat_2", times = 23),
                   rep("mo_PFTrDa_cat", times = 23),
                   rep("mo_PFBS_cat_2", times = 23),
                   rep("mo_PFOSA_cat_2", times = 23),
                   rep("mo_6_2diPAP_cat_2", times = 23),
                   rep("mo_8_2diPAP_cat_2", times = 23))) %>%
  select(Pollutants, everything())

# Visualisation des résultats significatifs
results_betadiv_univar_bray_curtis_pfas_t2 %>% 
  filter(`Pr(>F)` < 0.05) %>% 
  filter(`Explanatory variables` %in% c(explanatory_vars_pfas_t2)) %>% 
  View()

results_betadiv_multivar_bray_curtis_pfas_t2 %>% 
  filter(`Pr(>F)` < 0.05) %>% 
  filter(`Explanatory variables` %in% c(explanatory_vars_pfas_t2)) %>% 
  View()


### Trim 3. ----
explanatory_vars_t3 <- 
  bdd_final_t3 %>% 
  select("mo_MEPA_total_i_cor_t3_ter", 
         "mo_ETPA_total_i_cor_t3_ter", 
         "mo_PRPA_total_i_cor_t3_ter", 
         "mo_BUPA_total_cat_t3", 
         "mo_TRCS_total_i_cor_t3_ter", 
         "mo_BPA_total_i_cor_t3_ter", 
         "mo_BPS_total_cat_t3_2",
         "mo_OXBE_total_i_cor_t3_ter")  %>%
  colnames()

results_betadiv_univar_bray_curtis_t3 <- 
  lapply(explanatory_vars_t3, function(x) {
    formula <- reformulate(x, response = "all_dist_t3")
    adonis2(formula, data = bdd_final_t3, permutations = 999)
  })
results_betadiv_univar_bray_curtis_t3 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_t3) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("mo_MEPA_total_i_cor_t3_ter", times = 3), 
                   rep("mo_ETPA_total_i_cor_t3_ter", times = 3), 
                   rep("mo_PRPA_total_i_cor_t3_ter", times = 3),
                   rep("mo_BUPA_total_cat_t3", times = 3),
                   rep("mo_TRCS_total_i_cor_t3_ter", times = 3),
                   rep("mo_BPA_total_i_cor_t3_ter", times = 3),
                   rep("mo_BPS_total_cat_t3_2", times = 3),
                   rep("mo_OXBE_total_i_cor_t3_ter", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_t3 <-
  lapply(explanatory_vars_t3, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_t3")
    adonis2(formula, data = bdd_final_t3, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_t3 <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_t3) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(
    rep("mo_MEPA_total_i_cor_t3_ter", times = 23),
    rep("mo_ETPA_total_i_cor_t3_ter", times = 23),
    rep("mo_PRPA_total_i_cor_t3_ter", times = 23),
    rep("mo_BUPA_total_cat_t3", times = 23),
    rep("mo_TRCS_total_i_cor_t3_ter", times = 23),
    rep("mo_BPA_total_i_cor_t3_ter", times = 23),
    rep("mo_BPS_total_cat_t3_2", times = 23),
    rep("mo_OXBE_total_i_cor_t3_ter", times = 23)
  )) %>%
  select(Pollutants, everything())


# Visualisation des résultats significatifs
results_betadiv_univar_bray_curtis_t3 %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t3)) %>%
  View()

results_betadiv_multivar_bray_curtis_t3 %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t3)) %>%
  View()


### M2 ----
explanatory_vars_M2 <- 
  bdd_final_M2 %>% 
  select("ch_MEPA_total_i_cor_M2_ter", 
         "ch_ETPA_total_cat_M2_2", 
         "ch_PRPA_total_cat_M2_2",
         "ch_BUPA_total_cat_M2_2",
         "ch_TRCS_total_i_cor_M2_ter",
         "ch_BPA_total_i_cor_M2_ter", 
         "ch_BPS_total_cat_M2_2", 
         "ch_OXBE_total_i_cor_M2_ter")  %>%
  colnames()

results_betadiv_univar_bray_curtis_M2 <- 
  lapply(explanatory_vars_M2, function(x) {
    formula <- reformulate(x, response = "all_dist_M2")
    adonis2(formula, data = bdd_final_M2, permutations = 999)
  })
results_betadiv_univar_bray_curtis_M2 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_M2) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("ch_MEPA_total_i_cor_M2_ter", times = 3), 
                   rep("ch_ETPA_total_cat_M2_2", times = 3), 
                   rep("ch_PRPA_total_cat_M2_2", times = 3),
                   rep("ch_BUPA_total_cat_M2_2", times = 3),
                   rep("ch_TRCS_total_i_cor_M2_ter", times = 3),
                   rep("ch_BPA_total_i_cor_M2_ter", times = 3),
                   rep("ch_BPS_total_cat_M2_2", times = 3),
                   rep("ch_OXBE_total_i_cor_M2_ter", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_M2 <-
  lapply(explanatory_vars_M2, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_M2")
    adonis2(formula, data = bdd_final_M2, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_M2 <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_M2) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(
    rep("ch_MEPA_total_i_cor_M2_ter", times = 23),
    rep("ch_ETPA_total_cat_M2_2", times = 23),
    rep("ch_PRPA_total_cat_M2_2", times = 23),
    rep("ch_BUPA_total_cat_M2_2", times = 23),
    rep("ch_TRCS_total_i_cor_M2_ter", times = 23),
    rep("ch_BPA_total_i_cor_M2_ter", times = 23),
    rep("ch_BPS_total_cat_M2_2", times = 23),
    rep("ch_OXBE_total_i_cor_M2_ter", times = 23)
  )) %>%
  select(Pollutants, everything())


# Visualisation des résultats significatifs
results_betadiv_univar_bray_curtis_M2 %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_M2)) %>%
  View()

results_betadiv_multivar_bray_curtis_M2 %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_M2)) %>%
  View()


### Y1 ----
explanatory_vars_Y1 <- 
  bdd_final_Y1 %>% 
  select("ch_MEPA_total_i_cor_Y1_ter",
         "ch_ETPA_total_i_cor_Y1_ter",
         "ch_PRPA_total_i_cor_Y1_ter",
         "ch_BUPA_total_cat_Y1",      
         "ch_TRCS_total_i_cor_Y1_ter",
         "ch_BPA_total_i_cor_Y1_ter",
         "ch_BPS_total_cat_Y1",
         "ch_OXBE_total_i_cor_Y1_ter")  %>%
  colnames()

results_betadiv_univar_bray_curtis_Y1 <- 
  lapply(explanatory_vars_Y1, function(x) {
    formula <- reformulate(x, response = "all_dist_Y1")
    adonis2(formula, data = bdd_final_Y1, permutations = 999)
  })
results_betadiv_univar_bray_curtis_Y1 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_Y1) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("ch_MEPA_total_i_cor_Y1_ter", times = 3), 
                   rep("ch_ETPA_total_i_cor_Y1_ter", times = 3), 
                   rep("ch_PRPA_total_i_cor_Y1_ter", times = 3),
                   rep("ch_BUPA_total_cat_Y1", times = 3),
                   rep("ch_TRCS_total_i_cor_Y1_ter", times = 3),
                   rep("ch_BPA_total_i_cor_Y1_ter", times = 3),
                   rep("ch_BPS_total_cat_Y1", times = 3),
                   rep("ch_OXBE_total_i_cor_Y1_ter", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_Y1 <-
  lapply(explanatory_vars_Y1, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_Y1")
    adonis2(formula, data = bdd_final_Y1, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_Y1 <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_Y1) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(
    rep("ch_MEPA_total_i_cor_Y1_ter", times = 23),
    rep("ch_ETPA_total_i_cor_Y1_ter", times = 23),
    rep("ch_PRPA_total_i_cor_Y1_ter", times = 23),
    rep("ch_BUPA_total_cat_Y1", times = 23),
    rep("ch_TRCS_total_i_cor_Y1_ter", times = 23),
    rep("ch_BPA_total_i_cor_Y1_ter", times = 23),
    rep("ch_BPS_total_cat_Y1", times = 23),
    rep("ch_OXBE_total_i_cor_Y1_ter", times = 23)
  )) %>%
  select(Pollutants, everything())


# Visualisation des résultats significatifs
results_betadiv_univar_bray_curtis_Y1 %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_Y1)) %>%
  View()

results_betadiv_multivar_bray_curtis_Y1 %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_Y1)) %>%
  View()



### sous groupe OXBE M2 ----
#### 1st vs 2nd tertile ----
bdd_oxbe_M2_1_2 <- bdd_final_M2 %>% 
  filter(ch_OXBE_total_i_cor_M2_ter %in% c("1st tertile", "2nd tertile"))
oxbe_to_keep_1_2 <- bdd_oxbe_M2_1_2$ident
bdd_oxbe_M2_1_2 <- bdd_oxbe_M2_1_2 %>%
  select(ident, all_of(covariables), ch_OXBE_total_i_cor_M2_ter, all_of(oxbe_to_keep_1_2))
all_dist_oxbe_M2_1_2 <- bdd_oxbe_M2_1_2 %>% select(all_of(oxbe_to_keep_1_2)) %>% as.dist()

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

results_betadiv_multivar_bray_curtis_oxbe_M2_1st_2nd_ter <- 
  results_betadiv_multivar_bray_curtis_oxbe_M2_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_OXBE_total_i_cor_M2_ter", times = 23)), 
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

#### 1st vs 3rd tertile ----
bdd_oxbe_M2_1_3 <- bdd_final_M2 %>% 
  filter(ch_OXBE_total_i_cor_M2_ter %in% c("1st tertile", "3rd tertile"))
oxbe_to_keep_1_3 <- bdd_oxbe_M2_1_3$ident
bdd_oxbe_M2_1_3 <- bdd_oxbe_M2_1_3 %>%
  select(ident, all_of(covariables), ch_OXBE_total_i_cor_M2_ter, all_of(oxbe_to_keep_1_3))
all_dist_oxbe_M2_1_3 <- bdd_oxbe_M2_1_3 %>% select(all_of(oxbe_to_keep_1_3)) %>% as.dist()

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

results_betadiv_multivar_bray_curtis_oxbe_M2_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_oxbe_M2_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_OXBE_total_i_cor_M2_ter", times = 23)), 
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



#### 2nd vs 3rd tertile ----
bdd_oxbe_M2_2_3 <- bdd_final_M2 %>% 
  filter(ch_OXBE_total_i_cor_M2_ter %in% c("2nd tertile", "3rd tertile"))
oxbe_to_keep_2_3 <- bdd_oxbe_M2_2_3$ident
bdd_oxbe_M2_2_3 <- bdd_oxbe_M2_2_3 %>%
  select(ident, all_of(covariables), ch_OXBE_total_i_cor_M2_ter, all_of(oxbe_to_keep_2_3))
all_dist_oxbe_M2_2_3 <- bdd_oxbe_M2_2_3 %>% select(all_of(oxbe_to_keep_2_3)) %>% as.dist()

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

results_betadiv_multivar_bray_curtis_oxbe_M2_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_oxbe_M2_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_OXBE_total_i_cor_M2_ter", times = 23)), 
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_oxbe_M2_1_2, bdd_oxbe_M2_1_3, bdd_oxbe_M2_2_3, 
   oxbe_to_keep_1_2, oxbe_to_keep_1_3, oxbe_to_keep_2_3, 
   all_dist_oxbe_M2_1_2, all_dist_oxbe_M2_1_3, all_dist_oxbe_M2_2_3)


results_betadiv_multivar_bray_curtis_oxbe_M2 <-
  bind_rows(results_betadiv_multivar_bray_curtis_oxbe_M2_1st_2nd_ter, 
            results_betadiv_multivar_bray_curtis_oxbe_M2_1st_3rd_ter, 
            results_betadiv_multivar_bray_curtis_oxbe_M2_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("ch_OXBE_total_i_cor_M2_ter", "Residual", "Total"))



### sous groupe MEPA Y1 ----
#### 1st vs 2nd tertile ----
bdd_mepa_Y1_1_2 <- bdd_final_Y1 %>% 
  filter(ch_MEPA_total_i_cor_Y1_ter %in% c("1st tertile", "2nd tertile"))
mepa_to_keep_1_2 <- bdd_mepa_Y1_1_2$ident
bdd_mepa_Y1_1_2 <- bdd_mepa_Y1_1_2 %>%
  select(ident, all_of(covariables), ch_MEPA_total_i_cor_Y1_ter, all_of(mepa_to_keep_1_2))
all_dist_mepa_Y1_1_2 <- bdd_mepa_Y1_1_2 %>% select(all_of(mepa_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mepa_Y1_1st_2nd_ter <- 
  adonis2(all_dist_mepa_Y1_1_2 ~ 
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
          data = bdd_mepa_Y1_1_2, 
          permutations = 999)

results_betadiv_multivar_bray_curtis_mepa_Y1_1st_2nd_ter <- 
  results_betadiv_multivar_bray_curtis_mepa_Y1_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_MEPA_total_i_cor_Y1_ter", times = 23)), 
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

#### 1st vs 3rd tertile ----
bdd_mepa_Y1_1_3 <- bdd_final_Y1 %>% 
  filter(ch_MEPA_total_i_cor_Y1_ter %in% c("1st tertile", "3rd tertile"))
mepa_to_keep_1_3 <- bdd_mepa_Y1_1_3$ident
bdd_mepa_Y1_1_3 <- bdd_mepa_Y1_1_3 %>%
  select(ident, all_of(covariables), ch_MEPA_total_i_cor_Y1_ter, all_of(mepa_to_keep_1_3))
all_dist_mepa_Y1_1_3 <- bdd_mepa_Y1_1_3 %>% select(all_of(mepa_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mepa_Y1_1st_3rd_ter <- 
  adonis2(all_dist_mepa_Y1_1_3 ~ 
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
          data = bdd_mepa_Y1_1_3, 
          permutations = 999)

results_betadiv_multivar_bray_curtis_mepa_Y1_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_mepa_Y1_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_MEPA_total_i_cor_Y1_ter", times = 23)), 
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

#### 2nd vs 3rd tertile ----
bdd_mepa_Y1_2_3 <- bdd_final_Y1 %>% 
  filter(ch_MEPA_total_i_cor_Y1_ter %in% c("2nd tertile", "3rd tertile"))
mepa_to_keep_2_3 <- bdd_mepa_Y1_2_3$ident
bdd_mepa_Y1_2_3 <- bdd_mepa_Y1_2_3 %>%
  select(ident, all_of(covariables), ch_MEPA_total_i_cor_Y1_ter, all_of(mepa_to_keep_2_3))
all_dist_mepa_Y1_2_3 <- bdd_mepa_Y1_2_3 %>% select(all_of(mepa_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mepa_Y1_2nd_3rd_ter <- 
  adonis2(all_dist_mepa_Y1_2_3 ~ 
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
          data = bdd_mepa_Y1_2_3, 
          permutations = 999)

results_betadiv_multivar_bray_curtis_mepa_Y1_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_mepa_Y1_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_MEPA_total_i_cor_Y1_ter", times = 23)), 
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_mepa_Y1_1_2, bdd_mepa_Y1_1_3, bdd_mepa_Y1_2_3, 
   mepa_to_keep_1_2, mepa_to_keep_1_3, mepa_to_keep_2_3, 
   all_dist_mepa_Y1_1_2, all_dist_mepa_Y1_1_3, all_dist_mepa_Y1_2_3)


results_betadiv_multivar_bray_curtis_mepa_Y1 <-
  bind_rows(results_betadiv_multivar_bray_curtis_mepa_Y1_1st_2nd_ter, 
            results_betadiv_multivar_bray_curtis_mepa_Y1_1st_3rd_ter, 
            results_betadiv_multivar_bray_curtis_mepa_Y1_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("ch_MEPA_total_i_cor_Y1_ter", "Residual", "Total"))


### sous groupe PRPA Y1 ----
#### 1st vs 2nd tertile ----
bdd_prpa_Y1_1_2 <- bdd_final_Y1 %>% 
  filter(ch_PRPA_total_i_cor_Y1_ter %in% c("1st tertile", "2nd tertile"))
prpa_to_keep_1_2 <- bdd_prpa_Y1_1_2$ident
bdd_prpa_Y1_1_2 <- bdd_prpa_Y1_1_2 %>%
  select(ident, all_of(covariables), ch_PRPA_total_i_cor_Y1_ter, all_of(prpa_to_keep_1_2))
all_dist_prpa_Y1_1_2 <- bdd_prpa_Y1_1_2 %>% select(all_of(prpa_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_prpa_Y1_1st_2nd_ter <- 
  adonis2(all_dist_prpa_Y1_1_2 ~ 
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
          data = bdd_prpa_Y1_1_2, 
          permutations = 999)

results_betadiv_multivar_bray_curtis_prpa_Y1_1st_2nd_ter <- 
  results_betadiv_multivar_bray_curtis_prpa_Y1_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_PRPA_total_i_cor_Y1_ter", times = 23)), 
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

#### 1st vs 3rd tertile ----
bdd_prpa_Y1_1_3 <- bdd_final_Y1 %>% 
  filter(ch_PRPA_total_i_cor_Y1_ter %in% c("1st tertile", "3rd tertile"))
prpa_to_keep_1_3 <- bdd_prpa_Y1_1_3$ident
bdd_prpa_Y1_1_3 <- bdd_prpa_Y1_1_3 %>%
  select(ident, all_of(covariables), ch_PRPA_total_i_cor_Y1_ter, all_of(prpa_to_keep_1_3))
all_dist_prpa_Y1_1_3 <- bdd_prpa_Y1_1_3 %>% select(all_of(prpa_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_prpa_Y1_1st_3rd_ter <- 
  adonis2(all_dist_prpa_Y1_1_3 ~ 
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
          data = bdd_prpa_Y1_1_3, 
          permutations = 999)

results_betadiv_multivar_bray_curtis_prpa_Y1_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_prpa_Y1_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_PRPA_total_i_cor_Y1_ter", times = 23)), 
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



#### 2nd vs 3rd tertile ----
bdd_prpa_Y1_2_3 <- bdd_final_Y1 %>% 
  filter(ch_PRPA_total_i_cor_Y1_ter %in% c("2nd tertile", "3rd tertile"))
prpa_to_keep_2_3 <- bdd_prpa_Y1_2_3$ident
bdd_prpa_Y1_2_3 <- bdd_prpa_Y1_2_3 %>%
  select(ident, all_of(covariables), ch_PRPA_total_i_cor_Y1_ter, all_of(prpa_to_keep_2_3))
all_dist_prpa_Y1_2_3 <- bdd_prpa_Y1_2_3 %>% select(all_of(prpa_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_prpa_Y1_2nd_3rd_ter <- 
  adonis2(all_dist_prpa_Y1_2_3 ~ 
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
          data = bdd_prpa_Y1_2_3, 
          permutations = 999)

results_betadiv_multivar_bray_curtis_prpa_Y1_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_prpa_Y1_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_PRPA_total_i_cor_Y1_ter", times = 23)), 
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_prpa_Y1_1_2, bdd_prpa_Y1_1_3, bdd_prpa_Y1_2_3, 
   prpa_to_keep_1_2, prpa_to_keep_1_3, prpa_to_keep_2_3, 
   all_dist_prpa_Y1_1_2, all_dist_prpa_Y1_1_3, all_dist_prpa_Y1_2_3)


results_betadiv_multivar_bray_curtis_prpa_Y1 <-
  bind_rows(results_betadiv_multivar_bray_curtis_prpa_Y1_1st_2nd_ter, 
            results_betadiv_multivar_bray_curtis_prpa_Y1_1st_3rd_ter, 
            results_betadiv_multivar_bray_curtis_prpa_Y1_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("ch_PRPA_total_i_cor_Y1_ter", "Residual", "Total"))


### Assemblage ----
results_betadiv_univ <-
  list(
    results_betadiv_univar_bray_curtis_t2,
    results_betadiv_univar_bray_curtis_pfas_t2,
    results_betadiv_univar_bray_curtis_t3,
    results_betadiv_univar_bray_curtis_M2,
    results_betadiv_univar_bray_curtis_Y1)
results_betadiv_univ <- do.call(rbind, results_betadiv_univ, quote = FALSE)
results_betadiv_univ <- results_betadiv_univ %>%
  mutate(`Explanatory variables` = if_else(`Explanatory variables` %in% c(explanatory_vars_t2,
                                                                          explanatory_vars_pfas_t2,
                                                                          explanatory_vars_t3,
                                                                          explanatory_vars_M2,
                                                                          explanatory_vars_Y1),
                                           `Explanatory variables`,
                                           str_remove(`Explanatory variables`, "\\d+$")))

results_betadiv_multi <-
  list(
    results_betadiv_multivar_bray_curtis_t2,
    results_betadiv_multivar_bray_curtis_pfas_t2,
    results_betadiv_multivar_bray_curtis_t3,
    results_betadiv_multivar_bray_curtis_M2,
    results_betadiv_multivar_bray_curtis_Y1)
results_betadiv_multi <- do.call(rbind, results_betadiv_multi, quote = FALSE)
results_betadiv_multi <-
  results_betadiv_multi %>%
  mutate(
    `Explanatory variables` = if_else(`Explanatory variables` %in% c(explanatory_vars_t2,
                                                                     explanatory_vars_pfas_t2,
                                                                     explanatory_vars_t3,
                                                                     explanatory_vars_M2,
                                                                     explanatory_vars_Y1),
                                      `Explanatory variables`,
                                      str_remove(`Explanatory variables`, "\\d+$"))) %>%
  filter(
    `Explanatory variables` %in% c(
      explanatory_vars_t2,
      explanatory_vars_pfas_t2,
      explanatory_vars_t3,
      explanatory_vars_M2,
      explanatory_vars_Y1) |
      str_detect(`Explanatory variables`, "Residual") |
      str_detect(`Explanatory variables`, "Total"))

results_betadiv_detailled <- 
  bind_rows(results_betadiv_multivar_bray_curtis_oxbe_M2, 
            results_betadiv_multivar_bray_curtis_mepa_Y1, 
            results_betadiv_multivar_bray_curtis_prpa_Y1)

results_betadiv_complet <-
  list(
    univar = list(results_betadiv_univar_bray_curtis_t2 = results_betadiv_univar_bray_curtis_t2,
                  results_betadiv_univar_bray_curtis_pfas_t2 = results_betadiv_univar_bray_curtis_pfas_t2,
                  results_betadiv_univar_bray_curtis_t3 = results_betadiv_univar_bray_curtis_t3,
                  results_betadiv_univar_bray_curtis_M2 = results_betadiv_univar_bray_curtis_M2,
                  results_betadiv_univar_bray_curtis_Y1 = results_betadiv_univar_bray_curtis_Y1),
    multivar = list(results_betadiv_multivar_bray_curtis_t2 = results_betadiv_multivar_bray_curtis_t2,
                    results_betadiv_multivar_bray_curtis_pfas_t2 = results_betadiv_multivar_bray_curtis_pfas_t2,
                    results_betadiv_multivar_bray_curtis_t3 = results_betadiv_multivar_bray_curtis_t3,
                    results_betadiv_multivar_bray_curtis_M2 = results_betadiv_multivar_bray_curtis_M2,
                    results_betadiv_multivar_bray_curtis_Y1 = results_betadiv_multivar_bray_curtis_Y1), 
    detailled = list(
      oxbe_M2 = list(results_betadiv_multivar_bray_curtis_oxbe_M2_1st_2nd_ter = 
                       results_betadiv_multivar_bray_curtis_oxbe_M2_1st_2nd_ter, 
                     results_betadiv_multivar_bray_curtis_oxbe_M2_1st_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_oxbe_M2_1st_3rd_ter, 
                     results_betadiv_multivar_bray_curtis_oxbe_M2_2nd_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_oxbe_M2_2nd_3rd_ter), 
      mepa_Y1 = list(results_betadiv_multivar_bray_curtis_mepa_Y1_1st_2nd_ter = 
                       results_betadiv_multivar_bray_curtis_mepa_Y1_1st_2nd_ter, 
                     results_betadiv_multivar_bray_curtis_mepa_Y1_1st_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_mepa_Y1_1st_3rd_ter, 
                     results_betadiv_multivar_bray_curtis_mepa_Y1_2nd_3rd_ter =
                       results_betadiv_multivar_bray_curtis_mepa_Y1_2nd_3rd_ter), 
      prpa_Y1 = list(results_betadiv_multivar_bray_curtis_prpa_Y1_1st_2nd_ter = 
                       results_betadiv_multivar_bray_curtis_prpa_Y1_1st_2nd_ter, 
                     results_betadiv_multivar_bray_curtis_prpa_Y1_1st_3rd_ter =
                       results_betadiv_multivar_bray_curtis_prpa_Y1_1st_3rd_ter, 
                     results_betadiv_multivar_bray_curtis_prpa_Y1_2nd_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_prpa_Y1_2nd_3rd_ter)))

list <- list(results_betadiv_multi, 
             results_betadiv_univ, 
             results_betadiv_multivar_bray_curtis_oxbe_M2, 
             results_betadiv_multivar_bray_curtis_mepa_Y1, 
             results_betadiv_multivar_bray_curtis_prpa_Y1)
write_xlsx(list, path = "4_output/betadiv/results_betadiv_BC.xlsx")

rm(list, 
   
   results_betadiv_univar_bray_curtis_t2, 
   results_betadiv_univar_bray_curtis_pfas_t2, 
   results_betadiv_univar_bray_curtis_t3, 
   results_betadiv_univar_bray_curtis_M2, 
   results_betadiv_univar_bray_curtis_Y1, 
   
   results_betadiv_multivar_bray_curtis_t2, 
   results_betadiv_multivar_bray_curtis_pfas_t2, 
   results_betadiv_multivar_bray_curtis_t3, 
   results_betadiv_multivar_bray_curtis_M2, 
   results_betadiv_multivar_bray_curtis_Y1, 
   
   results_betadiv_multivar_bray_curtis_oxbe_M2_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_oxbe_M2_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_oxbe_M2_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_mepa_Y1_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_mepa_Y1_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_mepa_Y1_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_prpa_Y1_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_prpa_Y1_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_prpa_Y1_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_oxbe_M2, 
   results_betadiv_multivar_bray_curtis_mepa_Y1, 
   results_betadiv_multivar_bray_curtis_prpa_Y1)


results_betadiv_univ %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t2, explanatory_vars_t3, explanatory_vars_M2, explanatory_vars_Y1)) %>%
  View()

results_betadiv_multi %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t2, explanatory_vars_t3, explanatory_vars_M2, explanatory_vars_Y1)) %>%
  View()

## Boxplots ----
### Trim2. Phenols ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t2 aux données de dissimilarité
bdd_long_t2 <- bdd_final_t2 %>% select(ident, all_of(explanatory_vars_t2))
bdd_long_t2 <- full_join(bdd_long_t2, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_t2[,explanatory_vars_t2] <- lapply(bdd_long_t2[,explanatory_vars_t2], as.character)
bdd_long_t2 <- bdd_long_t2 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_t2)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_t2_bis <- bdd_final_t2 %>%
  select(ident, 
         all_of(explanatory_vars_t2)) %>%
  rename(ident_bis = ident, 
         mo_MEPA_total_i_cor_t2_ter_bis = mo_MEPA_total_i_cor_t2_ter, 
         mo_ETPA_total_i_cor_t2_ter_bis = mo_ETPA_total_i_cor_t2_ter, 
         mo_PRPA_total_i_cor_t2_ter_bis = mo_PRPA_total_i_cor_t2_ter, 
         mo_BUPA_total_cat_t2_bis = mo_BUPA_total_cat_t2, 
         mo_TRCS_total_i_cor_t2_ter_bis = mo_TRCS_total_i_cor_t2_ter, 
         mo_BPA_total_i_cor_t2_ter_bis = mo_BPA_total_i_cor_t2_ter, 
         mo_BPS_total_cat_t2_2_bis = mo_BPS_total_cat_t2_2, 
         mo_OXBE_total_i_cor_t2_ter_bis = mo_OXBE_total_i_cor_t2_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_t2 <- full_join(bdd_long_t2, bdd_long_t2_bis, by = "ident_bis") 
bdd_long_t2 <- bdd_long_t2 %>%
  filter(ident != 15804) %>%
  filter(ident_bis != 15804)
rm(bdd_long_t2_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_t2 <- bdd_long_t2 %>%
  mutate(
    MEPA = case_when(mo_MEPA_total_i_cor_t2_ter == "1st tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MEPA_total_i_cor_t2_ter == "2nd tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MEPA_total_i_cor_t2_ter == "3rd tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MEPA_total_i_cor_t2_ter == "1st tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MEPA_total_i_cor_t2_ter == "2nd tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MEPA_total_i_cor_t2_ter == "1st tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MEPA_total_i_cor_t2_ter == "3rd tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MEPA_total_i_cor_t2_ter == "2nd tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MEPA_total_i_cor_t2_ter == "3rd tertile" & mo_MEPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ETPA = case_when(mo_ETPA_total_i_cor_t2_ter == "1st tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_ETPA_total_i_cor_t2_ter == "2nd tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_ETPA_total_i_cor_t2_ter == "3rd tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_ETPA_total_i_cor_t2_ter == "1st tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_ETPA_total_i_cor_t2_ter == "2nd tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_ETPA_total_i_cor_t2_ter == "1st tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_ETPA_total_i_cor_t2_ter == "3rd tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_ETPA_total_i_cor_t2_ter == "2nd tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_ETPA_total_i_cor_t2_ter == "3rd tertile" & mo_ETPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    PRPA = case_when(mo_PRPA_total_i_cor_t2_ter == "1st tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_PRPA_total_i_cor_t2_ter == "2nd tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_PRPA_total_i_cor_t2_ter == "3rd tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_PRPA_total_i_cor_t2_ter == "1st tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_PRPA_total_i_cor_t2_ter == "2nd tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_PRPA_total_i_cor_t2_ter == "1st tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_PRPA_total_i_cor_t2_ter == "3rd tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_PRPA_total_i_cor_t2_ter == "2nd tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_PRPA_total_i_cor_t2_ter == "3rd tertile" & mo_PRPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BUPA = case_when(mo_BUPA_total_cat_t2 == "<LOD" & mo_BUPA_total_cat_t2_bis == "<LOD" ~ "<LOD", 
                     mo_BUPA_total_cat_t2 == "LOD-LOQ" & mo_BUPA_total_cat_t2_bis == "LOD-LOQ" ~ "LOD-LOQ", 
                     mo_BUPA_total_cat_t2 == ">LOQ" & mo_BUPA_total_cat_t2_bis == ">LOQ" ~ ">LOQ", 
                     
                     mo_BUPA_total_cat_t2 == "<LOD" & mo_BUPA_total_cat_t2_bis == "LOD-LOQ" ~ "<LOD - LOD-LOQ", 
                     mo_BUPA_total_cat_t2 == "LOD-LOQ" & mo_BUPA_total_cat_t2_bis == "<LOD" ~ "<LOD - LOD-LOQ", 
                     
                     mo_BUPA_total_cat_t2 == "<LOD" & mo_BUPA_total_cat_t2_bis == ">LOQ" ~ "<LOD - >LOQ", 
                     mo_BUPA_total_cat_t2 == ">LOQ" & mo_BUPA_total_cat_t2_bis == "<LOD" ~ "<LOD - >LOQ", 
                     
                     mo_BUPA_total_cat_t2 == "LOD-LOQ" & mo_BUPA_total_cat_t2_bis == ">LOQ" ~ "LOD-LOQ - >LOQ", 
                     mo_BUPA_total_cat_t2 == ">LOQ" & mo_BUPA_total_cat_t2_bis == "LOD-LOQ" ~ "LOD-LOQ - >LOQ"), 
    
    TRCS = case_when(mo_TRCS_total_i_cor_t2_ter == "1st tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_TRCS_total_i_cor_t2_ter == "2nd tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_TRCS_total_i_cor_t2_ter == "3rd tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_TRCS_total_i_cor_t2_ter == "1st tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_TRCS_total_i_cor_t2_ter == "2nd tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_TRCS_total_i_cor_t2_ter == "1st tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_TRCS_total_i_cor_t2_ter == "3rd tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_TRCS_total_i_cor_t2_ter == "2nd tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_TRCS_total_i_cor_t2_ter == "3rd tertile" & mo_TRCS_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPA = case_when(mo_BPA_total_i_cor_t2_ter == "1st tertile" & mo_BPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                    mo_BPA_total_i_cor_t2_ter == "2nd tertile" & mo_BPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    mo_BPA_total_i_cor_t2_ter == "3rd tertile" & mo_BPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    mo_BPA_total_i_cor_t2_ter == "1st tertile" & mo_BPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    mo_BPA_total_i_cor_t2_ter == "2nd tertile" & mo_BPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    mo_BPA_total_i_cor_t2_ter == "1st tertile" & mo_BPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    mo_BPA_total_i_cor_t2_ter == "3rd tertile" & mo_BPA_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    mo_BPA_total_i_cor_t2_ter == "2nd tertile" & mo_BPA_total_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    mo_BPA_total_i_cor_t2_ter == "3rd tertile" & mo_BPA_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPS = case_when(mo_BPS_total_cat_t2_2 == "<LOD" & mo_BPS_total_cat_t2_2_bis == "<LOD" ~ "<LOD", 
                    mo_BPS_total_cat_t2_2 == ">LOD" & mo_BPS_total_cat_t2_2_bis == ">LOD" ~ ">LOD",  
                    
                    mo_BPS_total_cat_t2_2 == "<LOD" & mo_BPS_total_cat_t2_2_bis == ">LOD" ~ "<LOD - >LOD", 
                    mo_BPS_total_cat_t2_2 == ">LOD" & mo_BPS_total_cat_t2_2_bis == "<LOD" ~ "<LOD - >LOD"), 
    
    OXBE = case_when(mo_OXBE_total_i_cor_t2_ter == "1st tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_OXBE_total_i_cor_t2_ter == "2nd tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_OXBE_total_i_cor_t2_ter == "3rd tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_OXBE_total_i_cor_t2_ter == "1st tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_OXBE_total_i_cor_t2_ter == "2nd tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_OXBE_total_i_cor_t2_ter == "1st tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_OXBE_total_i_cor_t2_ter == "3rd tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_OXBE_total_i_cor_t2_ter == "2nd tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_OXBE_total_i_cor_t2_ter == "3rd tertile" & mo_OXBE_total_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         MEPA, ETPA, PRPA, BUPA, TRCS, BPA, BPS, OXBE)



# on fait passer les données en long 
bdd_long_t2 <- bdd_long_t2 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile, <LOD, LOD-LOQ, >LOQ
test <- 
  bdd_long_t2 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile", "<LOD", "LOD-LOQ", ">LOQ")) %>%  
  filter(!Pollutant == "BPS") %>%
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile", 
                             " <LOD" = "<LOD", 
                             " LOD-LOQ" = "LOD-LOQ", 
                             " >LOQ" = ">LOQ"))
bdd_long_t2 <- bind_rows(bdd_long_t2, test)
rm(test)

bdd_long_t2 <- bdd_long_t2 %>%
  mutate(
    Pollutant = fct_recode(Pollutant,
                           "Methylparaben trim.2" = "MEPA",
                           "Ethylparaben trim.2" = "ETPA",
                           "Propylparaben trim.2" = "PRPA",
                           "Butylparaben trim.2" = "BUPA",
                           "Triclosan trim.2" = "TRCS",
                           "Bisphenol A trim.2" = "BPA",
                           "Bisphenol S trim.2" = "BPS",
                           "Benzophenone 3 trim.2" = "OXBE"),
    Pollutant = fct_relevel(Pollutant, 
                            "Methylparaben trim.2", "Ethylparaben trim.2", "Propylparaben trim.2", "Butylparaben trim.2", 
                            "Triclosan trim.2", "Bisphenol A trim.2", "Bisphenol S trim.2", "Benzophenone 3 trim.2"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile", " >LOQ", ">LOD",                             # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile", "LOD-LOQ - >LOQ",             # inter groupe moyenne et forte exposition
                         " 2nd tertile"," LOD-LOQ",                                   # intra groupe moyenne exposition 
                         
                         "3rd tertile", ">LOQ",                                    # intra groupe forte exposition 
                         "<LOD - >LOD", "<LOD - >LOQ", "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile", " <LOD",
                         
                         "2nd tertile","LOD-LOQ",                                   # intra groupe moyenne exposition 
                         "<LOD - LOD-LOQ", "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile", "<LOD"),                               # intra groupe faible exposition
    Groups_rec = fct_recode(Groups, 
                            "Medium vs High" = " 3rd tertile",
                            "Medium vs High" = " >LOQ",
                            "Low vs Medium" = ">LOD",
                            "Medium vs High" = "2nd tertile - 3rd tertile",
                            "Medium vs High" = "LOD-LOQ - >LOQ",
                            "Medium vs High" = " 2nd tertile",
                            "Medium vs High" = " LOD-LOQ",
                            "Low vs High" = "3rd tertile",
                            "Low vs High" = ">LOQ",
                            "Low vs Medium" = "<LOD - >LOD",
                            "Low vs High" = "<LOD - >LOQ",
                            "Low vs High" = "1st tertile - 3rd tertile",
                            "Low vs High" = " 1st tertile",
                            "Low vs Medium" = " <LOD",
                            "Low vs Medium" = "2nd tertile",
                            "Low vs Medium" = "LOD-LOQ",
                            "Low vs Medium" = "<LOD - LOD-LOQ",
                            "Low vs Medium" = "1st tertile - 2nd tertile",
                            "Low vs Medium" = "1st tertile",
                            "Low vs Medium" = "<LOD"))



#### boxplot version 1 ----
boxplot_t2_phenols <- bdd_long_t2 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `<LOD` = "#FFE0DE",
               ` 1st tertile` = "#FFE0DE", 
               ` <LOD` = "#FFE0DE",
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `<LOD - >LOD` = "gray70",
               `<LOD - >LOQ` = "gray70",
               `<LOD - LOD-LOQ` = "gray70",
               `LOD-LOQ - >LOQ`= "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `LOD-LOQ` = "#FF8F87",
               ` 2nd tertile` = "#FF8F87",   
               ` LOD-LOQ` = "#FF8F87",
               
               `3rd tertile` = "#FF4034",   
               `>LOD` = "#FF4034", 
               `>LOQ` = "#FF4034", 
               ` 3rd tertile` = "#FF4034", 
               ` >LOQ` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(vars(Pollutant), scales = "free", ncol = 2)

# #### boxplot version 2 
# bdd_long_t2 %>%
#   ggplot() +
#   aes(
#     x = Bray_curtis_dissimilarity,
#     y = Groups,
#     fill = Groups
#   ) +
#   geom_boxplot(position = position_nudge(y = case_when(bdd_long_t2$Groups == "1st tertile" ~ -0.2, 
#                                                        bdd_long_t2$Groups == " 1st tertile" ~ -0.2, 
#                                                        bdd_long_t2$Groups == " 2nd tertile" ~ -0.2, 
#                                                        bdd_long_t2$Groups == "<LOD" ~ -0.2, 
#                                                        bdd_long_t2$Groups == " <LOD" ~ -0.2, 
#                                                        bdd_long_t2$Groups == " LOD-LOQ" ~ -0.2, 
#                                                        
#                                                        bdd_long_t2$Groups == "2nd tertile" ~ 0.2, 
#                                                        bdd_long_t2$Groups == "3rd tertile" ~ 0.2, 
#                                                        bdd_long_t2$Groups == " 3rd tertile" ~ 0.2, 
#                                                        bdd_long_t2$Groups == "LOD_LOQ" ~ 0.2, 
#                                                        bdd_long_t2$Groups == ">LOQ" ~ 0.2, 
#                                                        bdd_long_t2$Groups == " >LOQ" ~ 0.2, 
#                                                        
#                                                        bdd_long_t2$Groups == "2nd tertile - 3rd tertile" ~ 0, 
#                                                        bdd_long_t2$Groups == "LOD-LOQ - >LOQ" ~ 0, 
#                                                        bdd_long_t2$Groups == "<LOD - >LOD" ~ 0, 
#                                                        bdd_long_t2$Groups == "<LOD - >LOQ" ~ 0, 
#                                                        bdd_long_t2$Groups == "1st tertile - 3rd tertile" ~ 0, 
#                                                        bdd_long_t2$Groups == "<LOD - LOD-LOQ" ~ 0, 
#                                                        bdd_long_t2$Groups == "1st tertile - 2nd tertile" ~ 0))) +  
#   scale_fill_manual(
#     values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
#                `<LOD` = "#FFE0DE",
#                ` 1st tertile` = "#FFE0DE", 
#                ` <LOD` = "#FFE0DE",
#                
#                `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
#                `1st tertile - 3rd tertile` = "gray70",
#                `2nd tertile - 3rd tertile` = "gray70",
#                `<LOD - >LOD` = "gray70",
#                `<LOD - >LOQ` = "gray70",
#                `<LOD - LOD-LOQ` = "gray70",
#                `LOD-LOQ - >LOQ`= "gray70",
#                
#                
#                `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
#                `LOD-LOQ` = "#FF8F87",
#                ` 2nd tertile` = "#FF8F87",   
#                ` LOD-LOQ` = "#FF8F87",
#                
#                `3rd tertile` = "#FF4034",   
#                `>LOD` = "#FF4034", 
#                `>LOQ` = "#FF4034", 
#                ` 3rd tertile` = "#FF4034", 
#                ` >LOQ` = "#FF4034")
#   ) +
#   labs(x = "Bray Curtis dissimilarity", 
#        y = "") +
#   theme_lucid() +
#   theme(legend.position = "none", 
#         axis.title = element_text(face = "bold")) +
#   facet_wrap(vars(Pollutant), scales = "free", ncol = 2)
# 
# 
# 
# 
# 
# #### boxplot version3 
# bdd_long_t2 %>%
#   ggplot() +
#   aes(
#     x = Bray_curtis_dissimilarity,
#     y = Groups,
#     fill = Groups
#   ) +
#   geom_boxplot()  +
#   facet_grid(Groups_rec ~ Pollutant, 
#              scales = "free")+  
#   scale_fill_manual(
#     values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
#                `<LOD` = "#FFE0DE",
#                ` 1st tertile` = "#FFE0DE", 
#                ` <LOD` = "#FFE0DE",
#                
#                `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
#                `1st tertile - 3rd tertile` = "gray70",
#                `2nd tertile - 3rd tertile` = "gray70",
#                `<LOD - >LOD` = "gray70",
#                `<LOD - >LOQ` = "gray70",
#                `<LOD - LOD-LOQ` = "gray70",
#                `LOD-LOQ - >LOQ`= "gray70",
#                
#                
#                `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
#                `LOD-LOQ` = "#FF8F87",
#                ` 2nd tertile` = "#FF8F87",   
#                ` LOD-LOQ` = "#FF8F87",
#                
#                `3rd tertile` = "#FF4034",   
#                `>LOD` = "#FF4034", 
#                `>LOQ` = "#FF4034", 
#                ` 3rd tertile` = "#FF4034", 
#                ` >LOQ` = "#FF4034")
#   ) +
#   labs(x = "Bray Curtis dissimilarity", 
#        y = "") +
#   theme_lucid() +
#   theme(legend.position = "none", 
#         axis.title = element_text(face = "bold"))



### Trim3. Phenols ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t3 aux données de dissimilarité
bdd_long_t3 <- bdd_final_t3 %>% select(ident, all_of(explanatory_vars_t3))
bdd_long_t3 <- full_join(bdd_long_t3, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_t3[,explanatory_vars_t3] <- lapply(bdd_long_t3[,explanatory_vars_t3], as.character)
bdd_long_t3 <- bdd_long_t3 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_t3)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_t3_bis <- bdd_final_t3 %>%
  select(ident, 
         all_of(explanatory_vars_t3)) %>%
  rename(ident_bis = ident, 
         mo_MEPA_total_i_cor_t3_ter_bis = mo_MEPA_total_i_cor_t3_ter, 
         mo_ETPA_total_i_cor_t3_ter_bis = mo_ETPA_total_i_cor_t3_ter, 
         mo_PRPA_total_i_cor_t3_ter_bis = mo_PRPA_total_i_cor_t3_ter, 
         mo_BUPA_total_cat_t3_bis = mo_BUPA_total_cat_t3, 
         mo_TRCS_total_i_cor_t3_ter_bis = mo_TRCS_total_i_cor_t3_ter, 
         mo_BPA_total_i_cor_t3_ter_bis = mo_BPA_total_i_cor_t3_ter, 
         mo_BPS_total_cat_t3_2_bis = mo_BPS_total_cat_t3_2, 
         mo_OXBE_total_i_cor_t3_ter_bis = mo_OXBE_total_i_cor_t3_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_t3 <- full_join(bdd_long_t3, bdd_long_t3_bis, by = "ident_bis") 
bdd_long_t3 <- bdd_long_t3 %>%
  filter(!ident %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  filter(!ident_bis %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  filter(!is.na(ident))
rm(bdd_long_t3_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_t3 <- bdd_long_t3 %>%
  mutate(
    MEPA = case_when(mo_MEPA_total_i_cor_t3_ter == "1st tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MEPA_total_i_cor_t3_ter == "2nd tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MEPA_total_i_cor_t3_ter == "3rd tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MEPA_total_i_cor_t3_ter == "1st tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MEPA_total_i_cor_t3_ter == "2nd tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MEPA_total_i_cor_t3_ter == "1st tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MEPA_total_i_cor_t3_ter == "3rd tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MEPA_total_i_cor_t3_ter == "2nd tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MEPA_total_i_cor_t3_ter == "3rd tertile" & mo_MEPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ETPA = case_when(mo_ETPA_total_i_cor_t3_ter == "1st tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_ETPA_total_i_cor_t3_ter == "2nd tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_ETPA_total_i_cor_t3_ter == "3rd tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_ETPA_total_i_cor_t3_ter == "1st tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_ETPA_total_i_cor_t3_ter == "2nd tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_ETPA_total_i_cor_t3_ter == "1st tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_ETPA_total_i_cor_t3_ter == "3rd tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_ETPA_total_i_cor_t3_ter == "2nd tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_ETPA_total_i_cor_t3_ter == "3rd tertile" & mo_ETPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    PRPA = case_when(mo_PRPA_total_i_cor_t3_ter == "1st tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_PRPA_total_i_cor_t3_ter == "2nd tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_PRPA_total_i_cor_t3_ter == "3rd tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_PRPA_total_i_cor_t3_ter == "1st tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_PRPA_total_i_cor_t3_ter == "2nd tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_PRPA_total_i_cor_t3_ter == "1st tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_PRPA_total_i_cor_t3_ter == "3rd tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_PRPA_total_i_cor_t3_ter == "2nd tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_PRPA_total_i_cor_t3_ter == "3rd tertile" & mo_PRPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BUPA = case_when(mo_BUPA_total_cat_t3 == "<LOD" & mo_BUPA_total_cat_t3_bis == "<LOD" ~ "<LOD", 
                     mo_BUPA_total_cat_t3 == "LOD-LOQ" & mo_BUPA_total_cat_t3_bis == "LOD-LOQ" ~ "LOD-LOQ", 
                     mo_BUPA_total_cat_t3 == ">LOQ" & mo_BUPA_total_cat_t3_bis == ">LOQ" ~ ">LOQ", 
                     
                     mo_BUPA_total_cat_t3 == "<LOD" & mo_BUPA_total_cat_t3_bis == "LOD-LOQ" ~ "<LOD - LOD-LOQ", 
                     mo_BUPA_total_cat_t3 == "LOD-LOQ" & mo_BUPA_total_cat_t3_bis == "<LOD" ~ "<LOD - LOD-LOQ", 
                     
                     mo_BUPA_total_cat_t3 == "<LOD" & mo_BUPA_total_cat_t3_bis == ">LOQ" ~ "<LOD - >LOQ", 
                     mo_BUPA_total_cat_t3 == ">LOQ" & mo_BUPA_total_cat_t3_bis == "<LOD" ~ "<LOD - >LOQ", 
                     
                     mo_BUPA_total_cat_t3 == "LOD-LOQ" & mo_BUPA_total_cat_t3_bis == ">LOQ" ~ "LOD-LOQ - >LOQ", 
                     mo_BUPA_total_cat_t3 == ">LOQ" & mo_BUPA_total_cat_t3_bis == "LOD-LOQ" ~ "LOD-LOQ - >LOQ"), 
    
    TRCS = case_when(mo_TRCS_total_i_cor_t3_ter == "1st tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_TRCS_total_i_cor_t3_ter == "2nd tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_TRCS_total_i_cor_t3_ter == "3rd tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_TRCS_total_i_cor_t3_ter == "1st tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_TRCS_total_i_cor_t3_ter == "2nd tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_TRCS_total_i_cor_t3_ter == "1st tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_TRCS_total_i_cor_t3_ter == "3rd tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_TRCS_total_i_cor_t3_ter == "2nd tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_TRCS_total_i_cor_t3_ter == "3rd tertile" & mo_TRCS_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPA = case_when(mo_BPA_total_i_cor_t3_ter == "1st tertile" & mo_BPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                    mo_BPA_total_i_cor_t3_ter == "2nd tertile" & mo_BPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    mo_BPA_total_i_cor_t3_ter == "3rd tertile" & mo_BPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    mo_BPA_total_i_cor_t3_ter == "1st tertile" & mo_BPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    mo_BPA_total_i_cor_t3_ter == "2nd tertile" & mo_BPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    mo_BPA_total_i_cor_t3_ter == "1st tertile" & mo_BPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    mo_BPA_total_i_cor_t3_ter == "3rd tertile" & mo_BPA_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    mo_BPA_total_i_cor_t3_ter == "2nd tertile" & mo_BPA_total_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    mo_BPA_total_i_cor_t3_ter == "3rd tertile" & mo_BPA_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPS = case_when(mo_BPS_total_cat_t3_2 == "<LOD" & mo_BPS_total_cat_t3_2_bis == "<LOD" ~ "<LOD", 
                    mo_BPS_total_cat_t3_2 == ">LOD" & mo_BPS_total_cat_t3_2_bis == ">LOD" ~ ">LOD",  
                    
                    mo_BPS_total_cat_t3_2 == "<LOD" & mo_BPS_total_cat_t3_2_bis == ">LOD" ~ "<LOD - >LOD", 
                    mo_BPS_total_cat_t3_2 == ">LOD" & mo_BPS_total_cat_t3_2_bis == "<LOD" ~ "<LOD - >LOD"), 
    
    OXBE = case_when(mo_OXBE_total_i_cor_t3_ter == "1st tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_OXBE_total_i_cor_t3_ter == "2nd tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_OXBE_total_i_cor_t3_ter == "3rd tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_OXBE_total_i_cor_t3_ter == "1st tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_OXBE_total_i_cor_t3_ter == "2nd tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_OXBE_total_i_cor_t3_ter == "1st tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_OXBE_total_i_cor_t3_ter == "3rd tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_OXBE_total_i_cor_t3_ter == "2nd tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_OXBE_total_i_cor_t3_ter == "3rd tertile" & mo_OXBE_total_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         MEPA, ETPA, PRPA, BUPA, TRCS, BPA, BPS, OXBE)



# on fait passer les données en long 
bdd_long_t3 <- bdd_long_t3 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile, <LOD, LOD-LOQ, >LOQ
test <- 
  bdd_long_t3 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile", "<LOD", "LOD-LOQ", ">LOQ")) %>%  
  filter(!Pollutant == "BPS") %>%
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile", 
                             " <LOD" = "<LOD", 
                             " LOD-LOQ" = "LOD-LOQ", 
                             " >LOQ" = ">LOQ"))
bdd_long_t3 <- bind_rows(bdd_long_t3, test)
rm(test)

bdd_long_t3 <- bdd_long_t3 %>%
  mutate(
    Pollutant = fct_recode(Pollutant,
                           "Methylparaben trim.3" = "MEPA",
                           "Ethylparaben trim.3" = "ETPA",
                           "Propylparaben trim.3" = "PRPA",
                           "Butylparaben trim.3" = "BUPA",
                           "Triclosan trim.3" = "TRCS",
                           "Bisphenol A trim.3" = "BPA",
                           "Bisphenol S trim.3" = "BPS",
                           "Benzophenone 3 trim.3" = "OXBE"),
    Pollutant = fct_relevel(Pollutant, 
                            "Methylparaben trim.3", "Ethylparaben trim.3", "Propylparaben trim.3", "Butylparaben trim.3", 
                            "Triclosan trim.3", "Bisphenol A trim.3", "Bisphenol S trim.3", "Benzophenone 3 trim.3"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile", " >LOQ", ">LOD",                             # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile", "LOD-LOQ - >LOQ",             # inter groupe moyenne et forte exposition
                         " 2nd tertile"," LOD-LOQ",                                   # intra groupe moyenne exposition 
                         
                         "3rd tertile", ">LOQ",                                    # intra groupe forte exposition 
                         "<LOD - >LOD", "<LOD - >LOQ", "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile", " <LOD",
                         
                         "2nd tertile","LOD-LOQ",                                   # intra groupe moyenne exposition 
                         "<LOD - LOD-LOQ", "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile", "<LOD"),                               # intra groupe faible exposition
    Groups_rec = fct_recode(Groups, 
                            "Medium vs High" = " 3rd tertile",
                            "Medium vs High" = " >LOQ",
                            "Low vs Medium" = ">LOD",
                            "Medium vs High" = "2nd tertile - 3rd tertile",
                            "Medium vs High" = "LOD-LOQ - >LOQ",
                            "Medium vs High" = " 2nd tertile",
                            "Medium vs High" = " LOD-LOQ",
                            "Low vs High" = "3rd tertile",
                            "Low vs High" = ">LOQ",
                            "Low vs Medium" = "<LOD - >LOD",
                            "Low vs High" = "<LOD - >LOQ",
                            "Low vs High" = "1st tertile - 3rd tertile",
                            "Low vs High" = " 1st tertile",
                            "Low vs Medium" = " <LOD",
                            "Low vs Medium" = "2nd tertile",
                            "Low vs Medium" = "LOD-LOQ",
                            "Low vs Medium" = "<LOD - LOD-LOQ",
                            "Low vs Medium" = "1st tertile - 2nd tertile",
                            "Low vs Medium" = "1st tertile",
                            "Low vs Medium" = "<LOD"))



#### boxplot version 1 ----
boxplot_t3_phenols <- bdd_long_t3 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `<LOD` = "#FFE0DE",
               ` 1st tertile` = "#FFE0DE", 
               ` <LOD` = "#FFE0DE",
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `<LOD - >LOD` = "gray70",
               `<LOD - >LOQ` = "gray70",
               `<LOD - LOD-LOQ` = "gray70",
               `LOD-LOQ - >LOQ`= "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `LOD-LOQ` = "#FF8F87",
               ` 2nd tertile` = "#FF8F87",   
               ` LOD-LOQ` = "#FF8F87",
               
               `3rd tertile` = "#FF4034",   
               `>LOD` = "#FF4034", 
               `>LOQ` = "#FF4034", 
               ` 3rd tertile` = "#FF4034", 
               ` >LOQ` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(vars(Pollutant), scales = "free", ncol = 2)




### 2 months Phenols ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t3 aux données de dissimilarité
bdd_long_M2 <- bdd_final_M2 %>% select(ident, all_of(explanatory_vars_M2))
bdd_long_M2 <- full_join(bdd_long_M2, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_M2[,explanatory_vars_M2] <- lapply(bdd_long_M2[,explanatory_vars_M2], as.character)
bdd_long_M2 <- bdd_long_M2 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_M2)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_M2_bis <- bdd_final_M2 %>%
  select(ident, 
         all_of(explanatory_vars_M2)) %>%
  rename(ident_bis = ident, 
         ch_MEPA_total_i_cor_M2_ter_bis = ch_MEPA_total_i_cor_M2_ter, 
         ch_ETPA_total_cat_M2_2_bis = ch_ETPA_total_cat_M2_2, 
         ch_PRPA_total_cat_M2_2_bis = ch_PRPA_total_cat_M2_2, 
         ch_BUPA_total_cat_M2_2_bis = ch_BUPA_total_cat_M2_2, 
         ch_TRCS_total_i_cor_M2_ter_bis = ch_TRCS_total_i_cor_M2_ter, 
         ch_BPA_total_i_cor_M2_ter_bis = ch_BPA_total_i_cor_M2_ter, 
         ch_BPS_total_cat_M2_2_bis = ch_BPS_total_cat_M2_2, 
         ch_OXBE_total_i_cor_M2_ter_bis = ch_OXBE_total_i_cor_M2_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_M2 <- full_join(bdd_long_M2, bdd_long_M2_bis, by = "ident_bis") 
bdd_long_M2 <- bdd_long_M2 %>%
  filter(!ident %in% c(14461, 27630, 25673, 25613, 22869, 17717, 15741, 15463, 23330, 22620, 18111, 22126)) %>%
  filter(!ident_bis %in% c(14461, 27630, 25673, 25613, 22869, 17717, 15741, 15463, 23330, 22620, 18111, 22126)) %>%
  filter(!is.na(ident))
rm(bdd_long_M2_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_M2 <- bdd_long_M2 %>%
  mutate(
    MEPA = case_when(ch_MEPA_total_i_cor_M2_ter == "1st tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MEPA_total_i_cor_M2_ter == "2nd tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MEPA_total_i_cor_M2_ter == "3rd tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MEPA_total_i_cor_M2_ter == "1st tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MEPA_total_i_cor_M2_ter == "2nd tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MEPA_total_i_cor_M2_ter == "1st tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MEPA_total_i_cor_M2_ter == "3rd tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MEPA_total_i_cor_M2_ter == "2nd tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MEPA_total_i_cor_M2_ter == "3rd tertile" & ch_MEPA_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ETPA = case_when(ch_ETPA_total_cat_M2_2 == "<LOQ" & ch_ETPA_total_cat_M2_2_bis == "<LOQ" ~ "<LOQ", 
                     ch_ETPA_total_cat_M2_2 == ">LOQ" & ch_ETPA_total_cat_M2_2_bis == ">LOQ" ~ ">LOQ",  
                     
                     ch_ETPA_total_cat_M2_2 == "<LOQ" & ch_ETPA_total_cat_M2_2_bis == ">LOQ" ~ "<LOQ - >LOQ", 
                     ch_ETPA_total_cat_M2_2 == ">LOQ" & ch_ETPA_total_cat_M2_2_bis == "<LOQ" ~ "<LOQ - >LOQ"), 
    
    PRPA = case_when(ch_PRPA_total_cat_M2_2 == "<LOD" & ch_PRPA_total_cat_M2_2_bis == "<LOD" ~ "<LOD", 
                     ch_PRPA_total_cat_M2_2 == ">LOD" & ch_PRPA_total_cat_M2_2_bis == ">LOD" ~ ">LOD",  
                     
                     ch_PRPA_total_cat_M2_2 == "<LOD" & ch_PRPA_total_cat_M2_2_bis == ">LOD" ~ "<LOD - >LOD", 
                     ch_PRPA_total_cat_M2_2 == ">LOD" & ch_PRPA_total_cat_M2_2_bis == "<LOD" ~ "<LOD - >LOD"), 
    
    BUPA = case_when(ch_BUPA_total_cat_M2_2 == "<LOD" & ch_BUPA_total_cat_M2_2_bis == "<LOD" ~ "<LOD", 
                     ch_BUPA_total_cat_M2_2 == ">LOD" & ch_BUPA_total_cat_M2_2_bis == ">LOD" ~ ">LOD",  
                     
                     ch_BUPA_total_cat_M2_2 == "<LOD" & ch_BUPA_total_cat_M2_2_bis == ">LOD" ~ "<LOD - >LOD", 
                     ch_BUPA_total_cat_M2_2 == ">LOD" & ch_BUPA_total_cat_M2_2_bis == "<LOD" ~ "<LOD - >LOD"), 
    
    TRCS = case_when(ch_TRCS_total_i_cor_M2_ter == "1st tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_TRCS_total_i_cor_M2_ter == "2nd tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_TRCS_total_i_cor_M2_ter == "3rd tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_TRCS_total_i_cor_M2_ter == "1st tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_TRCS_total_i_cor_M2_ter == "2nd tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_TRCS_total_i_cor_M2_ter == "1st tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_TRCS_total_i_cor_M2_ter == "3rd tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_TRCS_total_i_cor_M2_ter == "2nd tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_TRCS_total_i_cor_M2_ter == "3rd tertile" & ch_TRCS_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPA = case_when(ch_BPA_total_i_cor_M2_ter == "1st tertile" & ch_BPA_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile", 
                    ch_BPA_total_i_cor_M2_ter == "2nd tertile" & ch_BPA_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    ch_BPA_total_i_cor_M2_ter == "3rd tertile" & ch_BPA_total_i_cor_M2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    ch_BPA_total_i_cor_M2_ter == "1st tertile" & ch_BPA_total_i_cor_M2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    ch_BPA_total_i_cor_M2_ter == "2nd tertile" & ch_BPA_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    ch_BPA_total_i_cor_M2_ter == "1st tertile" & ch_BPA_total_i_cor_M2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    ch_BPA_total_i_cor_M2_ter == "3rd tertile" & ch_BPA_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    ch_BPA_total_i_cor_M2_ter == "2nd tertile" & ch_BPA_total_i_cor_M2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    ch_BPA_total_i_cor_M2_ter == "3rd tertile" & ch_BPA_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPS = case_when(ch_BPS_total_cat_M2_2 == "<LOD" & ch_BPS_total_cat_M2_2_bis == "<LOD" ~ "<LOD", 
                    ch_BPS_total_cat_M2_2 == ">LOD" & ch_BPS_total_cat_M2_2_bis == ">LOD" ~ ">LOD",  
                    
                    ch_BPS_total_cat_M2_2 == "<LOD" & ch_BPS_total_cat_M2_2_bis == ">LOD" ~ "<LOD - >LOD", 
                    ch_BPS_total_cat_M2_2 == ">LOD" & ch_BPS_total_cat_M2_2_bis == "<LOD" ~ "<LOD - >LOD"), 
    
    OXBE = case_when(ch_OXBE_total_i_cor_M2_ter == "1st tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_OXBE_total_i_cor_M2_ter == "2nd tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_OXBE_total_i_cor_M2_ter == "3rd tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_OXBE_total_i_cor_M2_ter == "1st tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_OXBE_total_i_cor_M2_ter == "2nd tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_OXBE_total_i_cor_M2_ter == "1st tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_OXBE_total_i_cor_M2_ter == "3rd tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_OXBE_total_i_cor_M2_ter == "2nd tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_OXBE_total_i_cor_M2_ter == "3rd tertile" & ch_OXBE_total_i_cor_M2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         MEPA, ETPA, PRPA, BUPA, TRCS, BPA, BPS, OXBE)



# on fait passer les données en long 
bdd_long_M2 <- bdd_long_M2 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile, <LOD, LOD-LOQ, >LOQ
test <- 
  bdd_long_M2 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile")) %>%  
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile"))
bdd_long_M2 <- bind_rows(bdd_long_M2, test)
rm(test)

bdd_long_M2 <- bdd_long_M2 %>%
  mutate(
    Pollutant = fct_recode(Pollutant,
                           "Methylparaben 2 months" = "MEPA",
                           "Ethylparaben 2 months" = "ETPA",
                           "Propylparaben 2 months" = "PRPA",
                           "Butylparaben 2 months" = "BUPA",
                           "Triclosan 2 months" = "TRCS",
                           "Bisphenol A 2 months" = "BPA",
                           "Bisphenol S 2 months" = "BPS",
                           "Benzophenone 3 2 months" = "OXBE"),
    Pollutant = fct_relevel(Pollutant, 
                            "Methylparaben 2 months", "Ethylparaben 2 months", "Propylparaben 2 months", "Butylparaben 2 months", 
                            "Triclosan 2 months", "Bisphenol A 2 months", "Bisphenol S 2 months", "Benzophenone 3 2 months"),
    Groups = fct_relevel(Groups,
                         ">LOD", ">LOQ",
                         "<LOD - >LOD", "<LOQ - >LOQ", 
                         "<LOD", "<LOQ",
                         
                         " 3rd tertile",                           # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile",           # inter groupe moyenne et forte exposition
                         " 2nd tertile",                                   # intra groupe moyenne exposition 
                         
                         "3rd tertile",                                    # intra groupe forte exposition 
                         "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile", 
                         
                         "2nd tertile",                                  # intra groupe moyenne exposition 
                         "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile"))                               # intra groupe faible exposition



#### boxplot version 1 ----
boxplot_M2_phenols <- bdd_long_M2 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `<LOD` = "#FFE0DE",
               ` 1st tertile` = "#FFE0DE", 
               `<LOQ` = "#FFE0DE",
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `<LOD - >LOD` = "gray70",
               `<LOQ - >LOQ` = "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               ` 2nd tertile` = "#FF8F87",   
               
               `3rd tertile` = "#FF4034",   
               `>LOD` = "#FF4034", 
               `>LOQ` = "#FF4034", 
               ` 3rd tertile` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(vars(Pollutant), scales = "free", ncol = 2)




### 1 year Phenols ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t3 aux données de dissimilarité
bdd_long_Y1 <- bdd_final_Y1 %>% select(ident, all_of(explanatory_vars_Y1))
bdd_long_Y1 <- full_join(bdd_long_Y1, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_Y1[,explanatory_vars_Y1] <- lapply(bdd_long_Y1[,explanatory_vars_Y1], as.character)
bdd_long_Y1 <- bdd_long_Y1 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_Y1)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_Y1_bis <- bdd_final_Y1 %>%
  select(ident, 
         all_of(explanatory_vars_Y1)) %>%
  rename(ident_bis = ident, 
         ch_MEPA_total_i_cor_Y1_ter_bis = ch_MEPA_total_i_cor_Y1_ter, 
         ch_ETPA_total_i_cor_Y1_ter_bis = ch_ETPA_total_i_cor_Y1_ter, 
         ch_PRPA_total_i_cor_Y1_ter_bis = ch_PRPA_total_i_cor_Y1_ter, 
         ch_BUPA_total_cat_Y1_bis = ch_BUPA_total_cat_Y1, 
         ch_TRCS_total_i_cor_Y1_ter_bis = ch_TRCS_total_i_cor_Y1_ter, 
         ch_BPA_total_i_cor_Y1_ter_bis = ch_BPA_total_i_cor_Y1_ter, 
         ch_BPS_total_cat_Y1_bis = ch_BPS_total_cat_Y1, 
         ch_OXBE_total_i_cor_Y1_ter_bis = ch_OXBE_total_i_cor_Y1_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_Y1 <- full_join(bdd_long_Y1, bdd_long_Y1_bis, by = "ident_bis") 
bdd_long_Y1 <- bdd_long_Y1 %>%
  filter(!ident %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  filter(!ident_bis %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  filter(!is.na(ident))
rm(bdd_long_Y1_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_Y1 <- bdd_long_Y1 %>%
  mutate(
    MEPA = case_when(ch_MEPA_total_i_cor_Y1_ter == "1st tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MEPA_total_i_cor_Y1_ter == "2nd tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MEPA_total_i_cor_Y1_ter == "3rd tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MEPA_total_i_cor_Y1_ter == "1st tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MEPA_total_i_cor_Y1_ter == "2nd tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MEPA_total_i_cor_Y1_ter == "1st tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MEPA_total_i_cor_Y1_ter == "3rd tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MEPA_total_i_cor_Y1_ter == "2nd tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MEPA_total_i_cor_Y1_ter == "3rd tertile" & ch_MEPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ETPA = case_when(ch_ETPA_total_i_cor_Y1_ter == "1st tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_ETPA_total_i_cor_Y1_ter == "2nd tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_ETPA_total_i_cor_Y1_ter == "3rd tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_ETPA_total_i_cor_Y1_ter == "1st tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_ETPA_total_i_cor_Y1_ter == "2nd tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_ETPA_total_i_cor_Y1_ter == "1st tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_ETPA_total_i_cor_Y1_ter == "3rd tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_ETPA_total_i_cor_Y1_ter == "2nd tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_ETPA_total_i_cor_Y1_ter == "3rd tertile" & ch_ETPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    PRPA = case_when(ch_PRPA_total_i_cor_Y1_ter == "1st tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_PRPA_total_i_cor_Y1_ter == "2nd tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_PRPA_total_i_cor_Y1_ter == "3rd tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_PRPA_total_i_cor_Y1_ter == "1st tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_PRPA_total_i_cor_Y1_ter == "2nd tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_PRPA_total_i_cor_Y1_ter == "1st tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_PRPA_total_i_cor_Y1_ter == "3rd tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_PRPA_total_i_cor_Y1_ter == "2nd tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_PRPA_total_i_cor_Y1_ter == "3rd tertile" & ch_PRPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BUPA = case_when(ch_BUPA_total_cat_Y1 == "<LOD" & ch_BUPA_total_cat_Y1_bis == "<LOD" ~ "<LOD", 
                     ch_BUPA_total_cat_Y1 == "LOD-LOQ" & ch_BUPA_total_cat_Y1_bis == "LOD-LOQ" ~ "LOD-LOQ", 
                     ch_BUPA_total_cat_Y1 == ">LOQ" & ch_BUPA_total_cat_Y1_bis == ">LOQ" ~ ">LOQ", 
                     
                     ch_BUPA_total_cat_Y1 == "<LOD" & ch_BUPA_total_cat_Y1_bis == "LOD-LOQ" ~ "<LOD - LOD-LOQ", 
                     ch_BUPA_total_cat_Y1 == "LOD-LOQ" & ch_BUPA_total_cat_Y1_bis == "<LOD" ~ "<LOD - LOD-LOQ", 
                     
                     ch_BUPA_total_cat_Y1 == "<LOD" & ch_BUPA_total_cat_Y1_bis == ">LOQ" ~ "<LOD - >LOQ", 
                     ch_BUPA_total_cat_Y1 == ">LOQ" & ch_BUPA_total_cat_Y1_bis == "<LOD" ~ "<LOD - >LOQ", 
                     
                     ch_BUPA_total_cat_Y1 == "LOD-LOQ" & ch_BUPA_total_cat_Y1_bis == ">LOQ" ~ "LOD-LOQ - >LOQ", 
                     ch_BUPA_total_cat_Y1 == ">LOQ" & ch_BUPA_total_cat_Y1_bis == "LOD-LOQ" ~ "LOD-LOQ - >LOQ"), 
    
    TRCS = case_when(ch_TRCS_total_i_cor_Y1_ter == "1st tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_TRCS_total_i_cor_Y1_ter == "2nd tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_TRCS_total_i_cor_Y1_ter == "3rd tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_TRCS_total_i_cor_Y1_ter == "1st tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_TRCS_total_i_cor_Y1_ter == "2nd tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_TRCS_total_i_cor_Y1_ter == "1st tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_TRCS_total_i_cor_Y1_ter == "3rd tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_TRCS_total_i_cor_Y1_ter == "2nd tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_TRCS_total_i_cor_Y1_ter == "3rd tertile" & ch_TRCS_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPA = case_when(ch_BPA_total_i_cor_Y1_ter == "1st tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                    ch_BPA_total_i_cor_Y1_ter == "2nd tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    ch_BPA_total_i_cor_Y1_ter == "3rd tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    ch_BPA_total_i_cor_Y1_ter == "1st tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    ch_BPA_total_i_cor_Y1_ter == "2nd tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    ch_BPA_total_i_cor_Y1_ter == "1st tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    ch_BPA_total_i_cor_Y1_ter == "3rd tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    ch_BPA_total_i_cor_Y1_ter == "2nd tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    ch_BPA_total_i_cor_Y1_ter == "3rd tertile" & ch_BPA_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    BPS = case_when(ch_BPS_total_cat_Y1 == "<LOD" & ch_BPS_total_cat_Y1_bis == "<LOD" ~ "<LOD", 
                    ch_BPS_total_cat_Y1 == "LOD-LOQ" & ch_BPS_total_cat_Y1_bis == "LOD-LOQ" ~ "LOD-LOQ", 
                    ch_BPS_total_cat_Y1 == ">LOQ" & ch_BPS_total_cat_Y1_bis == ">LOQ" ~ ">LOQ", 
                    
                    ch_BPS_total_cat_Y1 == "<LOD" & ch_BPS_total_cat_Y1_bis == "LOD-LOQ" ~ "<LOD - LOD-LOQ", 
                    ch_BPS_total_cat_Y1 == "LOD-LOQ" & ch_BPS_total_cat_Y1_bis == "<LOD" ~ "<LOD - LOD-LOQ", 
                    
                    ch_BPS_total_cat_Y1 == "<LOD" & ch_BPS_total_cat_Y1_bis == ">LOQ" ~ "<LOD - >LOQ", 
                    ch_BPS_total_cat_Y1 == ">LOQ" & ch_BPS_total_cat_Y1_bis == "<LOD" ~ "<LOD - >LOQ", 
                    
                    ch_BPS_total_cat_Y1 == "LOD-LOQ" & ch_BPS_total_cat_Y1_bis == ">LOQ" ~ "LOD-LOQ - >LOQ", 
                    ch_BPS_total_cat_Y1 == ">LOQ" & ch_BPS_total_cat_Y1_bis == "LOD-LOQ" ~ "LOD-LOQ - >LOQ"), 
    
    OXBE = case_when(ch_OXBE_total_i_cor_Y1_ter == "1st tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_OXBE_total_i_cor_Y1_ter == "2nd tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_OXBE_total_i_cor_Y1_ter == "3rd tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_OXBE_total_i_cor_Y1_ter == "1st tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_OXBE_total_i_cor_Y1_ter == "2nd tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_OXBE_total_i_cor_Y1_ter == "1st tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_OXBE_total_i_cor_Y1_ter == "3rd tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_OXBE_total_i_cor_Y1_ter == "2nd tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_OXBE_total_i_cor_Y1_ter == "3rd tertile" & ch_OXBE_total_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         MEPA, ETPA, PRPA, BUPA, TRCS, BPA, BPS, OXBE)



# on fait passer les données en long 
bdd_long_Y1 <- bdd_long_Y1 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile, <LOD, LOD-LOQ, >LOQ
test <- 
  bdd_long_Y1 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile", "<LOD", "LOD-LOQ", ">LOQ")) %>%  
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile", 
                             " <LOD" = "<LOD", 
                             " LOD-LOQ" = "LOD-LOQ", 
                             " >LOQ" = ">LOQ"))
bdd_long_Y1 <- bind_rows(bdd_long_Y1, test)
rm(test)

bdd_long_Y1 <- bdd_long_Y1 %>%
  mutate(
    Pollutant = fct_recode(Pollutant,
                           "Methylparaben 12 months" = "MEPA",
                           "Ethylparaben 12 months" = "ETPA",
                           "Propylparaben 12 months" = "PRPA",
                           "Butylparaben 12 months" = "BUPA",
                           "Triclosan 12 months" = "TRCS",
                           "Bisphenol A 12 months" = "BPA",
                           "Bisphenol S 12 months" = "BPS",
                           "Benzophenone 3 12 months" = "OXBE"),
    Pollutant = fct_relevel(Pollutant, 
                            "Methylparaben 12 months", "Ethylparaben 12 months", "Propylparaben 12 months", "Butylparaben 12 months", 
                            "Triclosan 12 months", "Bisphenol A 12 months", "Bisphenol S 12 months", "Benzophenone 3 12 months"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile", " >LOQ",                             # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile", "LOD-LOQ - >LOQ",             # inter groupe moyenne et forte exposition
                         " 2nd tertile"," LOD-LOQ",                                   # intra groupe moyenne exposition 
                         
                         "3rd tertile", ">LOQ",                                    # intra groupe forte exposition 
                         "<LOD - >LOQ", "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile", " <LOD",
                         
                         "2nd tertile","LOD-LOQ",                                   # intra groupe moyenne exposition 
                         "<LOD - LOD-LOQ", "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile", "<LOD"))                               # intra groupe faible exposition



#### boxplot version 1 ----
boxplot_Y1_phenols <- bdd_long_Y1 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `<LOD` = "#FFE0DE",
               ` 1st tertile` = "#FFE0DE", 
               ` <LOD` = "#FFE0DE",
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `<LOD - >LOQ` = "gray70",
               `<LOD - LOD-LOQ` = "gray70",
               `LOD-LOQ - >LOQ`= "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `LOD-LOQ` = "#FF8F87",
               ` 2nd tertile` = "#FF8F87",   
               ` LOD-LOQ` = "#FF8F87",
               
               `3rd tertile` = "#FF4034",   
               `>LOQ` = "#FF4034", 
               ` 3rd tertile` = "#FF4034", 
               ` >LOQ` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(vars(Pollutant), scales = "free", ncol = 2)

ggsave(plot = boxplot_t2_phenols, 
       filename = "4_output/betadiv/boxplot_t2_phenols.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)

ggsave(plot = boxplot_t3_phenols, 
       filename = "4_output/betadiv/boxplot_t3_phenols.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)

ggsave(plot = boxplot_M2_phenols, 
       filename = "4_output/betadiv/boxplot_M2_phenols.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)

ggsave(plot = boxplot_Y1_phenols, 
       filename = "4_output/betadiv/boxplot_Y1_phenols.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)


## Detailled (fig4 pour article) ----
### Version 1 : benzo M2 + mepa Y1 + prpa Y1 ----
#### 1st vs 2nd tertile ----
bdd_boxplot_article <- bdd_long_M2 %>%
  filter(Pollutant == "Benzophenone 3 2 months") 

bdd_boxplot_article_bis <- bdd_long_Y1 %>%
  filter(Pollutant %in% c("Methylparaben 12 months", "Propylparaben 12 months"))

bdd_boxplot_article <- 
  bind_rows(bdd_boxplot_article, bdd_boxplot_article_bis) %>%
  filter(!Groups %in% c(" 1st tertile", " 2nd tertile", " 3rd tertile")) %>%
  mutate(
    Groups = fct_relevel(Groups,       
                         ">LOD", ">LOQ", "<LOD - >LOD", "<LOQ - >LOQ", "<LOD", "<LOQ",
                         "3rd tertile", "2nd tertile - 3rd tertile", " 3rd tertile", " 2nd tertile",
         "1st tertile - 3rd tertile", " 1st tertile", "2nd tertile",
        "1st tertile - 2nd tertile", "1st tertile", " >LOQ", "LOD-LOQ - >LOQ",
        " LOD-LOQ", "<LOD - >LOQ", " <LOD", "LOD-LOQ", "<LOD - LOD-LOQ"),
    Pollutant = fct_recode(Pollutant,
                           "Benzophenone 3, 2-month exposure" = "Benzophenone 3 2 months",
                           "Methylparaben, 12-month exposure" = "Methylparaben 12 months",
                           "Propylparaben, 12-month exposure" = "Propylparaben 12 months"))

rm(bdd_boxplot_article_bis)

boxplot_1 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 2nd tertile", "2nd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) +
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes) 
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 2nd tertile", lbl)) {
        bquote("1"^st* " tertile - 2"^nd* " tertile")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile")
      } else {
        lbl
      }
    })
  }) 

my_tag <- c("0.16\n(1.26)", "0.21\n(1.19)", "0.14\n(1.29)") 
boxplot_1 <- tag_facet(boxplot_1, 
          x = 1.05, y = 2, 
          vjust = 0.5, hjust = 0.5,
          open = "", close = "",
          fontface = 1,
          size = 3.5,
          #family = "serif",
          tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 

my_tag_bis <- c("p-value\n(F statistic)", "p-value\n(F statistic)", "p-value\n(F statistic)")
boxplot_1 <- tag_facet(boxplot_1, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 


boxplot_1

#### 1st vs 3rd tertile ----
boxplot_2 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 3rd tertile", lbl)) {
        bquote("1"^st* " tertile - 3"^rd* " tertile")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile")
      } else {
        lbl
      }
    })
  })

my_tag <- c("0.19\n(1.19)", "0.03\n(1.62)", "0.002\n(2.01)") 

boxplot_2 <- tag_facet(boxplot_2, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())      # Supprime le texte des facettes) 

my_tag_bis <- c("       \n             ", "       \n             ", "       \n             ")
boxplot_2 <- tag_facet(boxplot_2, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())     # Supprime le texte des facettes) 


#### 2nd vs 3rd tertile ----
boxplot_3 <- bdd_boxplot_article %>%
  filter(Groups %in% c("2nd tertile", "2nd tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("2nd tertile - 3rd tertile", lbl)) {
        bquote("2"^nd* " tertile - 3"^rd* " tertile")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile")
      } else {
        lbl
      }
    })
  })
my_tag <- c("0.007\n(1.80)", "0.06\n(1.45)", "0.02\n(1.67)") 
boxplot_3 <- tag_facet(boxplot_3, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())    # Supprime le texte des facettes) 


my_tag_bis <- c("       \n             ", "       \n             ", "       \n             ")
boxplot_3 <- tag_facet(boxplot_3, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())    # Supprime le texte des facettes) 


boxplot <- boxplot_1 + boxplot_2 + boxplot_3 + plot_layout(nrow = 3)
boxplot
boxplot_2 + boxplot_3 + plot_layout(ncol = 2)

rm(boxplot_1, boxplot_2, boxplot_3, 
   bdd_long_t2, bdd_long_t3, bdd_long_M2, bdd_long_Y1)

ggsave(plot = boxplot, 
       filename = "4_output/betadiv/boxplot_figure4.tiff", 
       device = "tiff", 
       units = "cm", 
       dpi = 400,
       width = 40, 
       height = 12)

### Version 2 : mepa Y1 + prpa Y1 ----
#### 1st vs 2nd tertile ----
bdd_boxplot_article <- bdd_long_Y1 %>%
  filter(Pollutant %in% c("Methylparaben 12 months", "Propylparaben 12 months"))

bdd_boxplot_article <- bdd_boxplot_article %>%
  filter(!Groups %in% c(" 1st tertile", " 2nd tertile", " 3rd tertile")) %>%
  mutate(
    Groups = fct_relevel(Groups,       
                         "3rd tertile", "2nd tertile - 3rd tertile", " 3rd tertile", " 2nd tertile",
                         "1st tertile - 3rd tertile", " 1st tertile", "2nd tertile",
                         "1st tertile - 2nd tertile", "1st tertile"),
    Pollutant = fct_recode(Pollutant,
                           "Methylparaben, 12-month exposure" = "Methylparaben 12 months",
                           "Propylparaben, 12-month exposure" = "Propylparaben 12 months"))


boxplot_1 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 2nd tertile", "2nd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) +
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes) 
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 2nd tertile", lbl)) {
        bquote("1"^st* " tertile - 2"^nd* " tertile")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile")
      } else {
        lbl
      }
    })
  }) 

my_tag <- c("0.21\n(1.19)", "0.14\n(1.29)") 
boxplot_1 <- tag_facet(boxplot_1, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 

my_tag_bis <- c("p-value\n(F statistic)", "p-value\n(F statistic)")
boxplot_1 <- tag_facet(boxplot_1, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 


boxplot_1

#### 1st vs 3rd tertile ----
boxplot_2 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 3rd tertile", lbl)) {
        bquote("1"^st* " tertile - 3"^rd* " tertile")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile")
      } else {
        lbl
      }
    })
  })

my_tag <- c("0.03\n(1.62)", "0.002\n(2.01)") 

boxplot_2 <- tag_facet(boxplot_2, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())      # Supprime le texte des facettes) 

my_tag_bis <- c("       \n             ", "       \n             ")
boxplot_2 <- tag_facet(boxplot_2, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())     # Supprime le texte des facettes) 


#### 2nd vs 3rd tertile ----
boxplot_3 <- bdd_boxplot_article %>%
  filter(Groups %in% c("2nd tertile", "2nd tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("2nd tertile - 3rd tertile", lbl)) {
        bquote("2"^nd* " tertile - 3"^rd* " tertile")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile")
      } else {
        lbl
      }
    })
  })
my_tag <- c("0.06\n(1.45)", "0.02\n(1.67)") 
boxplot_3 <- tag_facet(boxplot_3, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())    # Supprime le texte des facettes) 


my_tag_bis <- c("       \n             ", "       \n             ")
boxplot_3 <- tag_facet(boxplot_3, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())    # Supprime le texte des facettes) 


boxplot <- boxplot_1 + boxplot_2 + boxplot_3 + plot_layout(nrow = 3)
boxplot


ggsave(plot = boxplot, 
       filename = "4_output/betadiv/boxplot_figure4.tiff", 
       device = "tiff", 
       units = "cm", 
       dpi = 400,
       width = 40, 
       height = 12)


### Version 3 : mepa Y1 + prpa Y1 (textes des axes modifiés) ----
#### 1st vs 2nd tertile ----
bdd_boxplot_article <- bdd_long_Y1 %>%
  filter(Pollutant %in% c("Methylparaben 12 months", "Propylparaben 12 months"))

bdd_boxplot_article <- bdd_boxplot_article %>%
  filter(!Groups %in% c(" 1st tertile", " 2nd tertile", " 3rd tertile")) %>%
  mutate(
    Groups = fct_relevel(Groups,       
                         "3rd tertile", "2nd tertile - 3rd tertile", " 3rd tertile", " 2nd tertile",
                         "1st tertile - 3rd tertile", " 1st tertile", "2nd tertile",
                         "1st tertile - 2nd tertile", "1st tertile"),
    Pollutant = fct_recode(Pollutant,
                           "Methylparaben, 12-month exposure" = "Methylparaben 12 months",
                           "Propylparaben, 12-month exposure" = "Propylparaben 12 months"))


boxplot_1 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 2nd tertile", "2nd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) +
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes) 
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 2nd tertile", lbl)) {
        bquote("1"^st* " tertile - 2"^nd* " tertile intergroup diversity")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile intragroup diversity")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile intragroup diversity")
      } else {
        lbl
      }
    })
  }) 

my_tag <- c("0.21\n(1.19)", "0.14\n(1.29)") 
boxplot_1 <- tag_facet(boxplot_1, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 

my_tag_bis <- c("p-value\n(F statistic)", "p-value\n(F statistic)")
boxplot_1 <- tag_facet(boxplot_1, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 


boxplot_1

#### 1st vs 3rd tertile ----
boxplot_2 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 3rd tertile", lbl)) {
        bquote("1"^st* " tertile - 3"^rd* " tertile intergroup diversity")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile intragroup diversity")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile intragroup diversity")
      } else {
        lbl
      }
    })
  })

my_tag <- c("0.03\n(1.62)", "0.002\n(2.01)") 

boxplot_2 <- tag_facet(boxplot_2, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())      # Supprime le texte des facettes) 

my_tag_bis <- c("       \n             ", "       \n             ")
boxplot_2 <- tag_facet(boxplot_2, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())     # Supprime le texte des facettes) 


#### 2nd vs 3rd tertile ----
boxplot_3 <- bdd_boxplot_article %>%
  filter(Groups %in% c("2nd tertile", "2nd tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  labs(x = "Bray Curtis β-diversity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(vars(Pollutant)) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("2nd tertile - 3rd tertile", lbl)) {
        bquote("2"^nd* " tertile - 3"^rd* " tertile intergroup diversity")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile intragroup diversity")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile intragroup diversity")
      } else {
        lbl
      }
    })
  })
my_tag <- c("0.06\n(1.45)", "0.02\n(1.67)") 
boxplot_3 <- tag_facet(boxplot_3, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())    # Supprime le texte des facettes) 


my_tag_bis <- c("       \n             ", "       \n             ")
boxplot_3 <- tag_facet(boxplot_3, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())    # Supprime le texte des facettes) 


boxplot <- boxplot_1 + boxplot_2 + boxplot_3 + plot_layout(nrow = 3)
boxplot


ggsave(plot = boxplot, 
       filename = "4_output/betadiv/boxplot_figure4_24_08_19.tiff", 
       device = "tiff", 
       units = "cm", 
       dpi = 400,
       width = 40, 
       height = 12)

# Test 1 seule fois vérif ----
set.seed(1996)
results <- 
  adonis2(
    all_dist_t2 ~ mo_MEPA_total_i_cor_t2_ter,
    data = bdd_final_t2,
    permutations = 999)

results <- 
  adonis2(
    all_dist_t2 ~ 
      mo_MEPA_total_i_cor_t2_ter +
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
      # mo_interpreg_3cat +
      mo_dipl_3cat_i +
      po_w_kg_3cat +
      po_he_3cat_i +
      ch_w_Y1_3cat_i +
      ch_he_Y1_3cat_i +
      po_gd +
      mo_age +
      mo_bmi_bepr_3cat_i +
      bf_duration_till48w_4cat_i,
    data = bdd_final_t2,
    permutations = 999
  )





# Corrélations entre les indices d'alpha diversité ----
cor_alpha <- 
  bdd_alpha %>% 
  select(contains("5000")) %>% 
  select(-ch_feces_seChao1_5000_ASV_Y1,-ch_feces_seACE_5000_ASV_Y1, 
         -ch_feces_invSimpson_5000_ASV_Y1) %>%
  rename("Specific richness" = "ch_feces_SpecRich_5000_ASV_Y1", 
         "Chao 1" = "ch_feces_Chao1_5000_ASV_Y1",
         "ACE" = "ch_feces_ACE_5000_ASV_Y1",
         "Shannon" = "ch_feces_Shannon_5000_ASV_Y1", 
         "Simpson" = "ch_feces_Simpson_5000_ASV_Y1", 
         "Fisher" = "ch_feces_Fisher_5000_ASV_Y1", 
         "Faith" = "ch_feces_Faith_5000_ASV_Y1")
heatmap_cor(cormat = cor_alpha, decimal = 2)


# Figure S7 ----
heatmap_alpha <- bdd_alpha %>% 
  select(ident, 
         "Specific richness" = ch_feces_SpecRich_5000_ASV_Y1, 
         "Shannon diversity" = ch_feces_Shannon_5000_ASV_Y1, 
         "Faith diversity" = ch_feces_Faith_5000_ASV_Y1)

heatmap_taxa <- bdd_taxa %>%
  select(ident, 
         "Firmicutes" = ch_feces_rel_p1_Y1, 
         "Actinobacteria" = ch_feces_rel_p2_Y1, 
         "Bacteroidetes" = ch_feces_rel_p3_Y1, 
         "Proteobacteria" = ch_feces_rel_p4_Y1)
heatmap_alpha <- left_join(heatmap_alpha, heatmap_taxa, by = "ident") %>%
  select(-ident)

cormat <- round(cor(heatmap_alpha, 
                    use = "pairwise.complete.obs", 
                    method = "spearman"), 3)

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){              # Utiliser la corrélation entre les variables comme mesure de distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


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
