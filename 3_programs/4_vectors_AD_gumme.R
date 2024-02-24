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
alpha_vec <- bdd_alpha %>% select(ch_feces_SpecRich_5000_ASV_Y1, ch_feces_Shannon_5000_ASV_Y1, ch_feces_Faith_5000_ASV_Y1) %>% colnames()
taxa_vec <- c("ch_feces_rel_p1_Y1", "ch_feces_rel_p2_Y1", "ch_feces_rel_p3_Y1", "ch_feces_rel_p4_Y1", 
              "ch_feces_rel_g1_Y1", "ch_feces_rel_g2_Y1", "ch_feces_rel_g3_Y1", "ch_feces_rel_g4_Y1")


pollutant_vec_t2 <- metadata %>% select(all_of(phenols_vec_2)) %>% select(contains("t2")) %>% colnames()
pollutant_vec_t2 <- c(pollutant_vec_t2, pfas_vec)
pollutant_vec_t3 <- metadata %>% select(all_of(phenols_vec_2)) %>% select(contains("t3")) %>% colnames()
pollutant_vec_M2 <- metadata %>% select(all_of(phenols_vec_2)) %>% select(contains("M2")) %>% colnames()
pollutant_vec_Y1 <- metadata %>% select(all_of(phenols_vec_2)) %>% select(contains("Y1")) %>% colnames()

