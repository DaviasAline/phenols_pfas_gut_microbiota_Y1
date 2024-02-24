# Figures Révision EMBO 
# Aline Davias 
# 10/0/2023

# Chargement des packages, des fonctions et des données ---- 
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
load("2_final_data/bdd_bkmr_phenols_pfas.RData")
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')  
library(vegan)
library(esquisse)

# Fig.1 ----
## Fig.1a ----
data_Fig.1a <- bdd_taxa %>% 
  select(ident, 
         starts_with("ch_feces_rel_p")) %>% 
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(ident, 
         ch_feces_rel_p1_Y1:ch_feces_rel_p10_Y1)

bar_order1 <- data_Fig.1a %>% select(ident, ch_feces_rel_p1_Y1)
bar_order2 <- data_Fig.1a %>% select(ident, ch_feces_rel_p2_Y1)
phylum.rel.long <- gather(data_Fig.1a, Variable, value, -ident)                 # = pivot_longer()
phylum.rel.long <- left_join(phylum.rel.long, bar_order1, by= "ident")
phylum.rel.long <- left_join(phylum.rel.long, bar_order2, by= "ident")
phylum.rel.long <- phylum.rel.long %>% 
  rename(Firmicute_order = ch_feces_rel_p1_Y1, 
         Actinobacteria_order = ch_feces_rel_p2_Y1)

phylum.rel.long$Variable <- fct_relevel(
  phylum.rel.long$Variable,
  "ch_feces_rel_p10_Y1", "ch_feces_rel_p9_Y1", "ch_feces_rel_p8_Y1",
  "ch_feces_rel_p7_Y1", "ch_feces_rel_p6_Y1",
  "ch_feces_rel_p5_Y1", "ch_feces_rel_p4_Y1", "ch_feces_rel_p3_Y1",
  "ch_feces_rel_p2_Y1", "ch_feces_rel_p1_Y1")

barplot.phylum.rel <- phylum.rel.long %>%
  ggplot() +
  aes(x = reorder(ident, -Firmicute_order), 
      y = value, 
      fill = Variable) + 
  geom_col(position = position_fill(), 
           width=1) +
  labs(
    x = "Feces of SEPAGES children",
    y = "Relative abundance (%)",
    title = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(
    name = "", 
    labels =  c("ch_feces_rel_p1_Y1" = "Firmicutes",
                "ch_feces_rel_p2_Y1" = "Actinobacteria",
                "ch_feces_rel_p3_Y1" = "Bacteroidetes",
                "ch_feces_rel_p4_Y1" = "Proteobacteria",
                "ch_feces_rel_p5_Y1" = "Verrucomicrobia",
                "ch_feces_rel_p10_Y1" = "Tenericutes",
                "ch_feces_rel_p6_Y1" = "Candidatus Saccharibacteria",
                "ch_feces_rel_p7_Y1" = "Cyanobacteria Chloroplast",
                "ch_feces_rel_p9_Y1" = "Fusobacteria",
                "ch_feces_rel_p8_Y1" = "Euryarchaeota"), 
    palette = "Paired", 
    direction = - 1)+
  scale_x_discrete(breaks = NULL)  +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),       # use  windowsFonts() to see how to set fonts
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=9, face = "bold"),
    axis.text.y = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.key.size = unit(0.5, "cm"))


## Fig.1b ----
data_Fig.1b <- bdd_taxa %>% 
  select(
    ident,
    ch_feces_rel_g1_Y1,
    ch_feces_rel_g2_Y1,
    ch_feces_rel_g3_Y1,
    ch_feces_rel_g4_Y1,
    ch_feces_rel_g5_Y1,
    ch_feces_rel_g6_Y1,
    ch_feces_rel_g7_Y1,
    ch_feces_rel_g8_Y1,
    ch_feces_rel_g9_Y1,
    ch_feces_rel_g10_Y1)  %>%
  filter(!is.na(ch_feces_rel_g1_Y1))%>%
  select(ident, 
         ch_feces_rel_g1_Y1:ch_feces_rel_g10_Y1)

bar_order3 <- data_Fig.1b %>% select(ident, ch_feces_rel_g1_Y1)        
bar_order4 <- data_Fig.1b %>% select(ident, ch_feces_rel_g2_Y1)
genus.rel.long <- gather(genus.rel, Variable, value, -ident) 
genus.rel.long <- left_join(genus.rel.long, bar_order3, by= "ident")
genus.rel.long <- left_join(genus.rel.long, bar_order4, by= "ident")
genus.rel.long <- genus.rel.long %>% 
  rename(
    Bifidobacterium_order = ch_feces_rel_g1_Y1, 
    Bacteroides_order = ch_feces_rel_g2_Y1) %>%
  mutate(
    ident = as.factor(ident), 
    Variable = fct_relevel(
      Variable,
      "ch_feces_rel_g1_Y1", "ch_feces_rel_g2_Y1", "ch_feces_rel_g3_Y1",
      "ch_feces_rel_g4_Y1", "ch_feces_rel_g5_Y1", "ch_feces_rel_g6_Y1",
      "ch_feces_rel_g7_Y1", "ch_feces_rel_g8_Y1", 
      "ch_feces_rel_g9_Y1","ch_feces_rel_g10_Y1"))

boxplot.genus.rel <- genus.rel.long %>%
  ggplot() +
  aes(x = "", 
      y = value, 
      fill = Variable) +
  geom_boxplot(shape = "circle") +
  scale_fill_brewer(
    palette = "Set3", 
    direction = 1, 
    labels = c("ch_feces_rel_g10_Y1" = "Anaerostipes", 
               "ch_feces_rel_g8_Y1" = "Clostridium XlVa", 
               "ch_feces_rel_g9_Y1" = "Lachnospiracea incertae sedis", 
               "ch_feces_rel_g7_Y1" = "Streptococcus",
               "ch_feces_rel_g6_Y1" = "Faecalibacterium",  
               "ch_feces_rel_g5_Y1" = "Akkermansia", 
               "ch_feces_rel_g4_Y1" = "Escherichia Shigella",  
               "ch_feces_rel_g3_Y1" = "Blautia", 
               "ch_feces_rel_g2_Y1" = "Bacteroides",   
               "ch_feces_rel_g1_Y1" = "Bifidobacterium")) +
  labs(
    x = "Most abundant genera",
    title = "",
    y = "Relative abundance (%)", 
    fill = ""
  ) +
  scale_y_continuous()+
  theme_classic() +
  theme(
    text = element_text(family = "serif"), 
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.text = element_text(size = 8, face = "italic"), 
    legend.position = "right", 
    legend.box = "vertical", 
    legend.key.size = unit(0.5, "cm"))

Fig1 <- 
  barplot.phylum.rel + 
  boxplot.genus.rel + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations écrites/4. Article_phenols_pfas_microbiote_Y1/3. Dossier soumission EMBO Mol Med/2. Review/Figures/EMM-2022-17207-V2_Fig1.tiff", 
       plot = Fig1, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 200,
       dpi = 300,
       limitsize = FALSE)


data_Fig.1a <- data_Fig.1a %>%
  select(ch_feces_rel_p1_Y1:ch_feces_rel_p10_Y1)

data_Fig.1b <- data_Fig.1b %>%
  select(ch_feces_rel_g1_Y1:ch_feces_rel_g10_Y1)




# Fig.2 ----
Fig.2 <- results_multi %>%    # from 6_review_EMBO_unipollutants.R
  filter(analysis == "confirmatory") %>% 
  select(model_type, 
         outcome, 
         exposure, 
         exposure_window, 
         term_rec, 
         estimate, 
         std.error, 
         statistic, 
         conf.low, 
         conf.high, 
         p.value)

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

ggsave("4_output/forestplot_alpha_confirmatory.tiff", 
       plot = forestplot_alpha_confirmatory, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 150,
       dpi = 300,
       limitsize = FALSE)


# Fig.3 ----
Fig.3 <- bkmr_risks_overall_results_specrich_num_window_5 %>%
  rename(est_specrich = est, 
         sd_specrich = sd)
Fig.3 <- left_join(Fig.3, 
                   bkmr_risks_overall_results_shannon_num_window_5, 
                   by = "quantile") %>%
  rename(est_shannon = est, 
         sd_shannon = sd)
Fig.3 <- left_join(Fig.3, 
                   bkmr_risks_overall_results_faith_num_window_5, 
                   by = "quantile") %>%
  rename(est_faith = est, 
         sd_faith = sd)

# Fig.EV1 ----

## Rarecurves ----
asv_raw <- 
  read_labelled_csv(
    "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_AD_20220504/1_intermediate_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
  select(
    ident, 
    ch_feces_ID_Y1, 
    starts_with("ch_feces_raw_ASV"))

subsample <- asv_raw %>% 
  na.omit() %>%
  sample_n(100) %>%
  select(-ch_feces_ID_Y1)

dev.off()
alpha_diversity_Y1_rarecurve5000 <- 
  rarecurve(subsample, 
            step = 100, 
            sample = 5000, 
            col = "blue", 
            ylab = "Observed ASVs", 
            xlab = "Depth of sequencing") 
text(x = par("usr")[1], y = par("usr")[4], labels = "(a)", pos = 2, font = 2)


## Scatterplot ----
Fig.EV1 <- bdd_alpha %>%
  select(ident, 
         ch_feces_SpecRich_cmin_ASV_Y1, 
         ch_feces_SpecRich_5000_ASV_Y1, 
         ch_feces_SpecRich_10000_ASV_Y1)
Fig.EV1_bis <- metadata %>%
  select(ident, 
         ch_feces_Nreads_ASVbased_Y1,
         ch_feces_NASV_Y1)
Fig.EV1 <- 
  left_join(Fig.EV1,
            Fig.EV1_bis, 
            by = "ident") %>%
  filter(!is.na(ch_feces_SpecRich_cmin_ASV_Y1))

var_label(Fig.EV1$ch_feces_SpecRich_cmin_ASV_Y1) <- NULL
var_label(Fig.EV1$ch_feces_SpecRich_5000_ASV_Y1) <- NULL
var_label(Fig.EV1$ch_feces_SpecRich_10000_ASV_Y1) <- NULL
var_label(Fig.EV1$ch_feces_Nreads_ASVbased_Y1) <- NULL
var_label(Fig.EV1$ch_feces_NASV_Y1) <- NULL

data_plot_long <- Fig.EV1 %>%
  select(ident,                                           # observed number of ASV before rarefaction
         ch_feces_SpecRich_cmin_ASV_Y1 ,                             # observed number of ASV after rarefaction at the minimal threshold (3580)
         ch_feces_SpecRich_5000_ASV_Y1,                              # observed number of ASV after rarefaction at the threshold 5000
         ch_feces_SpecRich_10000_ASV_Y1, 
         ch_feces_NASV_Y1) %>%                         # observed number of ASV after rarefaction at the threshold 10000
  pivot_longer(
    cols = c("ch_feces_SpecRich_cmin_ASV_Y1", "ch_feces_SpecRich_5000_ASV_Y1", "ch_feces_SpecRich_10000_ASV_Y1"), 
    names_to = "rich_after_rar", 
    values_to = "value") %>%
  rename(
    rich_before_rar = ch_feces_NASV_Y1) %>%
  mutate(
    rich_after_rar = fct_recode(rich_after_rar,
                                "Threshold 10 000" = "ch_feces_SpecRich_10000_ASV_Y1",
                                "Threshold 5 000" = "ch_feces_SpecRich_5000_ASV_Y1",
                                "Minimal threshold 3580" = "ch_feces_SpecRich_cmin_ASV_Y1"), 
    rich_after_rar = fct_relevel(rich_after_rar,
                                 "Minimal threshold 3580", 
                                 "Threshold 5 000", 
                                 "Threshold 10 000"))
scatterplot_1 <- data_plot_long %>%
  filter(rich_after_rar == "Minimal threshold 3580") %>%
ggplot() +
  aes(x = value, y = rich_before_rar) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()+
  labs(x = "Observed ASVs after rarefaction", 
       y = "Observed ASVs before rarefaction")

dev.off()
alpha_diversity_Y1_scatterplot <- data_plot_long %>%
  ggplot() +
  aes(x = value, y = rich_before_rar, colour = rich_after_rar) +
  geom_point(shape = "circle",
             size = 1.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_hue(direction = 1) +
  labs(x = "Observed ASVs after rarefaction", 
       y = "Observed ASVs before rarefaction",
       color = "Threshold of rarefation") +
  theme_classic()

alpha_diversity_Y1_scatterplot

dev.print(device = png, 
          file = "4_output/alpha_diversity_ASVbased_Y1_rarecurve5000.png", 
          width = 900, 
          height = 600)

dev.print(device = png, 
          file = "4_output/alpha_diversity_ASVbased_Y1_scatterplot.png", 
          width = 900, 
          height = 600)


Fig.EV1 %>%  select(-ident)


# Fig.EV2 ----
# DAG

# Fig.EV3 ----
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

ggsave("4_output/forestplot_alpha_exploratory.tiff", 
       plot = forestplot_alpha_exploratory, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 180,
       dpi = 300,
       limitsize = FALSE)

Fig.EV3 <- results_multi %>% 
  filter(analysis == "exploratory") %>% 
  select(model_type, 
         outcome, 
         exposure, 
         exposure_window, 
         term_rec, 
         estimate, 
         std.error, 
         statistic, 
         conf.low, 
         conf.high, 
         p.value)


# Fig.EV4 ----
Fig.EV4 <- bkmr_risks_singvar_specrich_num_window_5 %>%
  rename(est_specrich = est, 
         sd_specrich = sd)
Fig.EV4 <- left_join(Fig.EV6,
                     bkmr_risks_singvar_shannon_num_window_5, 
                     by = c("q.fixed", "variable")) %>%
  rename(est_shannon = est, 
         sd_shannon = sd)
Fig.EV4 <- left_join(Fig.EV6,
                     bkmr_risks_singvar_faith_num_window_5, 
                     by = c("q.fixed", "variable")) %>%
  rename(est_faith = est, 
         sd_faith = sd)


# Fig.S1 ----
phenols_cormat <- metadata[, phenols_vec_num]
colnames(phenols_cormat) <- colnames(phenols_cormat) %>%
  str_replace_all(c(
    "mo_" = "",
    "ch_" = "",
    "_total_i_cor_" = " ",
    "MEPA" = "Methylparaben", 
    "ETPA" = "Ethylparaben", 
    "BUPA" = "Butylparaben", 
    "PRPA" = "Propylparaben", 
    "OXBE" = "Benzophenone 3", 
    "TRCS" = "Triclosan", 
    "BPA" = "Bisphenol A", 
    "BPS" = "Bisphenol S"
  ))


pfas_cormat <- metadata[, pfas_vec_num]
colnames(pfas_cormat) <- colnames(pfas_cormat) %>%
  str_replace_all(c(
    "mo_" = "",
    "ch_" = "",
    "_i_cor" = "",
    "_cor" = ""
  ))
colnames(pfas_cormat) <- paste(colnames(pfas_cormat), "t2", " ")

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

Fig.S1 <- round(cor(phenols_pfas_cormat, 
                     use = "pairwise.complete.obs", 
                     method = "pearson"), 3)
Fig.S1 <- reorder_cormat(Fig.S1)                 # réordonner les coef de cor
Fig.S1 <- get_upper_tri(Fig.S1)               # obtenir que le triangle sup
Fig.S1 <- Fig.S1 %>%
  as.data.frame() %>%
  rownames_to_column(var = "exposure")

# Fig.S2 ----
pollutants_vec_num <- c(phenols_vec_num, pfas_vec_num)
Fig.S2 <- metadata %>% 
  select(all_of(covar_vec_num), 
         all_of(pollutants_vec_num))


Fig.S2 <- round(cor(Fig.S2,
                     use = "pairwise.complete.obs",
                     method = "spearman"), 1)

Fig.S2 <- Fig.S2 %>%
  as.data.frame() %>%
  select(all_of(pollutants_vec_num)) %>%
  t() %>%
  as.data.frame() %>%
  select(all_of(covar_vec_num)) %>%
  as.matrix()

colnames(Fig.S2) <- colnames(Fig.S2) %>%
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

rownames(Fig.S2) <- rownames(Fig.S2) %>%
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

Fig.S2 <- Fig.S2 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "exposure")



# Fig.S3 ----
Fig.S3 <- bdd_taxa %>% 
  select(all_of(taxa_vec)) %>%
  rename(p_firmicutes = ch_feces_rel_p1_Y1, 
         p_actinobacteria = ch_feces_rel_p2_Y1, 
         p_bacteroidetes = ch_feces_rel_p3_Y1, 
         p_proteobacteria = ch_feces_rel_p4_Y1, 
         g_bifidobacterium = ch_feces_rel_g1_Y1, 
         g_bacteroides = ch_feces_rel_g2_Y1,
         g_blautia = ch_feces_rel_g3_Y1, 
         g_escherichia_shigella = ch_feces_rel_g4_Y1)

Fig.S3 <- round(cor(Fig.S3, 
                    use = "pairwise.complete.obs", 
                    method = "spearman"), 3)
Fig.S3 <- reorder_cormat(Fig.S3)                 # réordonner les coef de cor
Fig.S3 <- get_upper_tri(Fig.S3)               # obtenir que le triangle sup

Fig.S3 <- Fig.S3 %>%
  as.data.frame() %>%
  rownames_to_column(var = "taxa")

# Export ----
data_fig <- list(
  Fig.1a, 
  Fig.1b, 
  Fig.2, 
  Fig.3,
  Fig.EV1, 
  Fig.EV3, 
  Fig.EV4, 
  Fig.S1, 
  Fig.S2,
  Fig.S3)

write.xlsx(data_fig, 
           "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/pollutants_gut_microbiota_Y1/4_output/data_fig.xlsx")


# Brochure SEPAGES ----
## barplot général phyla ----
barplot.phylum.rel <- phylum.rel.long %>%
  ggplot() +
  aes(x = reorder(ident, -Firmicute_order), 
      y = value, 
      fill = Variable) + 
  geom_col(position = position_fill(), 
           width=1) +
  labs(
    x = "Selles des enfants SEPAGES à un an",
    y = "Abondance relative (%)",
    title = "") +
  scale_y_continuous(labels = function(x) paste0(x*100,"")) +
  scale_fill_brewer(
    name = "Embranchements", 
    labels =  c("ch_feces_rel_p1_Y1" = "Firmicutes",
                "ch_feces_rel_p2_Y1" = "Actinobacteria",
                "ch_feces_rel_p3_Y1" = "Bacteroidetes",
                "ch_feces_rel_p4_Y1" = "Proteobacteria",
                "ch_feces_rel_p5_Y1" = "Verrucomicrobia",
                "ch_feces_rel_p10_Y1" = "Tenericutes",
                "ch_feces_rel_p6_Y1" = "Candidatus Saccharibacteria",
                "ch_feces_rel_p7_Y1" = "Cyanobacteria Chloroplast",
                "ch_feces_rel_p9_Y1" = "Fusobacteria",
                "ch_feces_rel_p8_Y1" = "Euryarchaeota"), 
    palette = "Paired", 
    direction = - 1)+
  scale_x_discrete(breaks = NULL)  +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=7, face = "bold"),
    axis.text.y = element_text(size =7),
    axis.text.x = element_text(size = 1),
    legend.text = element_text(size = 7), 
    legend.title = element_text(size = 7, face = "bold"),
    legend.key.size = unit(0.2, "cm"))

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations écrites/6. Brochure SEPAGES/plot_general_phyla.tiff", 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 200,
       height = 50)

## phyla ----
bdd_brochure_phyla <- bdd_taxa %>% select(ident, 
                            starts_with("ch_feces_rel_p"), 
                            all_of(covar_vec))

list_po_delmod <- list(
  mean_p1 = aggregate(ch_feces_rel_p1_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p2 = aggregate(ch_feces_rel_p2_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p3 = aggregate(ch_feces_rel_p3_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p4 = aggregate(ch_feces_rel_p4_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p5 = aggregate(ch_feces_rel_p5_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p6 = aggregate(ch_feces_rel_p6_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p7 = aggregate(ch_feces_rel_p7_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p8 = aggregate(ch_feces_rel_p8_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p9 = aggregate(ch_feces_rel_p9_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean),
  mean_p10 = aggregate(ch_feces_rel_p10_Y1 ~ po_delmod, data = bdd_brochure_phyla, FUN = mean)) %>%
  reduce(left_join, by = "po_delmod")%>%
  pivot_longer(cols = -po_delmod, names_to = "taxa") %>%
  mutate(
    po_delmod = fct_recode(
      po_delmod,
      "Césarienne" = "C-section",
      "Voie vaginale" = "Vaginal delivery"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_p10_Y1", "ch_feces_rel_p9_Y1", "ch_feces_rel_p8_Y1",
      "ch_feces_rel_p7_Y1", "ch_feces_rel_p6_Y1",
      "ch_feces_rel_p5_Y1", "ch_feces_rel_p4_Y1", "ch_feces_rel_p3_Y1",
      "ch_feces_rel_p2_Y1", "ch_feces_rel_p1_Y1"))


plot_phyla_delmod <- ggplot(list_po_delmod) +
  aes(x = po_delmod, y = value, fill = taxa) +
  geom_col() +
  labs(
    x = "Mode d'accouchement",
    y = "Abondance relative (%)",
    fill = "Embranchements"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1, 
    labels =  c("ch_feces_rel_p1_Y1" = "Firmicutes",
                "ch_feces_rel_p2_Y1" = "Actinobacteria",
                "ch_feces_rel_p3_Y1" = "Bacteroidetes",
                "ch_feces_rel_p4_Y1" = "Proteobacteria",
                "ch_feces_rel_p5_Y1" = "Verrucomicrobia",
                "ch_feces_rel_p10_Y1" = "Tenericutes",
                "ch_feces_rel_p6_Y1" = "Candidatus Saccharibacteria",
                "ch_feces_rel_p7_Y1" = "Cyanobacteria Chloroplast",
                "ch_feces_rel_p9_Y1" = "Fusobacteria",
                "ch_feces_rel_p8_Y1" = "Euryarchaeota"))+ 
  theme(
    legend.position = "none", 
    axis.title = element_text(size=9, face = "bold"),
    axis.text.y = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.key.size = unit(0.5, "cm"))

list_bf_duration <- list(
  mean_p1 = aggregate(ch_feces_rel_p1_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p2 = aggregate(ch_feces_rel_p2_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p3 = aggregate(ch_feces_rel_p3_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p4 = aggregate(ch_feces_rel_p4_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p5 = aggregate(ch_feces_rel_p5_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p6 = aggregate(ch_feces_rel_p6_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p7 = aggregate(ch_feces_rel_p7_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p8 = aggregate(ch_feces_rel_p8_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p9 = aggregate(ch_feces_rel_p9_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p10 = aggregate(ch_feces_rel_p10_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_phyla, FUN = mean)) %>%
  reduce(left_join, by = "bf_duration_till48w_4cat") %>%
  pivot_longer(cols = -bf_duration_till48w_4cat, names_to = "taxa") %>%
  mutate(
    bf_duration_till48w_2cat = fct_recode(
      bf_duration_till48w_4cat,
      "Allaités" = "<24 weeks",
      "Allaités" = "24-47 weeks",
      "Allaités" = "Still breastfeed at 48 weeks", 
      "Non allaités" = "Not breastfed"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_p10_Y1", "ch_feces_rel_p9_Y1", "ch_feces_rel_p8_Y1",
      "ch_feces_rel_p7_Y1", "ch_feces_rel_p6_Y1",
      "ch_feces_rel_p5_Y1", "ch_feces_rel_p4_Y1", "ch_feces_rel_p3_Y1",
      "ch_feces_rel_p2_Y1", "ch_feces_rel_p1_Y1"),
    value_2cat = ifelse(bf_duration_till48w_2cat == "Allaités", value / 3, value))


plot_phyla_breastfeeding <- ggplot(list_bf_duration) +
  aes(x = bf_duration_till48w_2cat, y = value_2cat, fill = taxa) +
  geom_col() +
  labs(
    x = "Allaitement",
    y = "Abondance relative (%)",
    fill = "Embranchements"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1,
    labels =  c("ch_feces_rel_p1_Y1" = "Firmicutes",
                "ch_feces_rel_p2_Y1" = "Actinobacteria",
                "ch_feces_rel_p3_Y1" = "Bacteroidetes",
                "ch_feces_rel_p4_Y1" = "Proteobacteria",
                "ch_feces_rel_p5_Y1" = "Verrucomicrobia",
                "ch_feces_rel_p10_Y1" = "Tenericutes",
                "ch_feces_rel_p6_Y1" = "Candidatus Saccharibacteria",
                "ch_feces_rel_p7_Y1" = "Cyanobacteria Chloroplast",
                "ch_feces_rel_p9_Y1" = "Fusobacteria",
                "ch_feces_rel_p8_Y1" = "Euryarchaeota"))+ 
  theme(
    legend.position = "none", 
    axis.title = element_text(size=9, face = "bold"),
    axis.text.y = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.key.size = unit(0.5, "cm"))


list_antibio <- list(
  mean_p1 = aggregate(ch_feces_rel_p1_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p2 = aggregate(ch_feces_rel_p2_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p3 = aggregate(ch_feces_rel_p3_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p4 = aggregate(ch_feces_rel_p4_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p5 = aggregate(ch_feces_rel_p5_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p6 = aggregate(ch_feces_rel_p6_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p7 = aggregate(ch_feces_rel_p7_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p8 = aggregate(ch_feces_rel_p8_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p9 = aggregate(ch_feces_rel_p9_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean),
  mean_p10 = aggregate(ch_feces_rel_p10_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_phyla, FUN = mean)) %>%
  reduce(left_join, by = "ch_antibio_Y1_2cat") %>%
  pivot_longer(cols = -ch_antibio_Y1_2cat, names_to = "taxa") %>%
  mutate(
    ch_antibio_Y1_2cat = fct_recode(
      ch_antibio_Y1_2cat,
      "Non" = "No",
      "Oui" = "Yes"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_p10_Y1", "ch_feces_rel_p9_Y1", "ch_feces_rel_p8_Y1",
      "ch_feces_rel_p7_Y1", "ch_feces_rel_p6_Y1",
      "ch_feces_rel_p5_Y1", "ch_feces_rel_p4_Y1", "ch_feces_rel_p3_Y1",
      "ch_feces_rel_p2_Y1", "ch_feces_rel_p1_Y1"))

plot_phyla_antibio <- ggplot(list_antibio) +
  aes(x = ch_antibio_Y1_2cat, y = value, fill = taxa) +
  geom_col() +
  labs(
    x = "Prise d'antibiotiques entre la naissane et l'âge de un an",
    y = "Abondance relative (%)",
    fill = "Embranchements"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1,
    labels =  c("ch_feces_rel_p1_Y1" = "Firmicutes",
                "ch_feces_rel_p2_Y1" = "Actinobacteria",
                "ch_feces_rel_p3_Y1" = "Bacteroidetes",
                "ch_feces_rel_p4_Y1" = "Proteobacteria",
                "ch_feces_rel_p5_Y1" = "Verrucomicrobia",
                "ch_feces_rel_p10_Y1" = "Tenericutes",
                "ch_feces_rel_p6_Y1" = "Candidatus Saccharibacteria",
                "ch_feces_rel_p7_Y1" = "Cyanobacteria Chloroplast",
                "ch_feces_rel_p9_Y1" = "Fusobacteria",
                "ch_feces_rel_p8_Y1" = "Euryarchaeota"))+ 
  theme(
    legend.position = "none", 
    axis.title = element_text(size=9, face = "bold"),
    axis.text.y = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.key.size = unit(0.5, "cm"))




list_sex <- list(
  mean_p1 = aggregate(ch_feces_rel_p1_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p2 = aggregate(ch_feces_rel_p2_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p3 = aggregate(ch_feces_rel_p3_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p4 = aggregate(ch_feces_rel_p4_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p5 = aggregate(ch_feces_rel_p5_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p6 = aggregate(ch_feces_rel_p6_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p7 = aggregate(ch_feces_rel_p7_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p8 = aggregate(ch_feces_rel_p8_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p9 = aggregate(ch_feces_rel_p9_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean),
  mean_p10 = aggregate(ch_feces_rel_p10_Y1 ~ ch_sex, data = bdd_brochure_phyla, FUN = mean)) %>%
  reduce(left_join, by = "ch_sex") %>%
  pivot_longer(cols = -ch_sex, names_to = "taxa") %>%
  mutate(
    ch_sex = fct_recode(
      ch_sex,
      "Garçons" = "Male",
      "Filles" = "Female"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_p10_Y1", "ch_feces_rel_p9_Y1", "ch_feces_rel_p8_Y1",
      "ch_feces_rel_p7_Y1", "ch_feces_rel_p6_Y1",
      "ch_feces_rel_p5_Y1", "ch_feces_rel_p4_Y1", "ch_feces_rel_p3_Y1",
      "ch_feces_rel_p2_Y1", "ch_feces_rel_p1_Y1"))

plot_phyla_sex <- ggplot(list_sex) +
  aes(x = ch_sex, y = value, fill = taxa) +
  geom_col() +
  labs(
    x = "Sexe de l'enfant",
    y = "Abondance relative (%)",
    fill = "Embranchements"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1,
    labels =  c("ch_feces_rel_p1_Y1" = "Firmicutes",
                "ch_feces_rel_p2_Y1" = "Actinobacteria",
                "ch_feces_rel_p3_Y1" = "Bacteroidetes",
                "ch_feces_rel_p4_Y1" = "Proteobacteria",
                "ch_feces_rel_p5_Y1" = "Verrucomicrobia",
                "ch_feces_rel_p10_Y1" = "Tenericutes",
                "ch_feces_rel_p6_Y1" = "Candidatus Saccharibacteria",
                "ch_feces_rel_p7_Y1" = "Cyanobacteria Chloroplast",
                "ch_feces_rel_p9_Y1" = "Fusobacteria",
                "ch_feces_rel_p8_Y1" = "Euryarchaeota")) +
  theme(
    axis.title = element_text(size=9, face = "bold"),
    axis.text.y = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.key.size = unit(0.5, "cm"))


plots_phyla_brochure <- plot_phyla_delmod + plot_phyla_breastfeeding + plot_phyla_antibio + plot_phyla_sex +
  plot_layout(nrow = 1)

## genera ----
bdd_brochure_genera <- 
  bdd_taxa %>% 
  select(ident, 
         starts_with("ch_feces_rel_g"),
         all_of(covar_vec)) 


list_po_delmod <- list(
  mean_g1 = aggregate(ch_feces_rel_g1_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g2 = aggregate(ch_feces_rel_g2_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g3 = aggregate(ch_feces_rel_g3_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g4 = aggregate(ch_feces_rel_g4_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g5 = aggregate(ch_feces_rel_g5_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g6 = aggregate(ch_feces_rel_g6_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g7 = aggregate(ch_feces_rel_g7_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g8 = aggregate(ch_feces_rel_g8_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g9 = aggregate(ch_feces_rel_g9_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean),
  mean_g10 = aggregate(ch_feces_rel_g10_Y1 ~ po_delmod, data = bdd_brochure_genera, FUN = mean)) %>%
  reduce(left_join, by = "po_delmod")%>%
  pivot_longer(cols = -po_delmod, names_to = "taxa") %>%
  mutate(
    po_delmod = fct_recode(
      po_delmod,
      "Césarienne" = "C-section",
      "Naturel" = "Vaginal delivery"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_g10_Y1", "ch_feces_rel_g9_Y1", "ch_feces_rel_g8_Y1",
      "ch_feces_rel_g7_Y1", "ch_feces_rel_g6_Y1",
      "ch_feces_rel_g5_Y1", "ch_feces_rel_g4_Y1", "ch_feces_rel_g3_Y1",
      "ch_feces_rel_g2_Y1", "ch_feces_rel_g1_Y1"))


plot_genera_delmod <- ggplot(list_po_delmod) +
  aes(x = po_delmod, y = value, fill = taxa) +
  geom_col() +
  labs(
    x = "Mode d'accouchement",
    y = "Abondance relative (%)",
    fill = "Genres bactériens"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1, 
    labels =  c("ch_feces_rel_g10_Y1" = "Anaerostipes", 
                "ch_feces_rel_g8_Y1" = "Clostridium XlVa", 
                "ch_feces_rel_g9_Y1" = "Lachnospiracea incertae sedis", 
                "ch_feces_rel_g7_Y1" = "Streptococcus",
                "ch_feces_rel_g6_Y1" = "Faecalibacterium",  
                "ch_feces_rel_g5_Y1" = "Akkermansia", 
                "ch_feces_rel_g4_Y1" = "Escherichia Shigella",  
                "ch_feces_rel_g3_Y1" = "Blautia", 
                "ch_feces_rel_g2_Y1" = "Bacteroides",   
                "ch_feces_rel_g1_Y1" = "Bifidobacterium"))+ 
  theme(
    legend.position = "none", 
    axis.title = element_text(size=7, face = "bold"),
    axis.text = element_text(size =7), 
    legend.text = element_text(size = 7), 
    legend.key.size = unit(0.2, "cm"), 
    legend.title = element_text(size = 7, face = "bold"))

list_bf_duration <- list(
  mean_g1 = aggregate(ch_feces_rel_g1_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g2 = aggregate(ch_feces_rel_g2_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g3 = aggregate(ch_feces_rel_g3_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g4 = aggregate(ch_feces_rel_g4_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g5 = aggregate(ch_feces_rel_g5_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g6 = aggregate(ch_feces_rel_g6_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g7 = aggregate(ch_feces_rel_g7_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g8 = aggregate(ch_feces_rel_g8_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g9 = aggregate(ch_feces_rel_g9_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean),
  mean_g10 = aggregate(ch_feces_rel_g10_Y1 ~ bf_duration_till48w_4cat, data = bdd_brochure_genera, FUN = mean)) %>%
  reduce(left_join, by = "bf_duration_till48w_4cat") %>%
  pivot_longer(cols = -bf_duration_till48w_4cat, names_to = "taxa") %>%
  mutate(
    bf_duration_till48w_2cat = fct_recode(
      bf_duration_till48w_4cat,
      "Allaités" = "<24 weeks",
      "Allaités" = "24-47 weeks",
      "Allaités" = "Still breastfeed at 48 weeks", 
      "Non allaités" = "Not breastfed"),
    bf_duration_till48w_2cat_bis = fct_recode(
      bf_duration_till48w_4cat,
      "Non allaités" = "<24 weeks",
      "Non allaités" = "24-47 weeks",
      "Allaités" = "Still breastfeed at 48 weeks", 
      "Non allaités" = "Not breastfed"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_g10_Y1", "ch_feces_rel_g9_Y1", "ch_feces_rel_g8_Y1",
      "ch_feces_rel_g7_Y1", "ch_feces_rel_g6_Y1",
      "ch_feces_rel_g5_Y1", "ch_feces_rel_g4_Y1", "ch_feces_rel_g3_Y1",
      "ch_feces_rel_g2_Y1", "ch_feces_rel_g1_Y1"),
    value_2cat = ifelse(bf_duration_till48w_2cat == "Allaités", value / 3, value), 
    value_2cat_bis = ifelse(bf_duration_till48w_2cat_bis == "Non allaités", value/3, value))



plot_genera_breastfeeding_bis <- ggplot(list_bf_duration) +
  aes(x = bf_duration_till48w_2cat_bis, y = value_2cat_bis, fill = taxa) +
  geom_col() +
  labs(
    x = "Allaitement",
    y = "Abondance relative (%)",
    fill = "Genres bactériens"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1,
    labels =  c("ch_feces_rel_g10_Y1" = "Anaerostipes", 
                "ch_feces_rel_g8_Y1" = "Clostridium XlVa", 
                "ch_feces_rel_g9_Y1" = "Lachnospiracea incertae sedis", 
                "ch_feces_rel_g7_Y1" = "Streptococcus",
                "ch_feces_rel_g6_Y1" = "Faecalibacterium",  
                "ch_feces_rel_g5_Y1" = "Akkermansia", 
                "ch_feces_rel_g4_Y1" = "Escherichia Shigella",  
                "ch_feces_rel_g3_Y1" = "Blautia", 
                "ch_feces_rel_g2_Y1" = "Bacteroides",   
                "ch_feces_rel_g1_Y1" = "Bifidobacterium"))+ 
  theme(
    legend.position = "none", 
    axis.title = element_text(size=7, face = "bold"),
    axis.text = element_text(size =7), 
    legend.text = element_text(size = 7), 
    legend.key.size = unit(0.2, "cm"), 
    legend.title = element_text(size = 7, face = "bold"))


list_antibio <- list(
  mean_g1 = aggregate(ch_feces_rel_g1_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g2 = aggregate(ch_feces_rel_g2_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g3 = aggregate(ch_feces_rel_g3_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g4 = aggregate(ch_feces_rel_g4_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g5 = aggregate(ch_feces_rel_g5_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g6 = aggregate(ch_feces_rel_g6_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g7 = aggregate(ch_feces_rel_g7_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g8 = aggregate(ch_feces_rel_g8_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g9 = aggregate(ch_feces_rel_g9_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean),
  mean_g10 = aggregate(ch_feces_rel_g10_Y1 ~ ch_antibio_Y1_2cat, data = bdd_brochure_genera, FUN = mean)) %>%
  reduce(left_join, by = "ch_antibio_Y1_2cat") %>%
  pivot_longer(cols = -ch_antibio_Y1_2cat, names_to = "taxa") %>%
  mutate(
    ch_antibio_Y1_2cat = fct_recode(
      ch_antibio_Y1_2cat,
      "Non" = "No",
      "Oui" = "Yes"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_g10_Y1", "ch_feces_rel_g9_Y1", "ch_feces_rel_g8_Y1",
      "ch_feces_rel_g7_Y1", "ch_feces_rel_g6_Y1",
      "ch_feces_rel_g5_Y1", "ch_feces_rel_g4_Y1", "ch_feces_rel_g3_Y1",
      "ch_feces_rel_g2_Y1", "ch_feces_rel_g1_Y1"))

plot_genera_antibio <- ggplot(list_antibio) +
  aes(x = ch_antibio_Y1_2cat, y = value, fill = taxa) +
  geom_col() +
  labs(
    x = "Prise d'antibiotiques 0-1an",
    y = "Abondance relative (%)",
    fill = "Genres bactériens"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1,
    labels =  c("ch_feces_rel_g10_Y1" = "Anaerostipes", 
                "ch_feces_rel_g8_Y1" = "Clostridium XlVa", 
                "ch_feces_rel_g9_Y1" = "Lachnospiracea incertae sedis", 
                "ch_feces_rel_g7_Y1" = "Streptococcus",
                "ch_feces_rel_g6_Y1" = "Faecalibacterium",  
                "ch_feces_rel_g5_Y1" = "Akkermansia", 
                "ch_feces_rel_g4_Y1" = "Escherichia Shigella",  
                "ch_feces_rel_g3_Y1" = "Blautia", 
                "ch_feces_rel_g2_Y1" = "Bacteroides",   
                "ch_feces_rel_g1_Y1" = "Bifidobacterium"))+ 
  theme(
    #legend.position = "none", 
    axis.title = element_text(size=7, face = "bold"),
    axis.text = element_text(size =7), 
    legend.text = element_text(size = 7), 
    legend.key.size = unit(0.2, "cm"), 
    legend.title = element_text(size = 7, face = "bold"))




list_sex <- list(
  mean_g1 = aggregate(ch_feces_rel_g1_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g2 = aggregate(ch_feces_rel_g2_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g3 = aggregate(ch_feces_rel_g3_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g4 = aggregate(ch_feces_rel_g4_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g5 = aggregate(ch_feces_rel_g5_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g6 = aggregate(ch_feces_rel_g6_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g7 = aggregate(ch_feces_rel_g7_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g8 = aggregate(ch_feces_rel_g8_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g9 = aggregate(ch_feces_rel_g9_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean),
  mean_g10 = aggregate(ch_feces_rel_g10_Y1 ~ ch_sex, data = bdd_brochure_genera, FUN = mean)) %>%
  reduce(left_join, by = "ch_sex") %>%
  pivot_longer(cols = -ch_sex, names_to = "taxa") %>%
  mutate(
    ch_sex = fct_recode(
      ch_sex,
      "Garçons" = "Male",
      "Filles" = "Female"),
    taxa = fct_relevel(
      taxa,
      "ch_feces_rel_g10_Y1", "ch_feces_rel_g9_Y1", "ch_feces_rel_g8_Y1",
      "ch_feces_rel_g7_Y1", "ch_feces_rel_g6_Y1",
      "ch_feces_rel_g5_Y1", "ch_feces_rel_g4_Y1", "ch_feces_rel_g3_Y1",
      "ch_feces_rel_g2_Y1", "ch_feces_rel_g1_Y1"))

plot_genera_sex <- ggplot(list_sex) +
  aes(x = ch_sex, y = value, fill = taxa) +
  geom_col() +
  labs(
    x = "Sexe de l'enfant",
    y = "Abondance relative (%)",
    fill = "Genres bactériens"
  ) +
  theme_classic()+
  scale_fill_brewer(
    type = "seq",
    palette = "Paired",
    direction = -1,
    labels =  c("ch_feces_rel_g10_Y1" = "Anaerostipes", 
                "ch_feces_rel_g8_Y1" = "Clostridium XlVa", 
                "ch_feces_rel_g9_Y1" = "Lachnospiracea incertae sedis", 
                "ch_feces_rel_g7_Y1" = "Streptococcus",
                "ch_feces_rel_g6_Y1" = "Faecalibacterium",  
                "ch_feces_rel_g5_Y1" = "Akkermansia", 
                "ch_feces_rel_g4_Y1" = "Escherichia Shigella",  
                "ch_feces_rel_g3_Y1" = "Blautia", 
                "ch_feces_rel_g2_Y1" = "Bacteroides",   
                "ch_feces_rel_g1_Y1" = "Bifidobacterium"))+
  theme(
    axis.title = element_text(size=7, face = "bold"),
    axis.text = element_text(size =7), 
    legend.text = element_text(size =7), 
    legend.key.size = unit(0.2, "cm"), 
    legend.title = element_text(size = 7, face = "bold"))


plots_genera_brochure <- plot_genera_delmod + plot_genera_breastfeeding_bis + plot_genera_antibio + 
  plot_layout(nrow = 1)


ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations écrites/6. Brochure SEPAGES/plots_genera_brochure.tiff", 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 200,
       height = 50)
