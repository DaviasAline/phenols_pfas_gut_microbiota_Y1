## Figures pour brochure SEPAGES 
## 10/07/2023
## Aline Davias 


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
library(expss)
library(patchwork)

test <- bdd_taxa %>%
  select(
    ident,
    ch_feces_rel_g1_Y1:ch_feces_rel_g10_Y1,
    po_delmod,
    ch_antibio_Y1_2cat,
    bf_duration_till48w_4cat) %>%
  na.omit() %>%
  mutate(
    bf_duration_till48w_2cat = fct_recode(bf_duration_till48w_4cat, 
                                          "Non" = "Not breastfed",
                                          "Non" = "<24 weeks",
                                          "Non" = "24-47 weeks",
                                          "Oui" = "Still breastfeed at 48 weeks"), 
    po_delmod = fct_recode(po_delmod,
                           "Césarienne" = "C-section",
                           "Naturel" = "Vaginal delivery"),
    ch_antibio_Y1_2cat = fct_recode(ch_antibio_Y1_2cat, 
                                    "Non" = "No",
                                    "Oui" = "Yes"))

num_var <- test %>% select(ch_feces_rel_g1_Y1:ch_feces_rel_g10_Y1) %>% colnames()

df_summary_delmod <- test %>%
  group_by(po_delmod) %>%
  summarise(
    across(starts_with("ch_feces_rel_"), mean)) %>%
  rename(
    Bifidobacterium  = ch_feces_rel_g1_Y1,
    Bacteroides = ch_feces_rel_g2_Y1, 
    Blautia = ch_feces_rel_g3_Y1,
    "Escherichia and Shigella" = ch_feces_rel_g4_Y1,
    Akkermansia = ch_feces_rel_g5_Y1, 
    Faecalibacterium = ch_feces_rel_g6_Y1, 
    Streptococcus = ch_feces_rel_g7_Y1, 
    "Clostridium XlVa" = ch_feces_rel_g8_Y1, 
    "Lachnospiracea incertae sedis" = ch_feces_rel_g9_Y1, 
    "Anaerostipes" = ch_feces_rel_g10_Y1) %>%
  pivot_longer(cols = -"po_delmod",
               names_to = "genera_name", 
               values_to = "genera_mean") %>%
  mutate(genera_name = fct_relevel(genera_name, 
                                   "Anaerostipes", "Lachnospiracea incertae sedis", "Clostridium XlVa",
                                   "Streptococcus", "Faecalibacterium", "Akkermansia", "Escherichia and Shigella",
                                   "Blautia", "Bacteroides", "Bifidobacterium"))


plot_1 <- ggplot(df_summary_delmod, aes(x = po_delmod, y = genera_mean, fill = genera_name)) +
  geom_bar(stat = "identity") +
  labs(x = "Mode d'accouchement", y = "Abondance relative (%)", fill = "Genres bactériens") +
  #scale_fill_manual(values = c("blue", "red")) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  theme_classic()+
  theme(legend.position = "none", 
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10))


df_summary_bf <- test %>%
  group_by(bf_duration_till48w_2cat) %>%
  summarise(
    across(starts_with("ch_feces_rel_"), mean)) %>%
  rename(
    Bifidobacterium  = ch_feces_rel_g1_Y1,
    Bacteroides = ch_feces_rel_g2_Y1, 
    Blautia = ch_feces_rel_g3_Y1,
    "Escherichia and Shigella" = ch_feces_rel_g4_Y1,
    Akkermansia = ch_feces_rel_g5_Y1, 
    Faecalibacterium = ch_feces_rel_g6_Y1, 
    Streptococcus = ch_feces_rel_g7_Y1, 
    "Clostridium XlVa" = ch_feces_rel_g8_Y1, 
    "Lachnospiracea incertae sedis" = ch_feces_rel_g9_Y1, 
    "Anaerostipes" = ch_feces_rel_g10_Y1) %>%
  pivot_longer(cols = -"bf_duration_till48w_2cat",
               names_to = "genera_name", 
               values_to = "genera_mean") %>%
  mutate(genera_name = fct_relevel(genera_name, 
                                   "Anaerostipes", "Lachnospiracea incertae sedis", "Clostridium XlVa",
                                   "Streptococcus", "Faecalibacterium", "Akkermansia", "Escherichia and Shigella",
                                   "Blautia", "Bacteroides", "Bifidobacterium"))

plot_2 <- ggplot(df_summary_bf, aes(x = bf_duration_till48w_2cat, y = genera_mean, fill = genera_name)) +
  geom_bar(stat = "identity") +
  labs(x = "Allaitement", y = "Abondance relative (%)", fill = "Genres bactériens") +
  #scale_fill_manual(values = c("blue", "red")) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  theme_classic()+
  theme(legend.position = "none", 
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10))



df_summary_antibio <- test %>%
  group_by(ch_antibio_Y1_2cat) %>%
  summarise(
    across(starts_with("ch_feces_rel_"), mean)) %>%
  rename(
    Bifidobacterium  = ch_feces_rel_g1_Y1,
    Bacteroides = ch_feces_rel_g2_Y1, 
    Blautia = ch_feces_rel_g3_Y1,
    "Escherichia and Shigella" = ch_feces_rel_g4_Y1,
    Akkermansia = ch_feces_rel_g5_Y1, 
    Faecalibacterium = ch_feces_rel_g6_Y1, 
    Streptococcus = ch_feces_rel_g7_Y1, 
    "Clostridium XlVa" = ch_feces_rel_g8_Y1, 
    "Lachnospiracea incertae sedis" = ch_feces_rel_g9_Y1, 
    "Anaerostipes" = ch_feces_rel_g10_Y1) %>%
  pivot_longer(cols = -"ch_antibio_Y1_2cat",
               names_to = "genera_name", 
               values_to = "genera_mean")%>%
  mutate(genera_name = fct_relevel(genera_name, 
                                   "Anaerostipes", "Lachnospiracea incertae sedis", "Clostridium XlVa",
                                   "Streptococcus", "Faecalibacterium", "Akkermansia", "Escherichia and Shigella",
                                   "Blautia", "Bacteroides", "Bifidobacterium"))

plot_3 <- ggplot(df_summary_antibio, aes(x = ch_antibio_Y1_2cat, y = genera_mean, fill = genera_name)) +
  geom_bar(stat = "identity") +
  labs(x = "Prise d'antibiotiques 0-1an", y = "Abondance relative (%)", fill = "Genres bactériens") +
  #scale_fill_manual(values = c("blue", "red")) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10), 
        legend.title = element_text(face = "bold"))


plot_brochure <- plot_1 + plot_2 + plot_3 + plot_annotation(tag_levels = "A", tag_suffix = ')')
plot_brochure


ggsave(filename = "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations écrites/6. Brochure SEPAGES/Microbiote_figure2_10072023.tiff",
       plot = plot_brochure, 
       device = "tiff",
       units = "cm",
       width = 22, 
       height = 10, 
       dpi = 400,
       limitsize = FALSE)

