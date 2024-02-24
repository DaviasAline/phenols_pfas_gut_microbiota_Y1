## Aline Davias 
## 2021/10/08
## Metadonnées SEPAGES utilisées dans le cadre Projet GUMME 



# 0_source_data_reading ----
## 1) Chargement des packages ----
library(haven)
library(readr)
library(tidyverse)
library(foreign)
library(questionr)

## 2) Chargement des métadonnées ----
metadata <- read_sas(
  "0_source_data/base_aline_211115.sas7bdat",            # base de données SEPAGES
  catalog_file = "0_source_data/formats.sas7bcat")

bdd_crb <- read_sas(
  "0_source_data/date_selle_aline210520.sas7bdat",       # base avec les dates de prélévement selles IAB (recu séparément)
  NULL)         


## 3) Fusion des métadonnées SEPAGES  reçues séparément ----
### data child (n=484) + data bdd_crb (n=360)
metadata <- metadata %>% rename(ch_feces_ID_Y1 = CodeEchantillon_selle_un_an)
bdd_crb <- bdd_crb %>% rename(ch_feces_ID_Y1 = Code_Echantillon)

metadata <-
  merge(metadata,
        bdd_crb,
        by = c("ident", "ch_feces_ID_Y1"),
        all.x = TRUE) %>%
  mutate(
    ch_feces_ID_Y1 =
      str_replace_all(
        ch_feces_ID_Y1,
        c(
          "SELLE000089.1" = "SELLE0089",
          "SELLE0711.1" = "SELLE0711",
          "SELLE1066.1" = "SELLE1066")),
    ch_feces_ID_Y1 = str_sub(ch_feces_ID_Y1, 1, 9), 
    ident = as.integer(ident))

