## Late integration Plasticell modelling ##
## 27/02/2025, Chris O'Grady ##

library(tidyverse)
library(Seurat)
library(glmnet)
library(caret)

source(file.path("01_scripts", "late_int_helpers.R"))

#### scRNA sequencing data  ####
load_and_process_pseudobulk()
generate_nk_proportion()
glm_cv_scRNA <- generate_glm_lasso(pseudo_counts_vg, proportion_NK, vessel_exclude = 1)
model_metrics <- evaluate_model(glm_cv_scRNA$glm_cv, glm_cv_scRNA$counts_test, glm_cv_scRNA$NK_proportion_test)
feature_importance <- get_feature_importance(glm_cv_scRNA$glm_cv)

## Diagnostic plots ##
# Plotting dataframes
r2_values <- pseudo_counts_vg_long %>%
  dplyr::filter(gene %in% rownames(feature_importance)) %>%
  group_by(gene) %>%
  summarise(
    r_squared = summary(lm(proportion_NK ~ expression))$r.squared
  )

# Merge R-squared values with the plot data
pseudo_counts_vg_long_r2 <- pseudo_counts_vg_long %>%
  dplyr::filter(gene %in% rownames(feature_importance)) %>%
  left_join(r2_values, by = "gene")

pseudo_counts_vg_long_r2$gene <- factor(pseudo_counts_vg_long_r2$gene, 
                                        levels = rownames(feature_importance))

## 
ggplot(feature_importance %>% dplyr::filter(Coefficient != 0),
       aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "black") +
  coord_flip() +
  standard_theme() +
  scale_fill_manual(values = c("#CC2936", "#08415C")) +
  labs(title = "Non-zero coefficients in LASSO regression model",
       x = "",
       y = "Coefficient")

# Create the plot
ggplot(pseudo_counts_vg_long_r2, aes(y = proportion_NK, x = expression)) +
  geom_point() +
  facet_wrap(~gene, scales = "free") + 
  geom_smooth(method = "lm") + 
  geom_text(aes(x = Inf, y = Inf, label = paste("R² =", round(r_squared, 2))),
            hjust = 1.1, vjust = 1.1, size = 3) +  # Adjust text placement
  ylab("Proportion of NK cells") + 
  xlab("scRNA expression") +
  standard_theme()

#### O-link ####
## Data formatting
processed_olink <- readRDS("00_data/raw/olink_processed.rds")

# Format O-link data
processed_olink_formatted <- processed_olink$data %>%
  select(vessel, Timepoint, NPX, Assay, label) %>%
  dplyr::filter(!is.na(label)) %>%
  mutate(vessel_day = paste0(str_replace(vessel, "V", "V0"), "_", str_replace(Timepoint, "ay ", ""))) %>%
  select(vessel_day, NPX, Assay) %>%
  group_by(vessel_day, Assay) %>%
  summarise(avg_NPX = mean(NPX, na.rm = TRUE))

# Convert to wide format
processed_olink_wide <- processed_olink_formatted %>%
  pivot_wider(names_from = vessel_day, values_from = avg_NPX) %>%
  column_to_rownames("Assay")

# Load NK proportion data
proportion_NK_olink <- readRDS("00_data/raw/unprocessed/scRNA_processing_files/run4_nk_proportion.rds") %>%
  dplyr::filter(vessel_day %in% colnames(processed_olink_wide))

glm_cv_olink <- generate_glm_lasso(processed_olink_wide, proportion_NK_olink, vessel_exclude = 2)
model_metrics_olink <- evaluate_model(glm_cv_olink$glm_cv, glm_cv_olink$counts_test, glm_cv_olink$NK_proportion_test)
feature_importance_olink <- get_feature_importance(glm_cv_olink$glm_cv)

# Plotting dataframes
r2_values_olink <- processed_olink_formatted %>% 
  left_join(., proportion_NK_olink, by = "vessel_day") %>%
  dplyr::filter(Assay %in% rownames(feature_importance_olink)) %>%
  group_by(Assay) %>%
  summarise(
    r_squared = summary(lm(proportion_NK ~ avg_NPX))$r.squared
  )

r2_values_olink$Assay <- factor(r2_values_olink$Assay, 
                                levels = rownames(feature_importance_olink))

# Merge R-squared values with the plot data
olink_long_formatted <- processed_olink_formatted %>%
  ungroup() %>% 
  dplyr::filter(Assay %in% rownames(feature_importance_olink)) %>%
  left_join(r2_values_olink, by = "Assay") %>% 
  left_join(., proportion_NK_olink, by = "vessel_day")

olink_long_formatted$Assay <- factor(olink_long_formatted$Assay, 
                                     levels = rownames(feature_importance_olink))

## 
ggplot(feature_importance_olink,
       aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "black") +
  coord_flip() +
  standard_theme() +
  scale_fill_manual(values = c("#CC2936", "#08415C")) +
  labs(title = "Non-zero coefficients in LASSO regression model",
       x = "",
       y = "Coefficient")

ggplot(olink_long_formatted, aes(y = proportion_NK, x = avg_NPX)) +
  geom_point() +
  facet_wrap(~Assay, scales = "free") + 
  geom_smooth(method = "lm") + 
  geom_text(data = r2_values_olink, 
            aes(x = Inf, y = Inf, label = paste("R² =", round(r_squared, 2))),
            hjust = 1.1, vjust = 1.1, size = 3) +  # Adjust text placement
  ylab("Proportion of NK cells") + 
  xlab("O-link abundance") +
  standard_theme()

#### Metabolite ####
processed_lcms <- readRDS("00_data/raw/processed_lcms.rds")
filtered_lcms_run4 <- processed_lcms[[2]] %>% 
  dplyr::filter(str_detect(label, "r4")) %>% 
  mutate(vessel_day = toupper(paste0("V0", vessel, "_D", day))) %>% 
  dplyr::filter(vessel_day %in% proportion_NK$vessel_day) %>% 
  left_join(., proportion_NK %>% select(vessel_day, proportion_NK))

filtered_lcms_run4_wide <- filtered_lcms_run4 %>%
  group_by(vessel_day, metabolite) %>%
  summarise(avg_scaled_normalised_peak_area = mean(avg_scaled_normalised_peak_area, na.rm = TRUE)) %>%
  pivot_wider(names_from = vessel_day, values_from = avg_scaled_normalised_peak_area) %>% 
  ungroup() %>% 
  column_to_rownames("metabolite")

proportion_NK_metab <- 
  readRDS(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_nk_proportion.rds")) %>% 
  dplyr::filter(vessel_day %in% unique(filtered_lcms_run4$vessel_day))

## Generate models
glm_cv_metab <- generate_glm_lasso(filtered_lcms_run4_wide, proportion_NK_metab, vessel_exclude = 2)
model_metrics_metab <- evaluate_model(glm_cv_metab$glm_cv, glm_cv_metab$counts_test, glm_cv_metab$NK_proportion_test)
feature_importance_metab <- get_feature_importance(glm_cv_metab$glm_cv)

# Metabolite plotting dataframes
r2_values_metab <- filtered_lcms_run4 %>%
  dplyr::filter(metabolite %in% rownames(feature_importance_metab)) %>%
  group_by(metabolite) %>%
  summarise(
    r_squared = summary(lm(proportion_NK ~ avg_scaled_normalised_peak_area))$r.squared
  )

r2_values_metab$metabolite <- factor(r2_values_metab$metabolite, 
                                     levels = rownames(feature_importance_metab))

# Merge R-squared values with the plot data
metab_long_formatted <- filtered_lcms_run4 %>%
  dplyr::filter(metabolite %in% rownames(feature_importance_metab)) %>%
  left_join(r2_values_metab, by = "metabolite")

metab_long_formatted$metabolite <- factor(metab_long_formatted$metabolite, 
                                          levels = rownames(feature_importance_metab))

##
ggplot(feature_importance_metab,
       aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "black") +
  coord_flip() +
  standard_theme() +
  scale_fill_manual(values = c("#CC2936", "#08415C")) +
  labs(title = "Non-zero coefficients in LASSO regression model",
       x = "",
       y = "Coefficient")

ggplot(metab_long_formatted, aes(y = proportion_NK, x = avg_scaled_normalised_peak_area)) +
  geom_point() +
  facet_wrap(~metabolite, scales = "free") + 
  geom_smooth(method = "lm") + 
  geom_text(data = r2_values_metab, 
            aes(x = Inf, y = Inf, label = paste("R² =", round(r_squared, 2))),
            hjust = 1.1, vjust = 1.1, size = 3) + 
  ylab("Proportion of NK cells") + 
  xlab("Metabolite abundance") +
  standard_theme()


#### miRNA ####
miRNA_long <- readRDS(file.path("00_data", "raw", "processed_miRNA.rds"))[[1]]%>% 
  as_tibble() %>% 
  mutate(vessel_day = paste0("V0", vessel, "_D", day)) %>% 
  dplyr::filter(run == "4") %>% 
  select(vessel_day, everything(), -label, -vessel, -day, -run, -sample_name) %>% 
  pivot_longer(cols = -vessel_day, names_to = "miRNA", values_to = "expression") 

processed_miRNA <- readRDS(file.path("00_data", "raw", "processed_miRNA.rds"))[[1]] %>% 
  as_tibble() %>% 
  mutate(vessel_day = paste0("V0", vessel, "_D", day)) %>% 
  dplyr::filter(run == "4") %>% 
  select(vessel_day, everything(), -label, -vessel, -day, -run, -sample_name) %>% 
  column_to_rownames("vessel_day")

proportion_NK_miRNA <- 
  readRDS(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_nk_proportion.rds")) %>% 
  dplyr::filter(vessel_day %in% rownames(processed_miRNA))


## Generate models
glm_cv_miRNA <- generate_glm_lasso(t(processed_miRNA), proportion_NK_miRNA, vessel_exclude = 1)
model_metrics_miRNA <- evaluate_model(glm_cv_miRNA$glm_cv, glm_cv_miRNA$counts_test, glm_cv_miRNA$NK_proportion_test)
feature_importance_miRNA <- get_feature_importance(glm_cv_miRNA$glm_cv)


## Plots
# Calculate R-squared for each gene
r2_values_miRNA <- miRNA_long %>%
  dplyr::filter(miRNA %in% rownames(feature_importance_miRNA)) %>%
  left_join(., proportion_NK_miRNA %>% dplyr::select(vessel_day, proportion_NK)) %>%
  group_by(miRNA) %>%
  summarise(
    r_squared = summary(lm(proportion_NK ~ expression))$r.squared
  )


# Merge R-squared values with the plot data
processed_miRNA_r2 <- miRNA_long %>%
  dplyr::filter(miRNA %in% rownames(feature_importance_miRNA)) %>%
  left_join(r2_values_miRNA, by = "miRNA") %>%
  left_join(., proportion_NK_miRNA %>% dplyr::select(vessel_day, proportion_NK))

processed_miRNA_r2$miRNA <- factor(processed_miRNA_r2$miRNA, 
                                   levels = rownames(feature_importance_miRNA))

r2_values_miRNA$miRNA <- factor(r2_values_miRNA$miRNA, 
                                levels = rownames(feature_importance_miRNA))

# Create the plot
ggplot(feature_importance_miRNA %>% dplyr::filter(Coefficient != 0),
       aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "black") +
  coord_flip() +
  standard_theme() +
  scale_fill_manual(values = c("#CC2936", "#08415C")) +
  labs(title = "Non-zero coefficients in LASSO regression model",
       x = "",
       y = "Coefficient")

ggplot(processed_miRNA_r2, aes(y = proportion_NK, x = expression)) +
  geom_point() +
  facet_wrap(~miRNA, scales = "free") + 
  geom_smooth(method = "lm") + 
  geom_text(data = r2_values_miRNA, 
            aes(x = Inf, y = Inf, label = paste("R² =", round(r_squared, 2))),
            hjust = 1.1, vjust = 1.1, size = 3) +  # Adjust text placement
  ylab("Proportion of NK cells") + 
  xlab("miRNA count") +
  standard_theme()





