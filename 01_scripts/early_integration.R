### Using the early integration approach to integrate the single-cell RNAseq data, OLink Proteomics data, miRNA data, and LCMS data.

## init

library(uwot)
conflicts_prefer(uwot::umap)
conflicts_prefer(base::intersect)
conflicts_prefer(base::unname)

## Combining all samples

#### Preliminary data wrangling

k_lcms <- processed_lcms$both_runs_data_wide |>
  select(-run, -day, -vessel, -sample_type) |>
  rename_all(\(x) ifelse(x == "label", x, paste0("LCMS_", x)))
k_mirna <- processed_miRNA$processed_miRNA |>
  select(-vessel, -sample_name, -day, -run) |>
  rename_all(\(x) ifelse(x == "label", x, paste0("miRNA_", x)))
k_olink <- processed_olink$data |>
  filter(DF == "Neat", Assay %in% processed_olink$m_assays) |>
  select(label, Assay, NPX) |>
  drop_na(label, NPX) |>
  pivot_wider(names_from = Assay, values_from = NPX, values_fn = \(x) mean(x, na.rm = TRUE)) |>
  rename_all(\(x) ifelse(x == "label", x, paste0("olink_", x)))
k_rnaseq <- combined_pseudobulk_df |>
  select(label, names(combined_run3_run4_variable_genes_union_df)[-(c(1:47, 3741, 3742))]) |>
  rename_all(\(x) ifelse(x == "label", x, paste0("scRNA_", x))) |>
  group_by(label) |>
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

#### Merging the datasets

early_data <- k_rnaseq |>
  full_join(k_olink, by = "label") |>
  full_join(k_mirna, by = "label") |>
  full_join(k_lcms, by = "label")

#### Scaling

early_data <- early_data |>
  mutate(across(c(everything(), -label), ~ scale(.x) |> as.numeric()))

#### Percentage missing data overall

early_data |>
  summarise(across(everything(), ~ sum(is.na(.x))/n())) |>
  pivot_longer(everything(), names_to = "feature", values_to = "missing") |>
  arrange(desc(missing)) |>
  mutate(
    omics = str_extract(feature, "^(.*?)(?=_)", group = 1)
  ) |>
  group_by(omics) |>
  summarise(
    missing = mean(missing) |> scales::percent(accuracy = 0.1)
  ) |>
  drop_na(omics)

#### Imputation
############## Don't
# library(mice)
# 
# early_data_imputed <- early_data |>
#   mice(m = 1, seed = 123) |>
#   complete()
#   

#### Dimensionality reduction

library(FactoMineR)
library(factoextra)
library(missMDA)
library(aqua)

k_days <- early_data$label |>
  str_extract("d([0-9]+)", group = 1) |>
  factor(
    levels = c(0:6, 8, 10, 12, 13, 15, 19, 22, 24, 26, 29, 31, 33, 35),
    labels = c(
      "0", "mesoderm commitment", "mesoderm commitment", "mesoderm commitment", "haemogenic endothelium", "haemogenic endothelium", "haemogenic endothelium", "early haematopoietic stem cells", "early haematopoietic stem cells", "early haematopoietic stem cells", "late haematopoietic stem cells", "late haematopoietic stem cells", "early Natural Killer cells", "early Natural Killer cells", "early Natural Killer cells", "early Natural Killer cells", "late Natural Killer cells", "late Natural Killer cells", "late Natural Killer cells", "late Natural Killer cells"
    )
  )
 # factor(levels = c(0:6, 8, 10, 12, 13, 15, 19, 22, 24, 26, 29, 31, 33, 35), labels = c("0", "1.5", "1.5", "3.5", "3.5", "5.5", "5.5", "~9", "~9", "12", "~14", "~14", "~21", "~21", "~25", "~25", "29", "~32", "~32", "35"))
  
early_data_PCA <- early_data |>
  select(-label) |>
  missMDA::imputePCA() |>
  magrittr::extract("completeObs") |>
  PCA(graph = FALSE)

fviz_pca_ind(early_data_PCA, habillage = k_days, addEllipses = TRUE, geom.ind = "point") +
  theme_catapult() +
  scale_colour_catapult_d(option = "H") +
  scale_fill_catapult_d(option = "H") +
  guides(
    colour = guide_legend(nrow = 3),
    shape = guide_legend(nrow = 3)
  )

fviz_pca_var(early_data_PCA, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), select.var = list("cos2" = 0.985), repel = TRUE) +
  theme_catapult()

fviz_pca_ind(early_data_PCA, axes = c(2,3), habillage = k_days, addEllipses = TRUE, geom.ind = "point") +
  theme_catapult() +
  scale_colour_catapult_d(option = "H") +
  scale_fill_catapult_d(option = "H")

early_data_umap <- early_data |>
  select(-label) |>
  missMDA::imputePCA() |>
  magrittr::extract2("completeObs")

u1 <- umap(early_data_umap |> as.matrix(), n_threads = 32, n_neighbors = 15, n_epochs = 1000, min_dist = 0.1, spread = 1, verbose = TRUE)

p0 <- bind_cols(k_days, u1, early_data$label) |>
  set_names(c("day", "u1", "u2", "label")) |>
  mutate(vessel = str_extract(label, "v[0-9]")) |>
  ggplot(aes(x = u1, y = u2, color = day, fill = day, text = label, shape = vessel)) +
  geom_point(alpha = 3/4, size = 4) +
  theme_catapult() +
 # stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/4) +
  scale_colour_catapult_d(option = "H") +
  scale_fill_catapult_d(option = "H") +
  labs(
    title = "UMAP of early integration data",
    subtitle = "Coloured by day"
  ) +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1), ncol = 3),
    fill = guide_legend(override.aes = list(alpha = 1), ncol = 3)
  )

plotly::ggplotly(p0)


# find common targets across omics
k_targets_olink_scrnaseq <- early_data |> 
  names() |> 
  as_tibble() |> 
  mutate(
    omic = str_extract(value, "^([a-zA-Z]+)_", group = 1), 
    target = str_extract(value, "_([A-Za-z0-9\\-\\._ \\']+)", group = 1)
  ) |> 
  select(omic, target) |> 
  group_by(omic) |> 
  nest() |> 
  filter(omic %in% c("scRNA", "olink"))

common_targets <- intersect(k_targets_olink_scrnaseq[[1, "data"]] |> unlist() |> unname(), k_targets_olink_scrnaseq[[2, "data"]] |> unlist() |> unname())

target_names <- c(paste0("scRNA_", common_targets), paste0("olink_", common_targets))

# correlate common targets
common_targets_cor <- early_data |>
  select(label, all_of(target_names)) |>
  mutate(days = k_days) |>
  pivot_longer(cols = -c(label, days), names_to = c("omics", "feature"), values_to = "value", names_pattern = "^([a-zA-Z]+)_([A-Za-z0-9\\-\\._ \\']+)$") |>
  pivot_wider(names_from = omics, values_from = value, values_fn = median) |>
  group_by(feature) |>
  summarise(correlation = cor(scRNA, olink, method = "kendall", use = "pairwise.complete.obs")) |>
  arrange(desc(abs(correlation)))


#' @export
plot_common_targets <- function(data, target) {
  
  target_names <- c(paste0("scRNA_", target), paste0("olink_", target))
  
  data |>
    select(label, all_of(target_names)) |>
    mutate(days = k_days) |>
    pivot_longer(cols = c(-label, -days), names_to = c("omics", "feature"), values_to = "value", names_pattern = "^([a-zA-Z]+)_([A-Za-z0-9\\-\\._ \\']+)$") |>
    pivot_wider(names_from = omics, values_from = value) |>
    ggplot(aes(x = scRNA, y = olink, colour = days)) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE, colour = "black") +
    theme_catapult() +
    labs(
      title = target,
      subtitle = "Correlation between scRNA and OLink proteomics data",
      x = "scRNA (standardised)",
      y = "OLink (standardised)"
    ) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 1), nrow = 3)
    )

}

plot_common_targets(early_data, "KRT18")

library(parsnip)
library(Boruta)

# Predictions
pred_data_1 <- processed_process_data$cytotoxicity_process_data_summarised_wide |>
  left_join(early_data, by = c("label" = "label")) |> # drop columns that only have missing data
  select(where(~ sum(!is.na(.x)) == 10))

k_b1 <- Boruta(pred_data_1[, -(1:9)], pred_data_1$combined_cytotoxicity, doTrace = 2)
k_b1_fixed <- TentativeRoughFix(k_b1)
plot(k_b1)


### iterations through output

pred_data_2 <- processed_process_data$all_process_data_unscaled_added_metrics |>
  left_join(early_data, by = c("label" = "label"))

variables_to_predict <- names(pred_data_2)[6:63]


get_important_features <- function(data, feature) {

  if (length(feature) != 1) {
    stop("Only one feature can be selected")
  }
  
  data_mini <- data |> 
    tidyr::drop_na(dplyr::all_of(feature))
  
  k_b2 <- Boruta::Boruta(data_mini[, -(1:63)], data_mini[[feature]], doTrace = 0)

  k_b2
  
}

featureSets <- vector(mode = "list", length = length(variables_to_predict))
names(featureSets) <- variables_to_predict

library(foreach)
library(doParallel)

cl <- makeCluster(32)
registerDoParallel(cl)

my_output <- foreach(i = seq_along(variables_to_predict)) %dopar% {
  cli::cat_line(paste0("Checking feature: ", variables_to_predict[i], "(", i, " out of ", length(variables_to_predict), ")"))
  k_f <- get_important_features(pred_data_2, variables_to_predict[i])
  if (length(k_f$finalDecision[which(k_f$finalDecision == "Confirmed")]) > 0) {
    cli::cat_bullet(paste0("Features confirmed for ", variables_to_predict[i], ": ", paste(names(k_f$finalDecision[which(k_f$finalDecision == "Confirmed")]), collapse = ", ")), bullet = "•", col = "green")
    # featureSets[[i]] <- names(k_f$finalDecision[which(k_f$finalDecision == "Confirmed")])
  } else {
    cli::cat_bullet(paste0("No features confirmed for ", variables_to_predict[i]), bullet = "•", col = "darkorange")
  }
  out <- list(names(k_f$finalDecision[which(k_f$finalDecision == "Confirmed")]))
  names(out) <- variables_to_predict[i]
  out
}

stopCluster(cl)

pred_data_2 |> 
  select(day, vessel, scRNA_CFH, `Nanog+ (%)`) |>
  drop_na(scRNA_CFH) |>
  mutate(
    day = factor(day,
      levels = c(0:6, 8, 10, 12, 13, 15, 19, 22, 24, 26, 29, 31, 33, 35),
      labels = c(
        "0", "mesoderm commitment", "mesoderm commitment", "mesoderm commitment", "haemogenic endothelium", "haemogenic endothelium", "haemogenic endothelium", "early haematopoietic stem cells", "early haematopoietic stem cells", "early haematopoietic stem cells", "late haematopoietic stem cells", "late haematopoietic stem cells", "early Natural Killer cells", "early Natural Killer cells", "early Natural Killer cells", "early Natural Killer cells", "late Natural Killer cells", "late Natural Killer cells", "late Natural Killer cells", "late Natural Killer cells"
      )
    )
  ) |>
  ggplot(
    aes(x = scRNA_CFH, y = `Nanog+ (%)`, colour = factor(day))
  ) + 
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  theme_catapult() +
  labs(
    x = "scRNAseq (standardised)",
    y = "CD45+/CD56+ productivity"
  ) +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1), nrow = 3)
  )
