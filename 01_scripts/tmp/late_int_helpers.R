#### Function Definitions ####
# Function to load and process pseudobulk data if not already saved
load_and_process_pseudobulk <- function() {
  if (!file.exists(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_pseudobulk.rds"))) {
    # Load Seurat object
    run4_combined_samples <- readRDS(file.path("/efs", "plasticell", "scRNA", "plasticell_scRNA",
                                               "00_data", "suerat", "subsamples", "run_4_combined_samples.rds"))
    
    # Aggregate expression and extract variable genes
    pseudo_seurat <- AggregateExpression(run4_combined_samples, assays = "RNA", return.seurat = TRUE, group.by = c("vessel", "day"))
    pseudo_counts <- GetAssayData(pseudo_seurat, slot = "data")
    variable_genes <- VariableFeatures(run4_combined_samples)
    pseudo_counts_vg <- pseudo_counts[variable_genes, ]
    
    # Prepare long-format data and merge with NK proportions
    pseudo_counts_vg_long <- pseudo_counts_vg %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      as_tibble() %>%
      pivot_longer(cols = -gene, names_to = "vessel_day", values_to = "expression") %>%
      left_join(proportion_NK %>% select(vessel_day, proportion_NK))
    
    # Save processed data
    saveRDS(pseudo_counts, file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_pseudobulk.rds"))
    saveRDS(pseudo_counts_vg, file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_pseudobulk_vg.rds"))
    saveRDS(variable_genes, file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_variable_genes.rds"))
    saveRDS(pseudo_counts_vg_long, file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_pseudobulk_long.rds"))
  }else{
    pseudo_counts <<- readRDS(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_pseudobulk.rds"))
    pseudo_counts_vg <<- readRDS(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_pseudobulk_vg.rds"))
    variable_genes <<- readRDS(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_variable_genes.rds"))
    pseudo_counts_vg_long <<- readRDS(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_pseudobulk_long.rds"))
  }
}

# Function to generate NK proportion dataframe if not already saved
generate_nk_proportion <- function() {
  if (!file.exists(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_nk_proportion.rds"))) {
    # Generate NK proportions from metadata
    proportion_NK <- run4_combined_samples@meta.data %>%
      mutate(vessel_day = paste0(vessel, "_", day)) %>%
      group_by(vessel_day) %>%
      mutate(total_cells = n(),
             nk_cell_count = sum(hpca.main == "NK_cell"),
             proportion_NK = nk_cell_count / total_cells) %>%
      distinct(nk_cell_count, total_cells, proportion_NK)
    
    # Save the proportion data
    saveRDS(proportion_NK, file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_nk_proportion.rds"))
  } else {
    proportion_NK <<- readRDS(file.path("00_data", "raw", "unprocessed", "scRNA_processing_files", "run4_nk_proportion.rds"))
  }
}

# Function to generate GLM Gaussian LASSO models if not already saved
generate_glm_lasso <- function(pseudo_counts_vg, proportion_NK, vessel_exclude) {
  set.seed(123)
  
  # Generate train/test datasets for vessel 1 excluded and vessel 1 only
  vessels_train <- proportion_NK %>% dplyr::filter(!str_detect(vessel_day, paste0("V0", vessel_exclude))) %>% pull(vessel_day)
  NK_proportion_train <- proportion_NK %>% dplyr::filter(vessel_day %in% vessels_train) %>% arrange(match(vessel_day, colnames(pseudo_counts_vg))) %>% column_to_rownames("vessel_day")
  counts_train <- t(pseudo_counts_vg[, match(rownames(NK_proportion_train), colnames(pseudo_counts_vg))])
  
  vessels_test <- proportion_NK %>% dplyr::filter(str_detect(vessel_day, paste0("V0", vessel_exclude))) %>% pull(vessel_day)
  NK_proportion_test <- proportion_NK %>% dplyr::filter(vessel_day %in% vessels_test) %>% arrange(match(vessel_day, colnames(pseudo_counts_vg))) %>% column_to_rownames("vessel_day")
  counts_test <- t(pseudo_counts_vg[, match(rownames(NK_proportion_test), colnames(pseudo_counts_vg))])
  
  
  if(all(rownames(counts_train) == rownames(NK_proportion_train))){
    print("Train data is correctly matched")
    glm_cv <- cv.glmnet(
      x = as.matrix(counts_train),
      y = as.numeric(NK_proportion_train$proportion_NK),
      alpha = 1, 
      family = "gaussian",
      nfolds = 10,
      standardize = FALSE
    )
  }
  
  return(list(glm_cv = glm_cv,
              vessels_train = vessels_train,
              NK_proportion_train = NK_proportion_train,
              counts_train = counts_train,
              vessels_test = vessels_test,
              NK_proportion_test = NK_proportion_test,
              counts_test = counts_test))
}

# Function to evaluate the model's performance
evaluate_model <- function(glm_cv, counts_test, NK_proportion_test) {
  predictions <- predict(glm_cv, newx = as.matrix(counts_test), s = "lambda.min")
  actuals <- as.numeric(NK_proportion_test$proportion_NK)
  
  # R-squared
  r_squared <- cor(actuals, predictions)^2
  
  # RMSE and MAE
  rmse <- sqrt(mean((actuals - predictions)^2))
  mae <- mean(abs(actuals - predictions))
  
  return(list(r_squared = r_squared, rmse = rmse, mae = mae))
}

# Function to extract and sort feature importance
get_feature_importance <- function(glm_cv, non_zero = TRUE) {
  coef_matrix <- coef(glm_cv, s = "lambda.min")
  coef_df <- as.data.frame(as.matrix(coef_matrix))
  colnames(coef_df) <- c("Coefficient")
  coef_df$Feature <- rownames(coef_df)
  
  coef_sorted <- coef_df %>%
    dplyr::filter(Feature != "(Intercept)") %>%
    arrange(desc(abs(Coefficient)))
  
  if(non_zero){
    feature_importance <- coef_sorted %>%
      dplyr::filter(Coefficient != 0) %>%
      arrange(desc(Coefficient))
  }else{
    feature_importance <- coef_sorted %>%
      arrange(desc(Coefficient))
  }
  
  return(feature_importance)
}
