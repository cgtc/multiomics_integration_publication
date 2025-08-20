removeZeroVar <- function(mat) {
  # Check variance for each column using apply
  zero_var_cols <- apply(mat, 2, function(x) min(x) == max(x))

  # Subset matrix to keep only columns with non-zero variance
  mat[, !zero_var_cols]
}

process_miRNA_data <- function(raw_miRNA, uncorrelated_miRNA_ids) {
  # filtering the matrix based on the unbiased miRNAs
  feature_counts_annotated_miRNA_low_counts <- raw_miRNA |>
    select(any_of(c("label", "vessel", "well", "sample_name", "day", "run", "file_name", uncorrelated_miRNA_ids)))

  # getting a matrix from the data frame
  feature_counts_matrix_miRNA_low_counts <- feature_counts_annotated_miRNA_low_counts |>
    select(!c(sample_name:file_name, label:vessel)) |>
    column_to_rownames(var = "well") |>
    removeZeroVar()

  # getting the metadata to accompany the matrix
  feature_counts_metadata_miRNA_low_counts <- feature_counts_annotated_miRNA_low_counts |>
    select(label:file_name)

  # getting the long version of the full dataset
  feature_counts_annotated_miRNA_long_low_counts <- feature_counts_annotated_miRNA_low_counts |>
    pivot_longer(
      cols = starts_with("ENS"), # Adjust as needed to select sample columns
      names_to = "gene_id",
      values_to = "count"
    )

  # finding miRNAs which don't have an average count over 3 on at least one day
  not_low_count_miRNAs <- feature_counts_annotated_miRNA_long_low_counts |>
    select(day, gene_id, count) |>
    group_by(gene_id, day) |>
    summarise(average_count = mean(count), .groups = "drop") |>
    group_by(gene_id) |>
    mutate(has_day_over_5 = any(average_count > 10)) |>
    ungroup() |>
    dplyr::filter(has_day_over_5 == TRUE) |>
    distinct(gene_id) |>
    pull(gene_id)

  # filtering the matrix based on the unbiased miRNAs
  feature_counts_annotated_miRNA <- feature_counts_annotated_miRNA_low_counts |>
    select(any_of(c("label", "vessel", "well", "sample_name", "day", "run", "file_name", not_low_count_miRNAs))) |>
    mutate(unique_sample_name = str_c(sample_name, day, sep = " day ")) |>
    relocate(unique_sample_name) |>
    # removing fresh media samples
    dplyr::filter(str_detect(unique_sample_name, "STR") | str_detect(unique_sample_name, "Stream")) |>
    dplyr::filter(!str_detect(unique_sample_name, "Fresh")) |>
    arrange(day) #|>
  # these days were originally removed to allow the DESeq transformation to take place (as these days are not present in both runs), however I instead changed the DESeq design from * to +
  # dplyr::filter(!day == 8) |>
  # dplyr::filter(!day == 22) |>
  # dplyr::filter(!day == 29)

  # getting a matrix from the data frame
  feature_counts_matrix_miRNA <- feature_counts_annotated_miRNA |>
    column_to_rownames(var = "unique_sample_name") |>
    select(!label:file_name) |>
    removeZeroVar()

  # getting the metadata to accompany the matrix
  feature_counts_metadata_miRNA <- feature_counts_annotated_miRNA |>
    select(unique_sample_name:file_name) |>
    column_to_rownames(var = "unique_sample_name")

  # getting the long version of the full dataset
  feature_counts_annotated_miRNA_long <- feature_counts_annotated_miRNA |>
    pivot_longer(
      cols = starts_with("ENS"), # Adjust as needed to select sample columns
      names_to = "gene_id",
      values_to = "count"
    )

  # creating deseq object
  # this is changed from the combined analysis by replacing the day * run with day + run to account for different sampling days between runs
  dds_data_day_6v3 <- DESeqDataSetFromMatrix(
    countData = t(feature_counts_matrix_miRNA),
    colData = feature_counts_metadata_miRNA |>
      mutate(day = as.factor(day)) |>
      mutate(day = relevel(as.factor(day), ref = "3")),
    design = ~ day + run
  )

  dds_day_6v3 <- DESeq(dds_data_day_6v3)

  resLFC_day_6v3 <- lfcShrink(dds_day_6v3, coef = "day_6_vs_3", type = "apeglm")

  normalised_counts_runs <- counts(dds_day_6v3, normalized = TRUE)

  vst_day_6v3 <- varianceStabilizingTransformation(dds_day_6v3, blind = FALSE)

  vst_counts_by_day <- assay(vst_day_6v3)

  feature_counts_annotated_miRNA_long_transformed <- vst_counts_by_day |>
    t() |>
    as.data.frame() |>
    rownames_to_column(var = "unique_sample_name") |>
    as_tibble() |>
    pivot_longer(
      cols = starts_with("ENS"), # Adjust as needed to select sample columns
      names_to = "gene_id",
      values_to = "count"
    ) |>
    left_join(feature_counts_metadata_miRNA |>
      rownames_to_column(var = "unique_sample_name"))

  feature_counts_matrix_miRNA_transformed <- vst_counts_by_day |>
    t() |>
    as.data.frame()

  processed_miRNA_wide <- feature_counts_annotated_miRNA_long_transformed |>
    pivot_wider(
      names_from = gene_id,
      values_from = count
    ) |>
    select(-unique_sample_name, -well, -file_name)

  return(list(
    processed_miRNA = processed_miRNA_wide
  ))
}
