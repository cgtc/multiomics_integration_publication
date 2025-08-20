process_process_data <- function(process_data_list) {
  ### PROCESS DATA PROCESSING
  r3_process_data_wide <- raw_process_data$r3_process_data |>
    pivot_longer(
      cols = contains("..."),
      names_to = "variable",
      values_to = "flow_value"
    ) |>
    mutate(base_name = str_extract(variable, "^[^...]+")) |>
    group_by(Label, `Cell portion`, base_name) |>
    mutate(mean_percentage = mean(flow_value, na.rm = TRUE)) |>
    ungroup() |>
    select(-variable, -flow_value) |>
    group_by(Label, `Cell portion`, base_name) |>
    slice_head(n = 1) |>
    ungroup() |>
    pivot_wider(
      names_from = base_name,
      values_from = mean_percentage
    ) |>
    mutate(
      vessel = as.numeric(str_extract(Label, "(?<=v)\\d+")),
      day = as.numeric(str_extract(Label, "(?<=d)\\d+")),
      run = as.numeric(str_extract(Label, "(?<=r)\\d+"))
    ) |>
    relocate(run, day, vessel) |>
    rename(
      aggregates = "Aggregates (%)",
      viability = "Viability (%)"
    ) |>
    mutate(across(ends_with("(%)"), ~ . * 100)) |>
    select(!Run)

  r4_process_data_wide <- raw_process_data$r4_process_data |>
    pivot_longer(
      cols = contains("..."),
      names_to = "variable",
      values_to = "flow_value"
    ) |>
    mutate(base_name = str_extract(variable, "^[^...]+")) |>
    group_by(Label, `Cell portion`, base_name) |>
    mutate(mean_percentage = mean(flow_value, na.rm = TRUE)) |>
    ungroup() |>
    select(-variable, -flow_value) |>
    group_by(Label, `Cell portion`, base_name) |>
    slice_head(n = 1) |>
    ungroup() |>
    pivot_wider(
      names_from = base_name,
      values_from = mean_percentage
    ) |>
    mutate(
      vessel = as.numeric(str_extract(Label, "(?<=v)\\d+")),
      day = as.numeric(str_extract(Label, "(?<=d)\\d+")),
      run = as.numeric(str_extract(Label, "(?<=r)\\d+"))
    ) |>
    relocate(run, day, vessel) |>
    rename(
      aggregates = "Aggregates (%)",
      viability = "Viability (%)"
    ) |>
    mutate(across(ends_with("(%)"), ~ . * 100)) |>
    select(!Run)

  r3_suspension_data <- r3_process_data_wide |>
    filter(`Cell portion` == "Suspension") |>
    arrange(day, vessel) |>
    select(!`Vessel type`:Day)

  r4_suspension_data <- r4_process_data_wide |>
    filter(`Cell portion` == "Suspension") |>
    arrange(day, vessel) |>
    select(!`Vessel type`:Day)

  unscaled_suspension_process_data <- r3_suspension_data |>
    rbind(r4_suspension_data) |>
    mutate(omics = "process") |>
    relocate(omics) |>
    rename(label = "Label")

  suspension_process_data_long_unscaled <- unscaled_suspension_process_data |>
    mutate(across(viability:`Tra1-60+ (%)`, ~ gsub("\\(.*?\\)", "", .))) |>
    mutate(across(viability:`Tra1-60+ (%)`, ~ gsub("([0-9])\\*", "\\1", .))) |>
    mutate(across(viability:`Tra1-60+ (%)`, ~ ifelse(is.nan(.), NA, .))) |>
    mutate(across(viability:`Tra1-60+ (%)`, as.numeric)) |>
    pivot_longer(
      cols = viability:`Tra1-60+ (%)`,
      values_to = "value",
      names_to = "feature"
    ) |>
    group_by(feature) |>
    filter(
      !(str_detect(feature, "%") & (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) <= 10)
    ) |>
    ungroup()

  suspension_process_data_long <- suspension_process_data_long_unscaled |>
    group_by(feature) |>
    mutate(
      scaled_value = scale(value),
      minmax_value = rescale(value)
    ) |>
    ungroup()

  suspension_process_data <- suspension_process_data_long |>
    select(-minmax_value, -scaled_value) |>
    pivot_wider(
      names_from = feature,
      values_from = value
    )


  ### CYTOTOXICITY DATA PROCESSING

  cytotoxicity_process_data <- raw_process_data$killing_assay_data |>
    pivot_longer(
      cols = contains("..."),
      names_to = "variable",
      values_to = "cytotoxicity_value"
    ) |>
    mutate(ratio = str_extract(variable, "^[^...]+")) |>
    select(-variable) |>
    clean_names() |>
    group_by(run, label, panel, ratio) |>
    summarise(
      cytotoxicity_value = mean(cytotoxicity_value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    ungroup() |>
    mutate(
      cell_line = str_extract(panel, "^[^ ]+"),
      timepoint = str_extract(panel, "(?<=\\s)(\\d+)")
    ) |>
    select(-panel) |>
    mutate(
      vessel = as.numeric(str_extract(label, "(?<=v)\\d+")),
      day = as.numeric(str_extract(label, "(?<=d)\\d+")),
      omics = "process"
    ) |>
    relocate(label, run, vessel, day, timepoint, cell_line) |>
    mutate(inverted_cytotoxicity = 100 - cytotoxicity_value) |>
    filter(
      !ratio == "Cancer only",
      !timepoint == 24
    ) |>
    mutate(cancer_cells_killed_per_NK = case_when(
      ratio == "1:1" ~ inverted_cytotoxicity / 100,
      ratio == "0:5:1" ~ (inverted_cytotoxicity / 100) * 2,
      ratio == "0:25:1" ~ (inverted_cytotoxicity / 100) * 4
    ))

  cytotoxicity_process_data_wide <- cytotoxicity_process_data |>
    pivot_wider(
      names_from = c(cell_line, timepoint, ratio),
      values_from = cancer_cells_killed_per_NK
    )

  cytotoxicity_process_data_summarised <- cytotoxicity_process_data |>
    group_by(label, run, vessel, day, cell_line, omics, timepoint) |>
    summarise(
      cancer_cells_killed_per_NK = mean(cancer_cells_killed_per_NK),
      ratio = "combined"
    ) |>
    ungroup()

  cytotoxicity_process_data_summarised_wide <- cytotoxicity_process_data_summarised |>
    pivot_wider(
      names_from = c(cell_line, timepoint, ratio),
      values_from = cancer_cells_killed_per_NK
    ) |>
    mutate(combined_cytotoxicity = rowMeans(across(c(A549_48_combined, PANC1_48_combined, T47D_48_combined)), na.rm = TRUE))

  cytotoxicity_process_data_summarised_wide_rescaled <- cytotoxicity_process_data_summarised |>
    pivot_wider(
      names_from = c(cell_line, timepoint, ratio),
      values_from = cancer_cells_killed_per_NK
    ) |>
    mutate(combined_cytotoxicity = rowMeans(across(c(A549_48_combined, PANC1_48_combined, T47D_48_combined)), na.rm = TRUE)) |>
    mutate(across(!label:omics, rescale))

  all_process_data_unscaled <- suspension_process_data |>
    full_join(cytotoxicity_process_data_summarised_wide,
      by = c("omics", "run", "day", "vessel", "label")
    )

  all_process_data_scaled <- suspension_process_data |>
    full_join(cytotoxicity_process_data_summarised_wide,
      by = c("omics", "run", "day", "vessel", "label")
    ) |>
    mutate(across(!label:omics, rescale))

  all_process_data_scaled_added_metrics <- all_process_data_scaled |>
    mutate(
      productivity_plus_a549_cytotoxicity = `CD45+/CD56+ productivity` + A549_48_combined,
      productivity_plus_panc1_cytotoxicity = `CD45+/CD56+ productivity` + PANC1_48_combined,
      productivity_plus_t47d_cytotoxicity = `CD45+/CD56+ productivity` + T47D_48_combined,
      productivity_plus_combined_cytotoxicity = `CD45+/CD56+ productivity` + combined_cytotoxicity
    )

  all_process_data_scaled_added_metrics_longer <- all_process_data_scaled_added_metrics |>
    pivot_longer(
      cols = !omics:label,
      names_to = "feature_name"
    )

  all_process_data_unscaled_added_metrics <- all_process_data_unscaled |>
    mutate(cancer_cells_killed_per_PSC = `CD45+/CD56+ productivity` * combined_cytotoxicity)

  all_process_data_unscaled_added_metrics_longer <- all_process_data_unscaled_added_metrics |>
    pivot_longer(
      cols = !omics:label,
      names_to = "feature_name"
    )

  return(list(
    suspension_process_data = suspension_process_data,
    all_process_data_unscaled_added_metrics = all_process_data_unscaled_added_metrics,
    cytotoxicity_process_data_summarised_wide = cytotoxicity_process_data_summarised_wide,
    all_process_data_scaled_added_metrics_longer = all_process_data_scaled_added_metrics_longer,
    all_process_data_scaled_added_metrics = all_process_data_scaled_added_metrics,
    all_process_data_unscaled_added_metrics_longer = all_process_data_unscaled_added_metrics_longer,
    all_process_data_unscaled_added_metrics = all_process_data_unscaled_added_metrics
  ))
}
