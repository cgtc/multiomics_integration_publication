process_lcms_data <- function(raw_lcms_r3.1, raw_lcms_r3.2, raw_lcms_r3.3, raw_lcms_r4, r3_metabolites_to_remove) {
  # Run 3 data processing
  # running data loading and wrangling that Chris has written

  #### Data tidying ####
  run_3_plotting_df_unfiltered <- bind_rows(raw_lcms_r3.1, raw_lcms_r3.2, raw_lcms_r3.3) |>
    filter(sample_type == "sample", str_detect(step, "D") & step != "DMEM") |>
    mutate(day = as.numeric(str_replace(step, "D", ""))) |>
    group_by(metabolite, day, run) |>
    mutate(
      avg_normalised_auc = mean(normalised_peak_area, na.rm = TRUE),
      mean_rt_deviation = mean(retention_time_deviation, na.rm = TRUE),
      std_normalised_auc = sd(retention_time_deviation, na.rm = TRUE)
    ) |>
    ungroup() |>
    group_by(metabolite, run) |>
    mutate(
      avg_scaled = scale(avg_normalised_auc),
      sd_scaled = scale(std_normalised_auc)
    ) |>
    ungroup()

  run_3_plotting_df <- run_3_plotting_df_unfiltered |>
    filter(!tolower(metabolite) %in% tolower(r3_metabolites_to_remove))

  ## Time series data ----
  ## Media
  run_3_media_bound_df <- bind_rows(raw_lcms_r3.1, raw_lcms_r3.2, raw_lcms_r3.3) |>
    filter(
      sample_type == "sample",
      str_detect(step, "Step"),
      !str_detect(step, fixed("+"))
    ) |>
    select(sample_name, dilution, stream, metabolite, step, normalised_peak_area, system, run) |>
    left_join(tibble(
      step = c("Step1", "Step2", "Step3", "Step3b", "Step4", "Step4b"),
      day = c(0, 3, 6, 8, 12, 15)
    )) |>
    mutate(
      ctrl_group = ifelse(str_detect(stream, "1"), "1X", "5X"),
      sample_type = "media"
    ) |>
    filter(!tolower(metabolite) %in% tolower(r3_metabolites_to_remove))

  run_3_max_day_media_bound_df <- run_3_media_bound_df |>
    filter(day == 15) |>
    mutate(day = 35) |>
    filter(!tolower(metabolite) %in% tolower(r3_metabolites_to_remove))

  run_3_media_bound_df <- run_3_media_bound_df |>
    full_join(run_3_max_day_media_bound_df) |>
    filter(!tolower(metabolite) %in% tolower(r3_metabolites_to_remove))

  ## Samples
  run_3_metabolite_bound_df <- bind_rows(raw_lcms_r3.1, raw_lcms_r3.2, raw_lcms_r3.3) |>
    filter(sample_type == "sample", str_detect(step, "D") & step != "DMEM") |>
    select(sample_name, dilution, stream, metabolite, step, normalised_peak_area, system, run) |>
    mutate(
      day = as.numeric(str_replace(step, "D", "")),
      ctrl_group = ifelse(str_detect(stream, "A|V1|V2|C"), "1X", "5X"),
      sample_type = "sample"
    ) |>
    filter(!tolower(metabolite) %in% tolower(r3_metabolites_to_remove))

  ## Combining, scaling and averaging
  run_3_combined_scaled_df <- bind_rows(run_3_metabolite_bound_df, run_3_media_bound_df) |>
    mutate(normalised_peak_area = ifelse(is.na(normalised_peak_area), 0, normalised_peak_area)) |>
    group_by(metabolite, ctrl_group, run) |>
    mutate(
      run_scaled = scale(normalised_peak_area)[, 1],
      n = n()
    ) |>
    filter(!tolower(metabolite) %in% tolower(r3_metabolites_to_remove))

  run_3_combined_scaled_df_summarised <- run_3_combined_scaled_df |>
    mutate(stream = case_when(
      # adjusted this based on Chris' scRNA sample naming
      str_detect(stream, "A") ~ "V5",
      str_detect(stream, "C") ~ "V6",
      TRUE ~ stream
    )) |>
    group_by(metabolite, day, stream, sample_type) |>
    summarise(
      avg_normalised_peak_area = mean(normalised_peak_area, na.rm = TRUE),
      sd_normalised_peak_area = sd(normalised_peak_area, na.rm = TRUE),
      sample_type = unique(sample_type),
      avg_scaled_normalised_peak_area = mean(run_scaled, na.RM = T),
      sd_scaled_normalised_peak_area = sd(run_scaled, na.rm = T)
    ) |>
    ungroup() |>
    filter(!tolower(metabolite) %in% tolower(r3_metabolites_to_remove))

  #### Run 4 processing ####

  #### processing following the original script ####

  run_4_data <- raw_lcms_r4 |>
    # remove QC samples
    filter(!is.na(sample_type)) |>
    mutate(area_ratio = as.numeric(area_ratio))

  run_4_data_filtered <- run_4_data |>
    # remove outlier
    filter(!sample_name == "X7_S_D35 V3") |>
    group_by(metabolite_name) |>
    # remove zero variance metabolites
    filter(!all(is.na(area_ratio) | area_ratio == 0)) |>
    ungroup() |>
    # remove internal standard
    filter(metabolite_name != "2-Isopropylmalic acid") |>
    arrange(metabolite_name, batch_name, day) |>
    mutate(stream = as_factor(stream))

  # Create a dataframe with the days you want to add for each metabolite_name
  new_fresh_media_days <- data.frame(day = c(10, 19, 22, 26, 29, 33, 35))

  fresh_samples_recalculated_added_rows <- run_4_data_filtered |>
    filter(sample_type == "F") |>
    group_by(metabolite_name) |>
    mutate(
      day = rep(new_fresh_media_days$day, length.out = n()),
      area_ratio = NA,
      area_exchanges = NA,
      ISTD_Area = NA,
      area = NA,
      LOD_alert = NA,
      sample_name = NA
    ) |>
    ungroup()

  fresh_samples_recalculated_initial <- run_4_data_filtered |>
    filter(sample_type == "F") |>
    filter(!day == 35) |>
    mutate(area_exchanges = area_ratio) |>
    relocate(metabolite_name, day, area_ratio, area_exchanges) |>
    arrange(metabolite_name, day)


  fresh_samples_recalculated_full <- rbind(fresh_samples_recalculated_initial, fresh_samples_recalculated_added_rows) |>
    relocate(metabolite_name, day, area_ratio, area_exchanges) |>
    arrange(metabolite_name, day)


  fresh_samples_recalculated <- fresh_samples_recalculated_full |>
    group_by(metabolite_name) |>
    mutate(
      area_exchanges = case_when(
        day == 8 ~ mean(area_ratio[day %in% c(6, 8)]),
        day == 15 ~ mean(area_ratio[day %in% c(12, 15)]),
        TRUE ~ area_ratio
      )
    ) |>
    mutate(
      area_exchanges = case_when(
        day == 10 ~ mean(c(area_exchanges[day == 8], area_ratio[day == 8])),
        TRUE ~ area_exchanges
      )
    ) |>
    mutate(
      area_exchanges = case_when(
        day == 19 ~ mean(c(area_exchanges[day == 15], area_ratio[day == 15])),
        TRUE ~ area_exchanges
      )
    ) |>
    mutate(
      area_exchanges = case_when(
        day == 22 ~ mean(c(area_exchanges[day == 19], area_ratio[day == 15])),
        TRUE ~ area_exchanges
      )
    ) |>
    mutate(
      area_exchanges = case_when(
        day == 26 ~ mean(c(area_exchanges[day == 22], area_ratio[day == 15])),
        TRUE ~ area_exchanges
      )
    ) |>
    mutate(
      area_exchanges = case_when(
        day == 29 ~ mean(c(area_exchanges[day == 26], area_ratio[day == 15])),
        TRUE ~ area_exchanges
      )
    ) |>
    mutate(
      area_exchanges = case_when(
        day == 33 ~ mean(c(area_exchanges[day == 29], area_ratio[day == 15])),
        TRUE ~ area_exchanges
      )
    ) |>
    filter(!day == 35)


  adj_fresh_media_time_plot_1.5 <- fresh_samples_recalculated |>
    mutate(stream = "1.5X") |>
    group_by(stream, sample_type, metabolite_name, day) |>
    summarise(
      avg_area = mean(area_exchanges),
      avg_sd = sd(area_exchanges)
    )

  adj_fresh_media_time_plot_2.25 <- fresh_samples_recalculated |>
    mutate(stream = "2.25X") |>
    group_by(stream, sample_type, metabolite_name, day) |>
    summarise(
      avg_area = mean(area_exchanges),
      avg_sd = sd(area_exchanges)
    )

  adj_fresh_media_time_plot_3 <- fresh_samples_recalculated |>
    mutate(stream = "3X") |>
    group_by(stream, sample_type, metabolite_name, day) |>
    summarise(
      avg_area = mean(area_exchanges),
      avg_sd = sd(area_exchanges)
    )

  adj_fresh_media_time_plot <- rbind(
    adj_fresh_media_time_plot_1.5,
    adj_fresh_media_time_plot_2.25,
    adj_fresh_media_time_plot_3
  )

  # remove standard (QC and fresh media samples already removed)
  spent_media_time_plot <- run_4_data_filtered |>
    # select spent media samples
    filter(!is.na(stream)) |>
    group_by(stream, sample_type, metabolite_name, day) |>
    summarise(
      avg_area = mean(area_ratio),
      avg_sd = sd(area_ratio)
    )

  adj_stream_temporal_plot_data <- rbind(
    spent_media_time_plot,
    adj_fresh_media_time_plot
  )

  produced_vs_supplemented <- r4_lcms_produced_vs_supplemented

  annotated_adj_stream_temporal_plot_data <- adj_stream_temporal_plot_data |>
    left_join(produced_vs_supplemented, by = "metabolite_name") |>
    mutate(avg_area = if_else(sample_type == "F" & Haem == "produced" & NK == "produced", 0, avg_area))

  annotated_adj_wider_stream_temporal_data <-
    annotated_adj_stream_temporal_plot_data |>
    select(!avg_sd) |>
    pivot_wider(
      names_from = sample_type,
      values_from = avg_area,
      values_fill = list(avg_area = NA)
    ) |>
    arrange(metabolite_name, day)

  fresh_spent_correlations_original <- annotated_adj_wider_stream_temporal_data |>
    filter(!is.na(`F`)) |>
    group_by(metabolite_name) |>
    cor_test(vars = `S`, vars2 = `F`, method = "kendall")

  fresh_spent_correlations_adj_p_values <- fresh_spent_correlations_original |>
    select(p) |>
    deframe() |>
    p.adjust(method = "BH")

  fresh_spent_correlations <- fresh_spent_correlations_original |>
    mutate(p = fresh_spent_correlations_adj_p_values)

  weak_fresh_spent_correlations <- fresh_spent_correlations |>
    filter(cor <= 0.4 | is.na(cor)) |>
    arrange(p)

  run_4_data_vessel_level <- run_4_data_filtered |>
    # select spent media samples
    filter(!is.na(stream)) |>
    group_by(vessel, sample_type, metabolite_name, day) |>
    summarise(
      avg_area = mean(area_ratio),
      avg_sd = sd(area_ratio)
    )

  fresh_and_spent_data <- rbind(
    run_4_data_vessel_level,
    adj_fresh_media_time_plot
  )

  annotated_fresh_and_spent_data <- fresh_and_spent_data |>
    left_join(produced_vs_supplemented, by = "metabolite_name") |>
    mutate(avg_area = if_else(sample_type == "F" & Haem == "produced" & NK == "produced", 0, avg_area))

  annotated_fresh_and_spent_data_wider <-
    annotated_fresh_and_spent_data |>
    select(!avg_sd) |>
    pivot_wider(
      names_from = sample_type,
      values_from = avg_area,
      values_fill = list(avg_area = NA)
    ) |>
    arrange(metabolite_name, day)

  run_4_data_fully_filtered <- annotated_fresh_and_spent_data |>
    filter(metabolite_name %in% weak_fresh_spent_correlations$metabolite_name) |>
    select(!Haem:NK) |>
    select(!avg_sd)

  run_3_data_fully_formatted <- run_3_combined_scaled_df_summarised |>
    mutate(run = 3) |>
    filter(metabolite %in% run_4_data_fully_filtered$metabolite_name) |>
    select(!sd_scaled_normalised_peak_area) |>
    rename(vessel = stream)

  run_4_data_fully_formatted <- run_4_data_fully_filtered |>
    mutate(run = 4) |>
    rename(
      metabolite = metabolite_name,
      avg_normalised_peak_area = avg_area
    ) |>
    group_by(metabolite) |>
    mutate(avg_scaled_normalised_peak_area = scale(avg_normalised_peak_area)[, 1]) |>
    ungroup() |>
    mutate(sample_type = case_when(
      sample_type == "F" ~ "media",
      sample_type == "S" ~ "sample"
    )) |>
    mutate(vessel = case_when(
      str_detect(stream, "1.5X") ~ "1.5X",
      str_detect(stream, "2.25X") ~ "2.25X",
      str_detect(stream, "3X") ~ "3X",
      TRUE ~ as.character(vessel)
    )) |>
    select(!stream)

  both_runs_data <- bind_rows(
    run_3_data_fully_formatted |>
      mutate(vessel = str_extract(vessel, "\\d+")),
    run_4_data_fully_formatted
  ) |>
    relocate(run, metabolite, sample_type, day, vessel) |>
    filter(sample_type == "sample") |>
    mutate(label = paste("r", run, "_v", vessel, "_d", day, sep = "")) |>
    relocate(label) |>
    select(!sd_normalised_peak_area)

  # create matrix with unique row names
  both_runs_data_wide <- both_runs_data |>
    select(label:vessel, avg_scaled_normalised_peak_area) |>
    pivot_wider(
      names_from = metabolite,
      values_from = avg_scaled_normalised_peak_area
    ) |>
    mutate(label = paste("r", run, "_v", vessel, "_d", day, sep = "")) |>
    relocate(label)

  return(
    list(
      both_runs_data_wide = both_runs_data_wide,
      both_runs_data = both_runs_data
    )
  )
}
