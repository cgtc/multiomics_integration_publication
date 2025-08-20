process_olink <- function(data) {
  
  days_for_fresh_media <- data.frame(Timepoint = forcats::fct_inorder(c("Day 3", "Day 6", "Day 8", "Day 12", "Day 15", "Day 19", "Day 22", "Day 26", "Day 29", "Day 33", "Day 35")))
  
  fresh_media_data <- data |>
    dplyr::filter(Timepoint == "Fresh Media") |>
    dplyr::select(vessel, Assay, NPX) |>
    dplyr::mutate(linNPX = 2^NPX) |>
    dplyr::select(-NPX) |>
    dplyr::group_by(
      Assay, vessel          # accounts for the 2 assays measured twice because they're twice in the panel
    ) |>
    dplyr::summarise(
      linNPX = mean(linNPX, na.rm = TRUE)
    ) |>
    dplyr::group_by(Assay, vessel) |>
    dplyr::reframe(
      linNPX = linNPX,
      Timepoint = days_for_fresh_media$Timepoint
    ) |>  # Dummy linNPX for each of the days
    dplyr::arrange(Assay, vessel, Timepoint) |>
    tidyr::pivot_wider(names_from = vessel, values_from = linNPX) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      Media = dplyr::case_when(
        Timepoint == "Day 3" ~ `Step 1`,
        Timepoint == "Day 6" ~ `Step 2`,
        Timepoint == "Day 8" ~ `Step 3`,
        Timepoint == "Day 12" ~ sum(sum(`Step 3`, `Step 3b`)/2, `Step 3b`)/2,
        Timepoint == "Day 15" ~ `Step 4`,
        Timepoint == "Day 19" ~ sum(`Step 4`, `Step 4b`)/2,
        Timepoint == "Day 22" ~ sum(sum(`Step 4`, `Step 4b`)/2, `Step 4b`)/2,
        Timepoint == "Day 26" ~ sum(sum(sum(`Step 4`, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2,
        Timepoint == "Day 29" ~ sum(sum(sum(sum(`Step 4`, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2,
        Timepoint == "Day 33" ~ sum(sum(sum(sum(sum(`Step 4`, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2,
        Timepoint == "Day 35" ~ sum(sum(sum(sum(sum(sum(`Step 4`, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2, `Step 4b`)/2,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::mutate(
      NPX = log2(Media)
    ) |>
    dplyr::select(-(`Step 1`:`Media`)) |>
    dplyr::mutate(
      vessel = "Fresh Media"
    ) |>
    dplyr::select(vessel, Assay, Timepoint, NPX)
  
  
  m_data <- data |>
    dplyr::filter(
      DF == "Neat", stream == "1.5X", Timepoint %in% c("Day 15", "Day 19", "Day 22", "Day 26", "Day 29", "Day 33", "Day 35"), vessel != "Fresh Media"
    ) |>
    dplyr::select(Assay, Timepoint, NPX) |>
    dplyr::left_join(fresh_media_data, by = c("Assay", "Timepoint")) |>
    dplyr::group_by(Assay) |>
    tidyr::nest() |>
    dplyr::mutate(
     # model = purrr::map(data, ~ lm(NPX.y ~ NPX.x, data = .x)),
     # model2 = purrr::map(data, ~ lm(NPX.y ~ poly(NPX.x, 2), data = .x)),
     # model_AIC = purrr::map_dbl(model, ~ AIC(.x)),
     # model2_AIC = purrr::map_dbl(model2, ~ AIC(.x)),
     # delta_AIC = model_AIC - model2_AIC,
     # model_slope = purrr::map_dbl(model, ~ coef(.x)[2]),
     # model2_slope = purrr::map_dbl(model2, ~ coef(.x)[2]),
      cor_kendall = purrr::map_dbl(data, ~ cor(.x$NPX.x, .x$NPX.y, method = "kendall")),
      sum_delta = purrr::map_dbl(data, ~ sum(.x$NPX.y - .x$NPX.x)),
      spread_delta = purrr::map_dbl(data, ~ sd(.x$NPX.y - .x$NPX.x)),
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-data)
  
  # TYPE 1: Identify proteins that have similar trends to the media but different value
  m_data_type1 <- m_data |>
    dplyr::filter(
      spread_delta < quantile(spread_delta, 1/3),
      (sum_delta > quantile(sum_delta, 0.8, type = 8)) | (sum_delta < quantile(sum_delta, 0.2, type = 8)),
      abs(sum_delta) > 30
    )
  
  # TYPE 2: Identify proteins that have different trends to the media and different value
  m_data_type2 <- m_data |>
    dplyr::filter(
      spread_delta > quantile(spread_delta, 1/3),
      (sum_delta > quantile(sum_delta, 0.8, type = 8)) | (sum_delta < quantile(sum_delta, 0.2, type = 8)),
      abs(sum_delta) > 40,
      cor_kendall < -0.3 | cor_kendall > 0.3
    )
  
  # Selected assays - final list
  m_assays <- sort(unique(c(m_data_type2$Assay, m_data_type1$Assay)))
  
  # Check with PCA selected assays
  # proteins that don't contribute to the PCA variance explained but differ from the media trend somehow
  # m_odd_assays <- m_assays[-which(m_assays %in% optim_var_assays)]
  
  m_label <- data |> dplyr::select(UniProt, vessel, Timepoint) |>
    dplyr::mutate(
      dplyr::across(c(2, 3), \(x) as.character(x) |> tolower()),
      label = paste0("r4_v", stringr::str_extract(vessel, "v([1-8])", group = 1), "_d", stringr::str_extract(Timepoint, "day ([0-9]+)", group = 1)),
      label = ifelse(stringr::str_detect(label, "NA_dNA"), NA_character_, label)
    ) |>
    dplyr::pull(label)
  
  processed_data <- dplyr::bind_cols(data, m_label)
  names(processed_data)[36] <- "label" #ffs
  processed_data$omics <- "olink"
  
  # Return all data as a list
  list(
    m_fresh_media_data = fresh_media_data,
    m_data = m_data,
    m_data_type1 = m_data_type1,
    m_data_type2 = m_data_type2,
    m_assays = m_assays,
    data = processed_data
  )
}
