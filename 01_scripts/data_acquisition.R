## Data acquisition script #
# 20/01/2025 #
# This script is designed to pull necessary processed data for performing integration ##
# "Raw" data in this context is the processed output from each -omics data set, that has been
# run using their given data analysis pipeline. For example, scRNA sequencing data has been processed
# using the Seurat pipeline.


# TODO - add data pipeline information for each -omics layer

#### Init ####
library(tidyverse)
library(here)
library(rstatix)
library(Seurat)
library(CGTCseq)
library(aws.s3)
library(conflicted)
library(DESeq2)
library(apeglm)
library(readxl)
library(aqua)
library(janitor)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(purrr::reduce)
conflicts_prefer(base::union)
conflicts_prefer(base::setdiff)
conflicts_prefer(dplyr::select)

setwd(here())

# Source helper functions
source("01_scripts/helpers.R")
source("01_scripts/scrnaseq_processing.R")
source("01_scripts/process_data_processing.R")

# Git credentials
token <- fetch_git_creds()
owner <- "cgtc"
branch <- "master"

# Data directories
if (!dir.exists("00_data")) dir.create("00_data")
if (!dir.exists("00_data/raw")) dir.create("00_data/raw")
if (!dir.exists("00_data/raw/processed")) dir.create("00_data/raw/processed")
if (!dir.exists("00_data/raw/unprocessed")) dir.create("00_data/raw/unprocessed")
if (!dir.exists("00_data/raw/unprocessed/lcms_processing_files")) dir.create("00_data/raw/unprocessed/lcms_processing_files")
if (!dir.exists("00_data/raw/unprocessed/miRNA_processing_files")) dir.create("00_data/raw/unprocessed/miRNA_processing_files")

#### scRNA sequencing data ####
# Loads processed scRNA sequencing data, following standard pipeline.
# Only selects top 3000 highly variable genes.
# s3://plasticell-transcriptomics/seurat/seurat/merged_seurat_normalised.rds
# Set AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY in your environment
# Also AWS_DEFAULT_REGION (important.)
# TODO: add local path for when running on the cloudy cloud

## Import from S3 (≈15GB) / change to local path if needed?
# cli::cat_rule("Downloading raw scRNA data (≈15GB, ≈38GB expanded in RAM).")
# merged_seurat_normalised <- fetch_rds("s3://plasticell-transcriptomics/seurat/seurat/merged_seurat_normalised.rds")
# scrna_variable_genes_df <- create_bound_object(merged_seurat_normalised)

## Annotations run 3/4
# cli::cat_rule("Downloading annotations.")
# annotations_scrnaseq_r3 <- fetch_rds("s3://plasticell-transcriptomics/seurat/seurat/singleR/run3_hpca_main_labels.rds")
# annotations_scrnaseq_r4 <- fetch_rds("s3://plasticell-transcriptomics/seurat/seurat/singleR/run4_hpca_main_labels.rds")

#### TODO: do we do unique(c(run3_variable_genes, run4_variable_genes)) or does each run take only its own variable genes?
## Grab variable genes and combine (both ways)
if (!file.exists("00_data/raw/combined_run3_run4_variable_genes_df.rds")) {
  cli::cat_rule("Downloading variable genes")
  run3_samples <- fetch_rds("s3://plasticell-transcriptomics/seurat/seurat/subsamples/run3_samples.rds")
  run4_samples <- fetch_rds("s3://plasticell-transcriptomics/seurat/seurat/subsamples/run4_combined_samples.rds")
  cli::cat_rule("Processing variable genes")

  ## Separate variable genes, then bind
  # run3_variable_genes_df <- create_bound_object(run3_samples)
  # run4_variable_genes_df <- create_bound_object(run4_samples)
  # combined_run3_run4_variable_genes_df <- bind_rows(run3_variable_genes_df, run4_variable_genes_df)

  # Union variable genes, then bind
  combined_run3_run4_variable_genes_union_df <- create_bound_object(list(run3_samples, run4_samples))

  saveRDS(combined_run3_run4_variable_genes_union_df, "00_data/raw/combined_run3_run4_variable_genes_df.rds")
} else {
  cli::cat_rule("File exists. Reading processed variable genes")

  # combined_run3_run4_variable_genes_df <- readRDS("00_data/raw/combined_run3_run4_variable_genes_df.rds")
  combined_run3_run4_variable_genes_union_df <- readRDS("00_data/raw/combined_run3_run4_variable_genes_df.rds")
}

#### scRNA pseudobulk data ####
# Loads processed scRNA pseudobulk data, following standard pipeline.
if (!file.exists("00_data/raw/combined_pseudobulk_df.rds")) {
  cli::cat_rule("Pseudobulk data object not found - processing scRNAseq into pseudobulk.")

  var_features <- setdiff(names(combined_run3_run4_variable_genes_union_df)[-c(1:47)],  
                              c("plasticell_sample_id", "label", "omics", "site"))

  pseudobulk_R3 <- create_pseudobulk_object(run3_samples, groups = c("day", "vessel")) |>
    mutate(
      omics = "scRNA",
      site = "CGTC",
      label = tolower(paste0("r3_v", str_extract(barcode, "V0([1-8])", group = 1), "_d", str_extract(barcode, "D([0-9]+)", group = 1)))
    ) |>
    select(label, omics, everything())


  pseudobulk_R4 <- create_pseudobulk_object(run4_samples, groups = c("day", "vessel")) |>
    mutate(
      omics = "scRNA",
      label = tolower(paste0("r4_v", str_extract(barcode, "V0([1-8])", group = 1), "_d", str_extract(barcode, "D([0-9]+)", group = 1)))
    ) |>
    select(label, omics, everything())

  combined_pseudobulk_df <- bind_rows(pseudobulk_R3, pseudobulk_R4)

  saveRDS(combined_pseudobulk_df, "00_data/raw/combined_pseudobulk_df.rds")
} else {
  cli::cat_rule("File exists. Reading processed pseudobulk data.")

  combined_pseudobulk_df <- readRDS("00_data/raw/combined_pseudobulk_df.rds")
}


#### loading lcms data ####
source(here("01_scripts/lcms_data_loading.R"))

## processing lcms data
if (!file.exists("00_data/raw/processed_lcms.rds")) {
  source("01_scripts/lcms_processing.R")

  cli::cat_rule("Processing lcms data")

  # Read the file
  processed_lcms <- process_lcms_data(
    raw_lcms_r3.1,
    raw_lcms_r3.2,
    raw_lcms_r3.3,
    raw_lcms_r4,
    r3_metabolites_to_remove
  )

  saveRDS(processed_lcms, "00_data/raw/processed_lcms.rds")
} else {
  cli::cat_rule("File exists. Reading processed lcms data")
  processed_lcms <- readRDS("00_data/raw/processed_lcms.rds")
}

#### O-link data ####
if (!file.exists("00_data/raw/olink_processed.rds")) {
  raw_olink <- fetch_rds("https://raw.githubusercontent.com/cgtc/olink_plasticell_cleanup/master/data/olink_data.rds")

  source("01_scripts/olink_processing.R")
  cli::cat_rule("Processing olink data")

  olink_processed <- process_olink(raw_olink)
  saveRDS(olink_processed, "00_data/raw/olink_processed.rds")
} else {
  cli::cat_rule("File exists. Reading processed olink data")

  processed_olink <- readRDS("00_data/raw/olink_processed.rds")
}

#### miRNA data ####
raw_miRNA <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_miRNA_reporting/master/data/multiomics_data/raw_miRNA_data.rds")

## loading unbiased miRNA IDs for filtering
uncorrelated_miRNA_ids <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_miRNA_reporting/master/data/multiomics_data/unbiased_miRNA_ids.rds", store_path = file.path("raw", "unprocessed", "miRNA_processing_files"))

## processing miRNA data
if (!file.exists("00_data/raw/processed_miRNA.rds")) {
  source("01_scripts/miRNA_processing.R")

  cli::cat_rule("Processing miRNA data")

  # Read the file
  processed_miRNA <- process_miRNA_data(
    raw_miRNA,
    uncorrelated_miRNA_ids
  )

  saveRDS(processed_miRNA, "00_data/raw/processed_miRNA.rds")
} else {
  cli::cat_rule("File exists. Reading processed miRNA data")
  processed_miRNA <- readRDS("00_data/raw/processed_miRNA.rds")
}

#### Process data ####
# Run 3
# >>> AGG FCM data: D26, V1 no aggregate data FCM due to low dissociated AGG count
# >>> BIOFLEX D12 V2 only gas was run, not chemistry = mistake
# >>> BIOFLEX D33 All vessels no Na+ Ca2+ K+
# >>> BIOFLEX D35 No data
# >>> Na+ Ca2+ K+ in same GDrive Plasticell Run3 folder in /DASbox
#
# Run 4
# >>> FCM D12 no day 12 FCM done based on low cell density, prioritised scRNAseq
# >>> BIOFLEX D0 all vessels all data
# >>> FCM D15 V1 SUS no duplicate data Not sure why? cell counts seemed to be comparable to other vessels
# >>> FCM D35 no extended NK panel Only used for

# raw_process_data <- list(
#   r3_process_data = read_excel("00_data/raw/processed/Plasticell Runs 3 and 4 process data summaries.xlsx", sheet = "Run 3"),
#   r4_process_data = read_excel("00_data/raw/processed/Plasticell Runs 3 and 4 process data summaries.xlsx", sheet = "Run 4"),
#   killing_assay_data = read_excel("00_data/raw/processed/Plasticell Runs 3 and 4 process data summaries.xlsx", sheet = "Killing assay")
# )
#
# saveRDS(raw_process_data, file = "raw_process_data.rds")

raw_process_data <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_flow_analysis/main/00_data/multiomics_process_data/raw_process_data.rds")

processed_process_data <- process_process_data(raw_process_data)
