## Run 3 Parts 1-2-3
raw_lcms_r3.1 <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_lcms_analysis/main/data/run_3/processed/pac_run1_processed.rds")
raw_lcms_r3.2 <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_lcms_analysis/main/data/run_3/processed/pac_run2_processed.rds")
raw_lcms_r3.3 <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_lcms_analysis/main/data/run_3/processed/pac_run3_processed.rds")

## Run 4
raw_lcms_r4 <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_lcms_analysis/main/data/run_4/multiomics/run_4_raw_data.rds")

## Run 3 metabolites to remove list
r3_metabolites_to_remove <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_lcms_analysis/main/data/run_3/run_3_metabs_to_remove.rds", store_path = file.path("raw", "unprocessed", "lcms_processing_files"))

## Run 4 produced vs supplemented list
r4_lcms_produced_vs_supplemented <- fetch_rds("https://raw.githubusercontent.com/cgtc/plasticell_lcms_analysis/main/data/run_4/multiomics/r4_lcms_produced_vs_supplemented.rds", store_path = file.path("raw", "unprocessed", "lcms_processing_files"))
