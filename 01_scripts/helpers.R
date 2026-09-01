#' Grab data from S3, GitHub, or local, prioritising local cache.
#'
#' @param file S3 link or local file name
#' @param store_path Path to store the file if it is downloaded from S3. 00_data is appended automagically.
#'
#' @returns Whatever you yoinked. Usually a dataframe from an RDS object.
#' @export
#' @noRd
fetch_rds <- function(file, store_path = file.path("raw", "unprocessed")) {
  if (!requireNamespace("fs", quietly = TRUE)) {
    install.packages("fs")
  }

  is_s3 <- FALSE
  is_gh <- FALSE
  f_name <- file

  if (stringr::str_detect(file, "^s3://")) {
    is_s3 <- TRUE
    f_name <- fs::path_file(file)
  }

  if (stringr::str_detect(file, "raw\\.githubusercontent\\.com")) {
    is_gh <- TRUE
    f_name <- fs::path_file(file)
  }

  if (tolower(fs::path_ext(f_name)) != "rds") stop("File must be an RDS file.", call. = FALSE)
  if ((is_s3 == TRUE) && (is_gh == TRUE)) stop("Internal parser error - File cannot be from both S3 and GitHub.", call. = FALSE)

  if (file.exists(file.path("00_data", store_path, f_name))) {
    cli::cat_rule("Reading file from local folder.")
    return(readRDS(file.path("00_data", store_path, f_name)))
  } 
  
  if (file.exists(file.path("00_data", "raw", f_name))) {
    cli::cat_rule("Reading processed data from standard folder.")
    return(readRDS(file.path("00_data", "raw", f_name)))
  } 
  
  if (is_s3) {
    cli::cat_rule("File not found locally - Downloading from S3")
    tryCatch(
      {
        aws.s3::save_object(file, file = file.path("00_data", store_path, f_name))
      },
      error = function(e) {
        stop("Error downloading file from S3: ", e)
      }
    )
    cli::cat_rule("Download Complete. Reading downloaded file.")
    return(readRDS(file.path("00_data", store_path, f_name)))
  } 
  
  if (is_gh) {
    cli::cat_rule("File not found locally - Downloading from GitHub")
    token <- fetch_git_creds()
    
    # Build headers with token if available
    headers <- if (is.null(token)) {
      c()
    } else {
      c(Authorization = paste("token", token))
    }
    
    tryCatch(
      {
        tmp <- readRDS(url(file, headers = headers))
      },
      error = function(e) {
        stop("Error downloading file from GitHub: ", e)
      }
    )
    cli::cat_rule("Download Complete. Saving downloaded file.")
    
    # Ensure directory exists
    dir.create(file.path("00_data", store_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(tmp, file.path("00_data", store_path, f_name))

    return(tmp)
  }
  
  stop("File not found locally nor on S3/GitHub.", call. = FALSE)
}


#' Grabs github token, wherever it may be
#'
#' @returns GitHub token
#' @export
#' @noRd
fetch_git_creds <- function() {
  gh_present <- TRUE
  
  if (!requireNamespace("gitcreds", quietly = TRUE)) {
    cli::cat_bullet("Gitcreds not installed. Falling back to GITHUB_PAT...")
    gh_present <- FALSE
  }

  if (Sys.getenv("GITHUB_PAT") == "" & !gh_present) {
    cli::cat_boxx("GITHUB_PAT not set. Please set the GITHUB_PAT environment variable.")
    return(NULL)
  } 
  
  if (Sys.getenv("GITHUB_PAT") != "") {
    cli::cat_bullet("GITHUB_PAT environment variable is set.")
    return(Sys.getenv("GITHUB_PAT"))
  } 
  
  if (gh_present) {
    cli::cat_bullet("Gitcreds installed.")

    if (!is.null(gitcreds::gitcreds_parse_output(gitcreds::gitcreds_get())[[4]])) {
      return(gitcreds::gitcreds_parse_output(gitcreds::gitcreds_get())[[4]])
    }
    
    cli::cat_boxx("No GitHub token found in gitcreds.")
    return(NULL)
  }
}
