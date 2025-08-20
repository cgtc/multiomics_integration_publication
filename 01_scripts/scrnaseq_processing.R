#' Creates a metadata object from a Seurat object, adds variable features to it, and returns it as a tibble.
#'
#' @param data Seurat object
#' @param .binding_id ID to bind the metadata to the data
#'
#' @returns A tibble with metadata + variable features normalised expression
#' @export
#' @noRd
create_bound_object <- function(data, .binding_id = "list_id") {
  if (is.list(data)) {
    n <- length(data)
    var_genes <- lapply(1:n, function(i) {
      data[[i]]@assays$RNA@meta.data$var.features[!data[[i]]@assays$RNA@meta.data$var.features |> is.na()]
    })
    
    # Getting the union of all variable genes across supplied runs.
    var_genes <- reduce(var_genes, union)
    
    metadata <- lapply(1:n, function(i) {
      create_metadata_object(data[[i]])
    })
    
    final <- lapply(1:n, function(x) {
      FetchData(data[[x]], vars = unique(var_genes), layer = "data") |> 
        rownames_to_column("barcode") |>
        as_tibble() |> 
        left_join(x = metadata[[x]], y = _) |>
        mutate(
          label = tolower(paste0(str_replace(run, "un ", ""), "_", str_replace(vessel, "0", ""), "_", day)),
          omics = "scRNA"
        )
    })
    
    return(bind_rows(final, .id = .binding_id))
  } else {
    var_genes <- data@assays$RNA@meta.data$var.features[!data@assays$RNA@meta.data$var.features |> is.na()]
    metadata <- create_metadata_object(data)
    final <- FetchData(data, vars = unique(var_genes), layer = "data") |> 
      rownames_to_column("barcode") |>
      as_tibble() |> 
      left_join(x = metadata, y = _) |>
      mutate(
        label = tolower(paste0(str_replace(run, "un ", ""), "_", str_replace(vessel, "0", ""), "_", day)),
        omics = "scRNA"
      )
    
    return(final)
  }
}



#' Creates a pseudobulk RNAseq tibble from a scRNAseq Seurat object by aggregating across specified groups.
#'
#' @param object Seurat object
#' @param groups Groups to aggregate by
#'
#' @returns A tibble with metadata + all features normalised expression
#' @export
#' @noRd
create_pseudobulk_object <- function(object, groups = c("day", "vessel")) {
  
  if (class(object) != "Seurat") stop("Object must be a Seurat object.")
  
  m_pb <- Seurat::AggregateExpression(
    object,
    return.seurat = TRUE,
    group.by = groups,
    assay = "RNA"
  )

  m_pb |>
    Seurat::FetchData(vars = rownames(m_pb)) |>
    tibble::rownames_to_column("barcode") |>
    tibble::as_tibble() 
}