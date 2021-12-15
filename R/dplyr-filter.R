#' @method filter pathway_database
#' @importFrom rlang quos !!!
#' @importFrom dplyr filter
#' @export
filter.pathway_database <- function(.data, ..., .preserve = FALSE) {
  dots <- quos(...)
  # browser()
  if (length(.data@gene_list) == 0) {
    .data@gene_list = vector(mode = "list", length = length(.data@pathway_id)) %>%
      purrr::map(function(x) {
        x = data.frame()
        x
      })
  }
  
  if (length(.data@compound_list) == 0) {
    .data@compound_list = vector(mode = "list", length = length(.data@pathway_id)) %>%
      purrr::map(function(x) {
        x = data.frame()
        x
      })
  }
  
  if (length(.data@protein_list) == 0) {
    .data@protein_list = vector(mode = "list", length = length(.data@pathway_id)) %>%
      purrr::map(function(x) {
        x = data.frame()
        x
      })
  }
  
  describtion = .data@describtion
  describtion =
    describtion %>%
    purrr::map(function(x) {
      x = paste(x, collapse = ";")
      if (is.null(x)) {
        return("")
      } else{
        x
      }
    }) %>%
    unlist
  
  pathway_class = .data@pathway_class
  pathway_class =
    pathway_class %>%
    purrr::map(function(x) {
      x = paste(x, collapse = ";")
      if (is.null(x)) {
        return("")
      } else{
        x
      }
    }) %>%
    unlist
  
  gene_list = .data@gene_list
  gene_number =
    gene_list %>%
    lapply(nrow) %>%
    unlist()
  
  compound_list = .data@compound_list
  compound_number =
    compound_list %>%
    lapply(nrow) %>%
    unlist()
  
  protein_list = .data@protein_list
  protein_number =
    protein_list %>%
    lapply(nrow) %>%
    unlist()
  
  if (is.null(protein_number)) {
    protein_number = rep(0, length(.data@pathway_id))
  }
  
  temp_data =
    data.frame(
      pathway_id = .data@pathway_id,
      pathway_name = .data@pathway_name,
      describtion,
      pathway_class,
      gene_number,
      compound_number,
      protein_number
    )
  
  
  temp_data =
    dplyr::filter(temp_data, !!!dots, .preserve = .preserve)
  
  idx =
    match(temp_data$pathway_id, .data@pathway_id)
  
  .data@pathway_id = .data@pathway_id[idx]
  .data@pathway_name = .data@pathway_name[idx]
  .data@describtion = .data@describtion[idx]
  .data@pathway_class = .data@pathway_class[idx]
  .data@gene_list = .data@gene_list[idx]
  .data@compound_list = .data@compound_list[idx]
  .data@protein_list = .data@protein_list[idx]
  .data@reference_list = .data@reference_list[idx]
  .data@related_disease = .data@related_disease[idx]
  .data@related_module = .data@related_module[idx]
  
  return(.data)
}

#' @importFrom dplyr filter
#' @export
dplyr::filter
