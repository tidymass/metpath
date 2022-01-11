#' An S4 class to represent pathways
#' @slot database_info database_info
#' @slot pathway_id pathway_id
#' @slot pathway_name pathway_name
#' @slot describtion describtion
#' @slot pathway_class pathway_class
#' @slot gene_list gene_list
#' @slot compound_list compound_list
#' @slot protein_list protein_list
#' @slot reference_list reference_list
#' @slot related_disease related_disease
#' @slot related_module related_module
#' @exportClass pathway_database

setClass(
  Class = "pathway_database",
  representation(
    database_info = "list",
    pathway_id = "vector",
    pathway_name = "vector",
    describtion = "list",
    pathway_class = "list",
    gene_list = "list",
    compound_list = "list",
    protein_list = "list",
    reference_list = "list",
    related_disease = "list",
    related_module = "list"
  ),
  prototype = list(
    database_info = list(),
    pathway_id = list(),
    pathway_name = list(),
    describtion = list(),
    pathway_class = list(),
    gene_list = list(),
    compound_list = list(),
    protein_list = list(),
    reference_list = list(),
    related_disease = list(),
    related_module = list()
  )
)


setMethod(
  f = "show",
  signature = "pathway_database",
  definition = function(object) {
    version <- try(object@database_info$version, silent = TRUE)
    source <- try(object@database_info$source, silent = TRUE)
    if (class(version) != "try-error") {
      cat(crayon::green("---------Pathway source&version---------\n"))
      cat(crayon::green(source, "&", version, "\n"))
    }
    cat(crayon::green("-----------Pathway information------------\n"))
    cat(crayon::green(length(object@pathway_id), "pathways", "\n"))
    cat(
      crayon::green(
        object@gene_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        "pathways have genes",
        "\n"
      )
    )
    
    cat(
      crayon::green(
        object@protein_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        "pathways have proteins",
        "\n"
      )
    )
    
    cat(
      crayon::green(
        object@compound_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        "pathways have compounds",
        "\n"
      )
    )
    
    cat(crayon::green("Pathway class (top 10):",
                      paste(unique(head(
                        unlist(object@pathway_class), 10
                      )), collapse = ";"),
                      "\n"))
  }
)


