#' @title Pathway Database Class
#' @description Defines the `pathway_database` class, which stores information about biological pathways, including genes, compounds, proteins, references, and related diseases.
#' 
#' @author Xiaotao Shen (\email{shenxt1990@outlook.com})
#'
#' @slot database_info A list containing metadata about the pathway database, such as source and version.
#' @slot pathway_id A vector of pathway IDs.
#' @slot pathway_name A vector of pathway names.
#' @slot describtion A list containing descriptions of the pathways.
#' @slot pathway_class A list categorizing pathways into different functional classes.
#' @slot gene_list A list of genes associated with each pathway.
#' @slot compound_list A list of compounds (metabolites) associated with each pathway.
#' @slot protein_list A list of proteins associated with each pathway.
#' @slot reference_list A list of references (e.g., publications) linked to each pathway.
#' @slot related_disease A list of diseases related to each pathway.
#' @slot related_module A list of pathway modules that are related or interconnected.
#' 
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
    if (!is(version, "try-error")) {
      message(crayon::green("---------Pathway source&version---------"))
      message(crayon::green(source, " & ", version))
    }
    message(crayon::green("-----------Pathway information------------"))
    message(crayon::green(length(object@pathway_id), " pathways"))
    message(
      crayon::green(
        object@gene_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        " pathways have genes"
      )
    )
    
    message(
      crayon::green(
        object@protein_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        " pathways have proteins"
      )
    )
    
    message(
      crayon::green(
        object@compound_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        " pathways have compounds"
      )
    )
    
    message(crayon::green("Pathway class (top 10):",
                          paste(unique(head(
                            unlist(object@pathway_class), 10
                          )), collapse = ";")))
  }
)
