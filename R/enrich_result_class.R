#'An S4 class to represent pathways
#' @docType class
#' @slot pathway_database pathway_database
#' @slot pathway_version pathway_version
#' @slot result result
#' @slot parameter tidymass
#' @exportClass enrich_result

setClass(
  Class = "enrich_result",
  representation(
    pathway_database = "character",
    pathway_version = "character",
    result = "data.frame",
    parameter = "tidymass_parameter"
  )
)

setMethod(
  f = "show",
  signature = "enrich_result",
  definition = function(object) {
    pathway_version <- try(object@pathway_version, silent = TRUE)
    pathway_database <- try(object@pathway_database, silent = TRUE)
    if (class(version) != "try-error") {
      cat(crayon::green("---------Pathway database&version---------\n"))
      cat(crayon::green(pathway_database, "&", pathway_version, "\n"))
    }
    cat(crayon::green("-----------Enrichment result------------\n"))
    cat(crayon::green(nrow(object@result), "pathways are enriched", "\n"))
    cat(crayon::green(
      nrow(object@result %>%
             dplyr::filter(p_value < 0.05)),
      "pathways p-values < 0.05",
      "\n"
    ))
    
    if (nrow(object@result) > 0) {
      if (nrow(object@result) > 5) {
        pathway_name = object@result$pathway_name[1:5]
        cat(crayon::green(
          paste(pathway_name, collapse = "\n"),
          "... (only top 5 shows)\n"
        ))
      } else{
        pathway_name = object@result$pathway_name
        cat(crayon::green(paste(pathway_name, collapse = ";")))
      }
    }
    
    cat(crayon::green("-----------Parameters------------\n"))
    if (.hasSlot(object = object, name = "parameter")) {
      data.frame(
        "Package" = object@parameter@pacakge_name,
        "Function used" = object@parameter@function_name,
        "Time" = object@parameter@time
      ) %>%
        print()
    } else{
      cat(crayon::yellow("-----------No parameters------------\n"))
    }
  }
)
