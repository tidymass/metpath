#' @title Extract Pathway Class
#' @description Extracts the classification of pathways from a `pathway_database` object.
#' @author Xiaotao Shen (\email{shenxt1990@outlook.com})
#' @param object A `pathway_database` class object.
#' @return A data frame with pathway classes and their counts.
#' @export
get_pathway_class <- function(object) {
  object@pathway_class %>%
    unlist() %>%
    data.frame(class = .) %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::ungroup()
}

#' @title [ method
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @method [ pathway_database
#' @param x x
#' @param i i
#' @param ... other parameters
#' @export
#' @rdname extract-pathway_database
#' @return pathway_database
`[.pathway_database` <-
  function(x, i, ...) {
    if (missing(i)) {
      return(x)
    }
    
    if (!missing(i)) {
      if (sum(duplicated(i)) > 0) {
        stop("No duplicated i allowed.\n")
      }
      if (is.character(i)) {
        i1 <- match(i, x@pathway_id)
        i2 <- match(i, x@pathway_name)
        if (any(is.na(i1))) {
          i = i2
        } else{
          i = i1
        }
      }
    } else{
      i = 1:length(x@pathway_id)
    }
    
    if (sum(is.na(i)) > 0) {
      stop("Some i are not in the pathway. Please check it.\n")
    }
    
    if (any(!i %in% 1:length(x@pathway_id))) {
      stop("Some variable index (i) are not in the object. Please check.")
    }
    
    x@pathway_id = x@pathway_id[i]
    x@pathway_name = x@pathway_name[i]
    x@describtion = x@describtion[i]
    x@pathway_class = x@pathway_class[i]
    x@gene_list = x@gene_list[i]
    x@compound_list = x@compound_list[i]
    x@protein_list = x@protein_list[i]
    x@reference_list = x@reference_list[i]
    x@related_disease = x@related_disease[i]
    x@related_module = x@related_module[i]
    
    return(x)
  }


#' @title Get Pathway Names
#' @description Returns pathway IDs from a `pathway_database` object.
#' @author Xiaotao Shen (\email{shenxt1990@outlook.com})
#' @param x A `pathway_database` class object.
#' @return A vector of pathway IDs.
#' @export
#' @rdname extract-pathway_database
names.pathway_database <-
  function(x) {
    x@pathway_id
  }


#' @title Get Pathway Database Length
#' @description Returns the number of pathways in a `pathway_database` object.
#' @author Xiaotao Shen (\email{shenxt1990@outlook.com})
#' @param x A `pathway_database` class object.
#' @return An integer indicating the number of pathways.
#' @export
#' @rdname extract-pathway_database
length.pathway_database <- function(x) {
  length(x@pathway_id)
}



#' @title Extract Database Information
#' @description Extracts database metadata from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A list containing database metadata.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "database_info",
  signature = "pathway_database",
  definition <- function(object) {
    object@database_info
  }
)



#' @title Extract Pathway IDs
#' @description Extracts pathway IDs from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A character vector of pathway IDs.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "pathway_id",
  signature = "pathway_database",
  definition <- function(object) {
    object@pathway_id
  }
)


#' @title Extract Pathway Names
#' @description Extracts pathway names from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A character vector of pathway names.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "pathway_name",
  signature = "pathway_database",
  definition <- function(object) {
    pathway_name = object@pathway_name
    names(pathway_name) = object@pathway_id
    pathway_name
  }
)

#' @title Extract Pathway Descriptions
#' @description Extracts pathway descriptions from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A character vector of pathway descriptions.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "describtion",
  signature = "pathway_database",
  definition <- function(object) {
    describtion = object@describtion
    names(describtion) = object@pathway_id
    describtion
  }
)


#' @title Extract Pathway Classes
#' @description Extracts pathway classification from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A character vector of pathway classes.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "pathway_class",
  signature = "pathway_database",
  definition <- function(object) {
    pathway_class = object@pathway_class
    names(pathway_class) = object@pathway_id
    pathway_class
  }
)


#' @title Extract Gene Lists
#' @description Extracts lists of genes involved in pathways from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A list of gene sets for each pathway.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "gene_list",
  signature = "pathway_database",
  definition <- function(object) {
    gene_list = object@gene_list
    names(gene_list) = object@pathway_id
    gene_list
  }
)


#' @title Extract Compound Lists
#' @description Extracts lists of compounds involved in pathways from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A list of compound sets for each pathway.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "compound_list",
  signature = "pathway_database",
  definition <- function(object) {
    compound_list = object@compound_list
    names(compound_list) = object@pathway_id
    compound_list
  }
)




#' @title Extract Reference Lists
#' @description Extracts literature references for pathways in a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A list of references for each pathway.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "reference_list",
  signature = "pathway_database",
  definition <- function(object) {
    reference_list = object@reference_list
    names(reference_list) = object@pathway_id
    reference_list
  }
)



#' @title Extract Related Diseases
#' @description Extracts disease associations from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A list of related diseases for each pathway.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "related_disease",
  signature = "pathway_database",
  definition <- function(object) {
    related_disease = object@related_disease
    names(related_disease) = object@pathway_id
    related_disease
  }
)


#' @title Extract Related Modules
#' @description Extracts related modules from a `pathway_database` object.
#' @param object A `pathway_database` class object.
#' @return A list of related modules for each pathway.
#' @export
#' @rdname extract-pathway_database
setMethod(
  f = "related_module",
  signature = "pathway_database",
  definition <- function(object) {
    related_module = object@related_module
    names(related_module) = object@pathway_id
    related_module
  }
)
