#' @title get_pathway_class
#' @description Extract the class of pathways.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object pathway_class object.
#' @export
#' @examples
#' library(metpath)
#' kegg_hsa_pathway = get_kegg_pathway(local = TRUE, organism = "hsa")
#' kegg_hsa_pathway
#'
#' get_pathway_class(kegg_hsa_pathway)

get_pathway_class = function(object) {
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
#' @rdname pathway_database-class
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


#' @title names method
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @method names pathway_database
#' @param x pathway_database class object
#' @export
#' @rdname pathway_database-class
#' @return ID of pathways
names.pathway_database <-
  function(x) {
    x@pathway_id
  }


#' @title length method
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @method length pathway_database
#' @param x x
#' @export
#' @rdname pathway_database-class
#' @return message
length.pathway_database <- function(x) {
  length(x@pathway_id)
}



#' @method database_info pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A list (database_info)
setMethod(
  f = "database_info",
  signature = "pathway_database",
  definition = function(object) {
    object@database_info
  }
)



##pathway_id
#' @method pathway_id pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (pathway_id)
setMethod(
  f = "pathway_id",
  signature = "pathway_database",
  definition = function(object) {
    object@pathway_id
  }
)


##pathway_name
#' @method pathway_name pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (pathway_name)
setMethod(
  f = "pathway_name",
  signature = "pathway_database",
  definition = function(object) {
    pathway_name = object@pathway_name
    names(pathway_name) = object@pathway_id
    pathway_name
  }
)

##describtion
#' @method describtion pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (describtion)
setMethod(
  f = "describtion",
  signature = "pathway_database",
  definition = function(object) {
    describtion = object@describtion
    names(describtion) = object@pathway_id
    describtion
  }
)


##pathway_class
#' @method pathway_class pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (pathway_class)
setMethod(
  f = "pathway_class",
  signature = "pathway_database",
  definition = function(object) {
    pathway_class = object@pathway_class
    names(pathway_class) = object@pathway_id
    pathway_class
  }
)


##gene_list
#' @method gene_list pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (gene_list)
setMethod(
  f = "gene_list",
  signature = "pathway_database",
  definition = function(object) {
    gene_list = object@gene_list
    names(gene_list) = object@pathway_id
    gene_list
  }
)


##compound_list
#' @method compound_list pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (compound_list)
setMethod(
  f = "compound_list",
  signature = "pathway_database",
  definition = function(object) {
    compound_list = object@compound_list
    names(compound_list) = object@pathway_id
    compound_list
  }
)




##reference_list
#' @method reference_list pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (reference_list)
setMethod(
  f = "reference_list",
  signature = "pathway_database",
  definition = function(object) {
    reference_list = object@reference_list
    names(reference_list) = object@pathway_id
    reference_list
  }
)



##related_disease
#' @method related_disease pathway_database
#' @docType methods
#' @rdname extract-pathway_database
#' @export
#' @param object (required) pathway_database class object
#' @return A vector (related_disease)
setMethod(
  f = "related_disease",
  signature = "pathway_database",
  definition = function(object) {
    related_disease = object@related_disease
    names(related_disease) = object@pathway_id
    related_disease
  }
)


##related_module
#' @method related_module pathway_database
#' @docType methods
#' @export
#' @rdname extract-pathway_database
#' @param object (required) pathway_database class object
#' @return A vector (related_module)
setMethod(
  f = "related_module",
  signature = "pathway_database",
  definition = function(object) {
    related_module = object@related_module
    names(related_module) = object@pathway_id
    related_module
  }
)
