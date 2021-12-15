
#' @title get_pathway_class
#' @description Extract the class of pathways.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
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

#' @method [ pathway_database
#' @param x x
#' @param i i
#' @param j j
#' @param drop drop
#' @param .. ..
#' @export
#' @rdname pathway_database-class
#' @return pathway_database
`[.pathway_database` <-
  function(x, i, ...) {
    if (missing(i)) {
      return(x)
    }
    
    if (!missing(i)) {
      if(sum(duplicated(i)) > 0){
        stop("No duplicated i allowed.\n")
      }
      if (is.character(i)) {
        i1 <- match(i, x@pathway_id)
        i2 <- match(i, x@pathway_name)
        if(any(is.na(i1))){
          i = i2
        }else{
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