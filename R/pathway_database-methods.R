#' @title filter_pathway
#' @description filter pathways according to pathway class
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param object pathway_database object.
#' @param class class pathway class you want to remain.
#' @param remain_idx Which pathways you want to remain (remain_idx).
#' @export

filter_pathway =
  function(object,
           class,
           remain_idx) {
    if (base::class(object) != "pathway_database"){
      stop(crayon::red('Only for pathway_database object.\n'))
    }
    
    if (missing(class) & missing(remain_idx)) {
      stop("Provide class or remain_idx.\n")
    }
    
    if (!missing(class) & !missing(remain_idx)) {
      message(crayon::yellow("Only use remain_idx."))
      remain_idx = remain_idx
    }
    
    if (!missing(class) & missing(remain_idx)) {
      pathway_class = object@pathway_class %>%
        unlist() %>%
        unique()
      class = class[class %in% pathway_class]
      if (length(class) == 0) {
        stop("All class are not in pathway object.\n")
      }
      
      remain_idx =
        object@pathway_class %>%
        purrr::map(function(x) {
          x %in% class
        }) %>%
        unlist() %>%
        which()
    }
    
    if (missing(class) & !missing(remain_idx)) {
      remain_idx = remain_idx
    }
    
    object@pathway_id =
      object@pathway_id[remain_idx]
    
    object@pathway_name =
      object@pathway_name[remain_idx]
    
    object@describtion =
      object@describtion[remain_idx]
    
    object@pathway_class =
      object@pathway_class[remain_idx]
    
    if (length(object@gene_list) > 0) {
      object@gene_list =
        object@gene_list[remain_idx]
    }
    
    if (length(object@compound_list) > 0) {
      object@compound_list =
        object@compound_list[remain_idx]
    }
    
    if (length(object@protein_list) > 0) {
      object@protein_list =
        object@protein_list[remain_idx]
    }
    
    if (length(object@reference_list) > 0) {
      object@reference_list =
        object@reference_list[remain_idx]
    }
    
    if (length(object@related_disease) > 0) {
      object@related_disease =
        object@related_disease[remain_idx]
    }
    
    if (length(object@related_module) > 0) {
      object@related_module =
        object@related_module[remain_idx]
    }
    return(object)
  }


#' @title get_pathway_class
#' @description Extract the class of pathways.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param object pathway_class object.
#' @export

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