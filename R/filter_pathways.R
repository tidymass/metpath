#' @title filter_pathway
#' @description filter pathways according to pathway class
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object pathway_database object.
#' @param class class pathway class you want to remain.
#' @param remain_idx Which pathways you want to remain (remain_idx).
#' @export

filter_pathway =
  function(object,
           class,
           remain_idx) {
    if (base::class(object) != "pathway_database") {
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
