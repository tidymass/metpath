#' @title get_hmdb_pathway
#' @description Get compound from HMDB (SMPDB)
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param threads threads
#' @export

get_hmdb_compound <- function(threads = 3) {
  data("hmdb_compound_database", envir = environment())
  message(
    crayon::yellow(
      "This database is downloaded in",
      hmdb_compound_database@database.info$Version
    )
  )
  cat("\n")
  return(hmdb_compound_database)
}


#' @title get_hmdb_pathway
#' @description Get pathways from HMDB (SMPDB)
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param threads threads
#' @export
# load("data/hmdb_pathway.rda")
# load("data/primary_pathway")
#
# idx = match(primary_pathway, hmdb_pathway@pathway_id)
# idx = idx[!is.na(idx)]
# pathway_class =
#   hmdb_pathway@pathway_class
#
# for(x in idx){
#   pathway_class[[x]] = paste(pathway_class[[x]], 
#   "primary_pathway", sep = ";")
# }
#
#
#
# pathway_class[idx]
#
# hmdb_pathway@pathway_class = pathway_class
#
# save(hmdb_pathway, file = "data/hmdb_pathway.rda", )

get_hmdb_pathway <- function(threads = 3) {
  data("hmdb_pathway", envir = environment())
  message(
    crayon::yellow(
      "This database is downloaded in",
      hmdb_pathway@database_info$version
    )
  )
  cat("\n")
  return(hmdb_pathway)
}
