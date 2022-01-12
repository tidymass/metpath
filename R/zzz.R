.onAttach <- function(libname, pkgname) {
  needed <- core[!is_attached(core)]
  if (length(needed) == 0)
    return()
  
  crayon::num_colors(TRUE)
  metpath_attach()
  
  
  if (!"package:conflicted" %in% search()) {
    x <- metpath_conflicts()
    msg(metpath_conflict_message(x), startup = TRUE)
  }
  
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}
