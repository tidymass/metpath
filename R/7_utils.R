

msg <- function(..., startup = FALSE) {
  if (startup) {
    if (!isTRUE(getOption("metpath.quiet"))) {
      packageStartupMessage(text_col(...))
    }
  } else {
    message(text_col(...))
  }
}

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  
  theme <- rstudioapi::getThemeInfo()
  
  if (isTRUE(theme$dark))
    crayon::white(x)
  else
    crayon::black(x)
  
}

#' List all packages in the metpath
#'
#' @param include_self Include metpath in the list?
#' @export
#' @examples
#' metpath_packages()
metpath_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("metpath")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <-
    vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  
  if (include_self) {
    names <- c(names, "metpath")
  }
  
  names
}

invert <- function(x) {
  if (length(x) == 0)
    return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(paste0(...),
                crayon::make_style(grDevices::grey(level), grey = TRUE))
}



#' @title set_label
#' @description set_label
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param query_id query_id
#' @param database database
#' @param threads threads
#' @return  The MSE analysis result.

set_label <-
  function (query_id,
            database,
            threads = parallel::detectCores() - 2) {
    future::plan(strategy = future::multisession, workers = threads)
    return_result =
      furrr::future_map(
        .x = database,
        .f = function(x) {
          temp = match(query_id, x)
          temp[!is.na(temp)] = 1
          temp[is.na(temp)] = 0
          temp
        }
      ) %>%
      do.call(cbind, .) %>%
      as.data.frame()
    
    rownames(return_result) = query_id
    return_result
  }
