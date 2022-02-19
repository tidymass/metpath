#' @title metpath_logo
#' @description Get the detailed of metpath package.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @importFrom KEGGREST keggList keggGet
#' @importFrom pbapply pblapply
#' @importFrom future plan multisession
#' @importFrom furrr future_map
#' @importFrom metid construct_database
#' @importFrom stringr str_replace str_split str_replace_all str_trim
#' @importFrom crayon yellow red num_colors blue col_align col_nchar
#' @importFrom dplyr filter mutate select everything case_when
#' @importFrom purrr map
#' @importFrom openxlsx write.xlsx
#' @importFrom rstudioapi isAvailable hasFun getThemeInfo
#' @importFrom utils packageDescription head
#' @importFrom rlang quos !!!
#' @import ggplot2
#' @import ggraph
#' @importFrom tidygraph tbl_graph
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods .hasSlot new
#' @importFrom stats p.adjust fisher.test phyper
#' @importFrom utils data str
#' @importFrom magrittr %>%
#' @return logo
#' @export

metpath_logo <- function() {
  cat(crayon::green("Thank you for using metpath!\n"))
  message(crayon::green("Version", metpath_version, "(", update_date, ')\n'))
  cat(crayon::green("More information: google tidymass metpath.\n"))
  cat(crayon::green(
    c(
      "                 _   _____      _   _     ",
      "                | | |  __ \\    | | | |    ",
      "  _ __ ___   ___| |_| |__) |_ _| |_| |__  ",
      " | '_ ` _ \\ / _ \\ __|  ___/ _` | __| '_ \\ ",
      " | | | | | |  __/ |_| |  | (_| | |_| | | |",
      " |_| |_| |_|\\___|\\__|_|   \\__,_|\\__|_| |_|",
      "                                          ",
      "                                     "
    )
    
  ), sep = "\n")
}


metpath_version = "0.99.2"
update_date = as.character(Sys.time())
