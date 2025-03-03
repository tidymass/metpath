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
#' @importClassesFrom massdataset tidymass_parameter
#' @return logo
#' @export

metpath_logo <- function() {
  message("Thank you for using metpath!")
  message("Version ", metpath_version, " (", update_date, ')')
  message("More information: metpath.tidymass.org")
  cat(
    c(
      "                 _   _____      _   _     ",
      "                | | |  __ \\    | | | |    ",
      "  _ __ ___   ___| |_| |__) |_ _| |_| |__  ",
      " | '_ ` _ \\ / _ \\ __|  ___/ _` | __| '_ \\ ",
      " | | | | | |  __/ |_| |  | (_| | |_| | | |",
      " |_| |_| |_|\\___|\\__|_|   \\__,_|\\__|_| |_|",
      "                                          ",
      "                                     "
    ),
    sep = "\n"
  )
}

metpath_version <-
  as.character(utils::packageVersion(pkg = "metpath"))
update_date <-
  as.character(Sys.time())
