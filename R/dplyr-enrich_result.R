#' @method filter enrich_result
#' @docType methods
#' @importFrom rlang quos !!!
#' @importFrom dplyr filter
#' @export
filter.enrich_result <- function(.data, ..., .preserve = FALSE) {
  dots <- quos(...)
  
  .data@result =
    dplyr::filter(.data@result, !!!dots, .preserve = .preserve)
  return(.data)
}

#' @importFrom dplyr filter
#' @export
dplyr::filter


#' @method arrange enrich_result
#' @docType methods
#' @importFrom rlang quos !!!
#' @importFrom dplyr arrange
#' @export
arrange.enrich_result <- function(.data, ...) {
  dots <- rlang::quos(...)
  
  .data@result =
    dplyr::arrange(.data@result, !!!dots)
  return(.data)
}

#' @importFrom dplyr arrange
#' @export
dplyr::arrange
