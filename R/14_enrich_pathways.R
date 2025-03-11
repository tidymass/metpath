# ##wikipathways
# library(r4projects)
# setwd(get_project_wd())
# load("data/wikipathway_hsa_pathway.rda")
# load("data/query_id_hmdb.rda")
# load("data/query_id_kegg.rda")
# query_id <-
#   query_id_hmdb
# p_cutoff = 0.05
# p_adjust_method <- "fdr"
# method <- "hypergeometric"
# threads = 3
# pathway_database <- wikipathway_hsa_pathway
# only_primary_pathway = FALSE
#
# result1 <-
#   enrich_pathways(
#     query_id = query_id_hmdb,
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = pathway_database,
#     only_primary_pathway = only_primary_pathway,
#     p_cutoff = p_cutoff,
#     p_adjust_method = p_adjust_method,
#     method = method,
#     threads = threads
#   )
#
# result2 <-
#   enrich_pathways(
#     query_id = query_id_kegg,
#     query_type = "compound",
#     id_type = "KEGG",
#     pathway_database = pathway_database,
#     only_primary_pathway = only_primary_pathway,
#     p_cutoff = p_cutoff,
#     p_adjust_method = p_adjust_method,
#     method = method,
#     threads = threads
#   )
#
# enrich_bar_plot(object = result1, top = 10)
# enrich_bar_plot(object = result2, top = 10)



#' @title Pathway Enrichment Analysis
#' @description This function performs pathway enrichment analysis using
#' metabolic pathway databases. It supports compound and metabolite enrichment
#' based on HMDB or KEGG identifiers.
#'
#' @param query_id A vector of compound or metabolite identifiers.
#' @param query_type Character. Type of query, either `"compound"` or `"protein"`. Default: `"compound"`.
#' @param id_type Character. Type of identifier, either `"HMDB"` or `"KEGG"`. Default: `"HMDB"`.
#' @param pathway_database The pathway database object containing pathway information.
#' @param only_primary_pathway Logical. Whether to filter only primary pathways. Default: `FALSE`.
#' @param p_cutoff Numeric. The significance threshold for pathway enrichment. Default: `0.05`.
#' @param p_adjust_method Character. P-value adjustment method. Options include `"holm"`, `"hochberg"`,
#' `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`, `"none"`. Default: `"holm"`.
#' @param method Character. Statistical method used for enrichment analysis. Options: `"hypergeometric"`
#' or `"fisher"`. Default: `"hypergeometric"`.
#' @param threads Integer. Number of threads to use for computation. Default: `3`.
#'
#' @return A data frame containing the pathway enrichment results, including pathway ID,
#' pathway name, description, p-values, adjusted p-values, and mapping statistics.
#'
#' @importFrom dplyr filter mutate pull arrange
#' @importFrom purrr map map_dfr
#' @importFrom stringr str_detect str_split
#' @importFrom furrr future_map
#' @importFrom plyr dlply .
#' @importFrom stats p.adjust fisher.test phyper
#'
#' @examples
#' \dontrun{
#'   data("sample_pathway_database") # Load example pathway database
#'   query_metabolites <- c("HMDB0000122", "HMDB0000532", "HMDB0000889")
#'   enrichment_results <- enrich_pathways(
#'     query_id = query_metabolites,
#'     query_type = "compound",
#'     id_type = "HMDB",
#'     pathway_database = sample_pathway_database,
#'     p_cutoff = 0.05,
#'     method = "fisher"
#'   )
#'   print(enrichment_results)
#' }
#'
#' @export

enrich_pathways <-
  function(query_id,
           query_type = c("compound", "protein"),
           id_type = c("HMDB", "KEGG"),
           pathway_database,
           only_primary_pathway = FALSE,
           p_cutoff = 0.05,
           p_adjust_method = c("holm",
                               "hochberg",
                               "hommel",
                               "bonferroni",
                               "BH",
                               "BY",
                               "fdr",
                               "none"),
           method = c("hypergeometric", "fisher"),
           threads = 3) {
    query_type <- match.arg(query_type)
    id_type <- match.arg(id_type)
    method <- match.arg(method)
    p_adjust_method <- match.arg(p_adjust_method)
    
    if (missing(pathway_database)) {
      stop("Please provide pathway database.\n")
    }
    # browser()
    if (query_type == "compound") {
      result <-
        enrich_metabolic_pathway(
          query_id = query_id,
          id_type = id_type,
          pathway_database = pathway_database,
          only_primary_pathway = only_primary_pathway,
          p_cutoff = p_cutoff,
          p_adjust_method = p_adjust_method,
          method = method,
          threads = threads
        )
      return(result)
    }
  }



#' @title Metabolic Pathway Enrichment Analysis
#' @description Performs enrichment analysis for metabolic pathways using HMDB or KEGG identifiers.
#' The function applies statistical tests to assess the overrepresentation of queried
#' compounds or proteins in metabolic pathways.
#'
#' @param query_id A vector of compound or protein identifiers.
#' @param id_type Character. Type of identifier, either `"HMDB"` or `"KEGG"`. Default: `"HMDB"`.
#' @param pathway_database The pathway database object containing metabolic pathway information.
#' @param only_primary_pathway Logical. Whether to retain only primary pathways. Default: `FALSE`.
#' @param p_cutoff Numeric. P-value cutoff for statistical significance. Default: `0.05`.
#' @param p_adjust_method Character. Method to adjust p-values. Options include `"holm"`, `"hochberg"`,
#' `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`, `"none"`. Default: `"holm"`.
#' @param method Character. Statistical method for enrichment analysis. Options: `"hypergeometric"` or `"fisher"`.
#' Default: `"hypergeometric"`.
#' @param threads Integer. Number of computational threads. Default: `3`.
#'
#' @return A list containing pathway enrichment results, including pathway ID, pathway name,
#' pathway description, and enrichment statistics.
#'
#' @importFrom dplyr filter pull mutate arrange
#' @importFrom purrr map
#' @importFrom stringr str_detect str_split
#' @importFrom furrr future_map
#' @importFrom stats fisher.test p.adjust phyper
#' @importFrom plyr dlply .
#'
#' @examples
#' \dontrun{
#'   data("sample_pathway_database")
#'   query_compounds <- c("HMDB0000122", "HMDB0000532")
#'   metabolic_enrichment <- enrich_metabolic_pathway(
#'     query_id = query_compounds,
#'     id_type = "HMDB",
#'     pathway_database = sample_pathway_database,
#'     p_cutoff = 0.05,
#'     method = "fisher"
#'   )
#'   print(metabolic_enrichment)
#' }
#'
#' @export

enrich_metabolic_pathway <-
  function(query_id,
           id_type = c("HMDB", "KEGG"),
           pathway_database,
           only_primary_pathway = FALSE,
           p_cutoff = 0.05,
           p_adjust_method = c("holm",
                               "hochberg",
                               "hommel",
                               "bonferroni",
                               "BH",
                               "BY",
                               "fdr",
                               "none"),
           method = c("hypergeometric", "fisher"),
           threads = 3) {
    # browser()
    id_type <- match.arg(id_type)
    method <- match.arg(method)
    p_adjust_method <- match.arg(p_adjust_method)
    
    remain_idx <-
      pathway_database@compound_list %>%
      purrr::map(function(x) {
        nrow(x)
      }) %>%
      unlist() %>%
      `>`(0) %>%
      which()
    
    pathway_database <-
      filter_pathway(object = pathway_database, remain_idx = remain_idx)
    
    has_metabolite_id <-
      lapply(pathway_database@compound_list, function(x) {
        any(colnames(x) == paste(id_type, ".ID", sep = ""))
      }) %>%
      unlist() %>%
      all()
    
    if (!has_metabolite_id) {
      stop("The compound list should contains ",
           paste0(id_type, ".ID"),
           ".\n")
    }
    
    remain_idx <-
      pathway_database@compound_list %>%
      purrr::map(function(x) {
        nrow(x)
      }) %>%
      unlist() %>%
      `>`(0) %>%
      which()
    
    pathway_database <-
      filter_pathway(object = pathway_database, remain_idx = remain_idx)
    
    ##only remain primary class pathway
    if (only_primary_pathway) {
      remain_idx =
        pathway_database@pathway_class %>%
        purrr::map(function(x) {
          stringr::str_detect(x, "primary_pathway")
        }) %>%
        unlist() %>%
        which()
      if (length(remain_idx) == 0) {
        stop(
          "No pathways left if you set 'only_primary_pathway' as TRUE.
             Try to set 'only_primary_pathway' as FALSE.\n"
        )
      }
      
      pathway_database =
        filter_pathway(object = pathway_database, remain_idx = remain_idx)
    }
    
    for (i in 1:length(pathway_database@pathway_id)) {
      x <-
        pathway_database@compound_list[[i]]
      colnames(x)[colnames(x) == paste0(id_type, ".ID")] <- "temp.ID"
      pathway_database@compound_list[[i]] <-
        x %>%
        dplyr::filter(!is.na(temp.ID)) %>%
        dplyr::filter(temp.ID != "")
    }
    
    ##remove pathways without compound
    remain_idx <-
      pathway_database@compound_list %>%
      purrr::map(function(x) {
        nrow(x)
      }) %>%
      unlist() %>%
      `>`(0) %>%
      which()
    
    pathway_database <-
      filter_pathway(object = pathway_database, remain_idx = remain_idx)
    
    if (length(pathway_database@pathway_id) == 0) {
      stop("Pathway database length is zero.\n")
    } else{
      message(crayon::green(length(pathway_database@describtion), "pathways."))
    }
    
    database <- pathway_database@compound_list
    names(database) <- pathway_database@pathway_id
    
    database <-
      database %>%
      purrr::map(function(x) {
        x %>%
          dplyr::pull(temp.ID) %>%
          unique() %>%
          unname() %>%
          as.character() %>%
          stringr::str_split("\\{\\}") %>%
          unlist() %>%
          unique()
      })
    
    describtion <-
      pathway_database@describtion %>%
      lapply(function(x) {
        if (is.null(x)) {
          return(NA)
        }
        if (length(x) == 0) {
          return(NA)
        }
        return(paste0(x, collapse = ";"))
      }) %>%
      unlist()
    
    pathway_class <-
      pathway_database@pathway_class %>%
      lapply(function(x) {
        if (is.null(x)) {
          return(NA)
        }
        if (length(x) == 0) {
          return(NA)
        }
        return(x)
      }) %>%
      unlist()
    
    if (is.null(describtion)) {
      describtion = rep(NA, length(pathway_database@pathway_id))
    }
    
    if (is.null(pathway_class)) {
      pathway_class = rep(NA, length(pathway_database@pathway_id))
    }
    
    pathway_info = data.frame(
      pathway_id = pathway_database@pathway_id,
      pathway_name = pathway_database@pathway_name,
      describtion = describtion,
      pathway_class = pathway_class,
      stringsAsFactors = FALSE
    )
    
    
    all_id <-
      database %>%
      unlist() %>%
      unique() %>%
      unname() %>%
      as.character()
    
    ##remove the ID which is not in the all_id
    query_id <- query_id[which(is.element(query_id, all_id))]
    query_id <- unique(query_id[!is.na(query_id)])
    
    if (length(query_id) == 0) {
      return(NULL)
    }
    
    sig_id <- as.character(query_id)
    
    num_all <- length(all_id)
    num_sig <- length(sig_id)
    
    # --------------------------------------------------------------
    #Generating label matrix for detected metabolites
    # --------------------------------------------------------------
    
    all_matrix <-
      set_label(query_id = all_id,
                database = database,
                threads = threads)
    
    # delete metabolite set
    all_matrix2 <-
      all_matrix[, colSums(all_matrix) != 0, drop = FALSE]
    
    # error handling
    if (ncol(all_matrix2) < 2) {
      # stop function
      return(NULL)
    }
    
    # --------------------------------------------------------------
    #Generating label matrix for significant metabolites
    # --------------------------------------------------------------
    
    sig_matrix <-
      set_label(query_id = sig_id,
                database = database,
                threads = threads)
    
    sig_matrix <-
      sig_matrix[, colSums(all_matrix) != 0, drop = FALSE]
    
    # -------------------------------
    #Calculating  ORA
    # -------------------------------
    #for each pathway
    p_value =
      furrr::future_map(
        .x = 1:ncol(all_matrix2),
        .f = function(i) {
          # ------------------------------------
          #Generating 2~2 table
          # -------------------------------------
          a1 <-
            sum(sig_matrix[, i])# significant and including pathway
          a2 <-
            sum(all_matrix2[, i]) - sum(sig_matrix[, i])
          # not significant and including pathway
          a3 <-
            length(sig_id) - a1
          # significant and not including pathway
          a4 <-
            (length(all_id) - length(sig_id)) - a2
          # not significant and not including pathway
          
          tab <- t(matrix(c(a1, a2, a3, a4), 2))
          if (method == "hypergeometric") {
            phyper(
              q = a1 - 1,
              m = sum(all_matrix2[, i]),
              n = num_all - sum(all_matrix2[, i]),
              k = num_sig,
              lower.tail = FALSE
            )
          } else{
            # ----------------------------------
            # Fisher's exact test
            # ----------------------------------
            check <- tryCatch({
              resfish <- fisher.test(tab, alternative = "greater")
            }, error = function(e) {
              NA
            })
            
            if (!is(object = check, class2 = "htest")) {
              return(1)
            } else{
              resfish <- fisher.test(tab, alternative = "greater")
              return(resfish$p.value)
            }
          }
        },
        .progress = TRUE
      ) %>%
      unlist()
    
    # -----------------------
    #q-value
    # -----------------------
    p_value_adjust <- p.adjust(p_value, method = p_adjust_method)
    
    # ----------------------
    #Result
    # ----------------------
    result <-
      data.frame(pathway_info, p_value, p_value_adjust)
    
    result$all_id <-
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          paste(database[[x]], collapse = ";")
        }
      ) %>%
      unlist()
    
    result$all_number <-
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          length(database[[x]])
        }
      ) %>%
      unlist()
    
    result$mapped_id <-
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          paste(query_id[query_id %in% database[[x]]], collapse = ";")
        }
      ) %>%
      unlist()
    
    result$mapped_number <-
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          sum(query_id %in% database[[x]])
        }
      ) %>%
      unlist()
    
    result$mapped_percentage <-
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          sum(query_id %in% database[[x]]) * 100 / length(database[[x]])
        }
      ) %>%
      unlist()
    
    # result =
    #   result %>%
    #   dplyr::arrange(p_value) %>%
    #   dplyr::filter(p_value <= p_cutoff)
    
    ##remove the duplicated pathways
    result <-
      result %>%
      plyr::dlply(.variables = plyr::.(pathway_name)) %>%
      purrr::map(function(x) {
        if (nrow(x) == 1) {
          return(x)
        } else{
          x =
            x %>%
            dplyr::filter(p_value_adjust == min(p_value_adjust)) %>%
            dplyr::filter(mapped_number == max(mapped_number)) %>%
            dplyr::filter(all_number == max(all_number))
          x[1, , drop = FALSE]
        }
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    result <-
      result %>%
      dplyr::arrange(p_value_adjust)
    
    return_result <-
      new(
        Class = "enrich_result",
        pathway_database = pathway_database@database_info$source,
        pathway_version = metpath_version,
        result = result,
        parameter = new(
          Class = "tidymass_parameter",
          pacakge_name = "metpath",
          function_name = "enrich_hmdb()",
          parameter = list(
            query_id = query_id,
            query_type = "compound",
            id_type = id_type,
            pathway_database = paste(
              pathway_database@database_info$source,
              pathway_database@database_info$version,
              sep = ","
            ),
            p_cutoff = p_cutoff,
            p_adjust_method = p_adjust_method,
            method = method,
            threads = threads
          ),
          time = Sys.time()
        )
      )
    
    return(return_result)
    
  }
