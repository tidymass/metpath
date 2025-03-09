#' @title enrich_kegg
#' @description Pathway enrichment for SMPDB database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param query_id The vector of query IDs.
#' @param query_type "compound" or "gene"
#' @param id_type HMDB
#' @param pathway_database KEGG or other metabolomics pathway database.
#' @param p_cutoff p_cutoff
#' @param p_adjust_method p_adjust_method
#' @param method Hypergeometric or fisher test.
#' @param threads threads
#' @return  The MSE analysis result.
#' @export

enrich_kegg <-
  function(query_id,
           query_type = c("compound", "gene"),
           id_type = c("KEGG"),
           pathway_database,
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
    query_type = match.arg(query_type)
    id_type = match.arg(id_type)
    method <- match.arg(method)
    p_adjust_method = match.arg(p_adjust_method)
    
    if (query_type == "compound") {
      if (pathway_database@database_info$source != "KEGG") {
        stop("pathway_database must from KEGG.\n")
      }
      
      remain_idx =
        pathway_database@compound_list %>%
        purrr::map(function(x) {
          nrow(x)
        }) %>%
        unlist() %>%
        `>`(0) %>%
        which()
      
      pathway_database =
        filter_pathway(object = pathway_database, remain_idx = remain_idx)
      
      for (i in 1:length(pathway_database@pathway_id)) {
        pathway_database@compound_list[[i]] =
          pathway_database@compound_list[[i]] %>%
          dplyr::filter(!is.na(KEGG.ID)) %>%
          dplyr::filter(KEGG.ID != "")
      }
      
      ##remove pathways without compound
      remain_idx =
        pathway_database@compound_list %>%
        purrr::map(function(x) {
          nrow(x)
        }) %>%
        unlist() %>%
        `>`(0) %>%
        which()
      
      pathway_database =
        filter_pathway(object = pathway_database, remain_idx = remain_idx)
      
      if (length(pathway_database@pathway_id) == 0) {
        stop("Pathway database length is zero.\n")
      } else{
        message(crayon::green(length(pathway_database@describtion), "pathways."))
      }
      
      database = pathway_database@compound_list
      names(database) = pathway_database@pathway_id
      
      database =
        database %>%
        purrr::map(function(x) {
          x$KEGG.ID
        })
      
      for (idx in 1:length(pathway_database@describtion)) {
        if (is.null(pathway_database@describtion[[idx]])) {
          pathway_database@describtion[[idx]] = NA
        }
        pathway_database@describtion[[idx]] =
          paste(pathway_database@describtion[[idx]], collapse = "{}")
      }
      
      pathway_info <- data.frame(
        pathway_id = pathway_database@pathway_id,
        pathway_name = pathway_database@pathway_name,
        describtion = unlist(pathway_database@describtion),
        pathway_class = unlist(pathway_database@pathway_class),
        stringsAsFactors = FALSE
      )
    }
    
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
    result =
      data.frame(pathway_info, p_value, p_value_adjust)
    
    result$all_id =
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          paste(database[[x]], collapse = ";")
        }
      ) %>%
      unlist()
    
    result$all_number =
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          length(database[[x]])
        }
      ) %>%
      unlist()
    
    result$mapped_id =
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          paste(query_id[query_id %in% database[[x]]], collapse = ";")
        }
      ) %>%
      unlist()
    
    result$mapped_number =
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          sum(query_id %in% database[[x]])
        }
      ) %>%
      unlist()
    
    result$mapped_percentage =
      furrr::future_map(
        .x = result$pathway_id,
        .f = function(x) {
          sum(query_id %in% database[[x]]) * 100 / length(database[[x]])
        }
      ) %>%
      unlist()
    
    result <-
      result %>%
      dplyr::arrange(p_value_adjust)
    
    return_result =
      new(
        Class = "enrich_result",
        pathway_database = pathway_database@database_info$source,
        pathway_version = metpath_version,
        result = result,
        parameter = new(
          Class = "tidymass_parameter",
          pacakge_name = "metpath",
          function_name = "enrich_kegg()",
          parameter = list(
            query_id = query_id,
            query_type = query_type,
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
    
  }

