# library(r4projects)
# setwd(get_project_wd())
# source("R/enrich_pathways.R")
# source("R/7_utils.R")
#
# setwd("demo_data/2_smartd_project")
# load("peak_marker")
# load("urine_metabolomics_data.rda")
# # load("hmdb_compound_ms1.rda")
# load("kegg_compound_ms1.rda")
# # load("metabolic_network.rda")
# # edge_data <- tidygraph::as_tibble(metabolic_network, active = "edges")
# # node_data <- tidygraph::as_tibble(metabolic_network, active = "nodes")
# #
# # node_data <-
# #   node_data %>%
# #   dplyr::filter(!is.na(KEGG_ID))
# #
# # edge_data <-
# #   edge_data %>%
# #   dplyr::filter(!is.na(from_compound_KEGG_ID) &
# #                   !is.na(to_compound_KEGG_ID)) %>%
# #   dplyr::filter(
# #     from_compound_KEGG_ID %in% node_data$KEGG_ID &
# #       to_compound_KEGG_ID %in% node_data$KEGG_ID
# #   )
# #
# # node_data <-
# #   node_data %>%
# #   dplyr::filter(KEGG_ID %in% unique(
# #     c(
# #       edge_data$from_compound_KEGG_ID,
# #       edge_data$to_compound_KEGG_ID
# #     )
# #   ))
# #
# # node_data <-
# #   node_data %>%
# #   dplyr::mutate(name = KEGG_ID) %>%
# #   dplyr::distinct(name, .keep_all = TRUE)
# #
# # edge_data <-
# #   edge_data %>%
# #   dplyr::mutate(from = from_compound_KEGG_ID, to = to_compound_KEGG_ID) %>%
# #   dplyr::filter(from %in% node_data$name &
# #                   to %in% node_data$name) %>%
# #   dplyr::distinct(from, to, .keep_all = TRUE)
# #
# # node_data <-
# #   node_data %>%
# #   dplyr::filter(name %in% unique(c(edge_data$from, edge_data$to)))
# #
# # metabolic_network <-
# #   tidygraph::tbl_graph(nodes = node_data,
# #                        edges = edge_data,
# #                        directed = FALSE)
# #
# # save(metabolic_network, file = "metabolic_network.rda")
# # load("metabolic_network.rda")
# #
# ###remove duplicated edges
# # node_data <- tidygraph::as_tibble(metabolic_network, active = "nodes")
# # edge_data <- tidygraph::as_tibble(metabolic_network, active = "edges")
# # edge_data <-
# #   edge_data %>%
# #   dplyr::mutate(from = from_compound_KEGG_ID, to = to_compound_KEGG_ID) %>%
# #   dplyr::filter(from != to) %>%
# #   dplyr::distinct(from, to, .keep_all = TRUE)
# #
# # edge_data <-
# #   edge_data %>%
# #   dplyr::mutate(rpair_id = apply(cbind(from, to), 1, function(x)
# #     paste(sort(as.character(
# #       x
# #     )), collapse = "_"))) %>%
# #   dplyr::distinct(rpair_id, .keep_all = TRUE)
# #
# #
# # metabolic_network <-
# #   tidygraph::tbl_graph(nodes = node_data,
# #                        edges = edge_data,
# #                        directed = FALSE)
# #
# # save(metabolic_network, file = "metabolic_network.rda")
# load("metabolic_network.rda")
# # load("hmdb_pathway.rda")
# load("kegg_hsa_pathway.rda")
#
# # ###for HMDB, only remain the primary_pathway
# # hmdb_pathway <-
# #   hmdb_pathway %>%
# #   dplyr::filter(stringr::str_detect(pathway_class, "primary_pathway"))
#
# library(tidymass)
# library(tidyverse)
# library(igraph)
# library(ggraph)
#
# up_marker <-
#   peak_marker %>%
#   dplyr::filter(score > 0 & correlation1 > 0.3 & p1_adj < 0.05)  %>%
#   dplyr::select("Gene ID", score, p1_adj)
#
# down_marker <-
#   peak_marker %>%
#   dplyr::filter(score < 0 &
#                   correlation1 < -0.3 & p1_adj < 0.05)  %>%
#   dplyr::select("Gene ID", score, p1_adj)
#
# up_marker <-
#   up_marker %>%
#   dplyr::left_join(urine_metabolomics_data@variable_info,
#                    by = c("Gene ID" = "variable_id"))
#
# down_marker <-
#   down_marker %>%
#   dplyr::left_join(urine_metabolomics_data@variable_info,
#                    by = c("Gene ID" = "variable_id"))
#
# feature_table_marker_up <-
#   up_marker %>%
#   dplyr::select("Gene ID", mz, rt, score, p1_adj) %>%
#   dplyr::rename(variable_id = "Gene ID",
#                 degree = score,
#                 p_value = p1_adj) %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(variable_id, "POS") ~ "positive",
#     stringr::str_detect(variable_id, "NEG") ~ "negative"
#   ))
#
# feature_table_marker_down <-
#   down_marker %>%
#   dplyr::select("Gene ID", mz, rt, score, p1_adj) %>%
#   dplyr::rename(variable_id = "Gene ID",
#                 degree = score,
#                 p_value = p1_adj) %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(variable_id, "POS") ~ "positive",
#     stringr::str_detect(variable_id, "NEG") ~ "negative"
#   ))
#
# feature_table_marker <-
#   rbind(feature_table_marker_up, feature_table_marker_down)
#
# feature_table_all <-
#   urine_metabolomics_data@variable_info %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(variable_id, "POS") ~ "positive",
#     stringr::str_detect(variable_id, "NEG") ~ "negative"
#   )) %>%
#   dplyr::select(variable_id, mz, rt, polarity)
#
# metabolite_database <-
#   kegg_compound_ms1
#
# pathway_database <-
#   kegg_hsa_pathway
#
# column <- "rp"
# adduct.table <- NULL
# adduct.table = NULL
# ms1.match.ppm = 15
# rt.match.tol = 5
# mz.ppm.thr = 400
# threads = 3
# include_hidden_metabolites = FALSE
#
# fpa_result <-
#   perform_fpa(
#     feature_table_marker = feature_table_marker,
#     feature_table_all = feature_table_all,
#     metabolite_database = metabolite_database,
#     column = column,
#     adduct.table = adduct.table,
#     ms1.match.ppm = ms1.match.ppm,
#     rt.match.tol = rt.match.tol,
#     mz.ppm.thr = mz.ppm.thr,
#     threads = threads,
#     include_hidden_metabolites = include_hidden_metabolites,
#     metabolic_network = metabolic_network,
#     pathway_database = pathway_database
#   )
#
#
# #####network visualization
# result$enriched_pathways
#
# plot_metabolic_network_fpa(
#   fpa_result,
#   feature_table_marker = feature_table_marker,
#   include_feature = FALSE,
#   include_hidden_metabolites = FALSE,
#   add_compound_name = TRUE
# )
#
# plot_metabolic_module_fpa(
#   fpa_result,
#   feature_table_marker = feature_table_marker,
#   include_feature = FALSE,
#   include_hidden_metabolites = FALSE,
#   add_compound_name = TRUE,
#   metabolic_module_index = 6
# )
#

#' Perform Feature-based Pathway Analysis (FPA)
#'
#' This function performs feature-based pathway analysis. It contains
#' metabolite annotation based on ms1, network-based metabolite filtering,
#' metabolic module identification.
#'
#' @param feature_table_marker A data frame of marker features.
#' @param feature_table_all A data frame containing all features.
#' @param metabolite_database A metabolite database object containing metabolite spectra.
#' @param column Character vector specifying the chromatography column type. Default is `c("rp", "hilic")`.
#' @param adduct.table Optional data frame containing adduct information. Default is `NULL`.
#' @param ms1.match.ppm Numeric; mass tolerance for MS1 matching in parts per million (ppm). Default is `10`.
#' @param rt.match.tol Numeric; retention time tolerance threshold in seconds. Default is `5`.
#' @param mz.ppm.thr Numeric; M/Z tolerance threshold for filtering in ppm. Default is `400`.
#' @param threads Integer; number of threads to use for parallel processing. Default is `3`.
#' @param include_hidden_metabolites Logical; whether to include hidden metabolites. Default is `FALSE`.
#' @param metabolic_network A metabolic network object representing metabolic reactions.
#' @param pathway_database A pathway database object containing metabolic pathways.
#'
#' @return A list containing:
#'   - `metabolic_module_result`: A data frame with metabolic module information.
#'   - `final_metabolic_network`: A `tbl_graph` object representing the final metabolic network.
#'   - `annotation_table_final`: A data frame with final annotated features.
#'
#' @export
#'

perform_fpa <-
  function(feature_table_marker,
           feature_table_all,
           metabolite_database,
           column = c("rp", "hilic"),
           adduct.table = NULL,
           ms1.match.ppm = 25,
           rt.match.tol = 5,
           ##unit is second
           mz.ppm.thr = 400,
           threads = 3,
           include_hidden_metabolites = FALSE,
           metabolic_network,
           pathway_database) {
    ###check the feature_table_marker and feature_table_all
    check_feature_table_marker(feature_table_marker)
    check_feature_table_all(feature_table_all)
    ###metabolite annotation
    message("Annotating metabolites for feature table...\n")
    
    # Extract metabolic network data: nodes (metabolites) and edges (connections)
    edge_data <-
      tidygraph::as_tibble(metabolic_network, active = "edges")
    
    node_data <-
      tidygraph::as_tibble(metabolic_network, active = "nodes")
    
    # Remove entries with missing KEGG IDs from the metabolite database
    metabolite_database@spectra.info <-
      metabolite_database@spectra.info %>%
      dplyr::filter(!is.na(KEGG.ID)) %>%
      dplyr::filter(KEGG.ID %in% node_data$KEGG_ID)
    
    # Perform metabolite annotation based on ms1
    annotation_table_all <-
      annotate_metabolites_fpa(
        feature_table = feature_table_all,
        metabolite_database = metabolite_database,
        column = column,
        adduct.table = adduct.table,
        ms1.match.ppm = ms1.match.ppm,
        rt.match.tol = rt.match.tol,
        mz.ppm.thr = mz.ppm.thr,
        threads = threads
      )
    
    # Remove pathways that contain no compound
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
    
    # Remove pathway databases that have no pathways left after filtering
    
    if (length(pathway_database) == 0) {
      stop("All pathway database length is zero.\n")
    }
    
    # Combine peaks into metabolite classes based on retention time similarity
    ####
    message("Identifying metablite classes...\n")
    annotation_table_all <-
      unique(annotation_table_all$Lab.ID) %>%
      furrr::future_map(
        .f = function(temp_id) {
          x = annotation_table_all %>%
            dplyr::filter(Lab.ID == temp_id)
          
          x <- x %>%
            dplyr::mutate(rt = as.numeric(rt)) %>%
            dplyr::arrange(rt)
          
          rt_class <-
            group_peaks_rt(rt = x$rt, rt.tol = rt.match.tol) %>%
            dplyr::arrange(rt)
          
          rt_class <- paste(x$Lab.ID[1], rt_class$class, sep = "_")
          
          x <-
            data.frame(x,
                       compound_class = rt_class,
                       stringsAsFactors = FALSE)
          
          x <-
            unique(x$compound_class) %>%
            purrr::map(function(y) {
              z <-
                x[x$compound_class == y, , drop = FALSE]
              score <- score_peak_group(z)
              z <- data.frame(z, score, stringsAsFactors = FALSE)
              z
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame()
          x
        },
        .progress = TRUE
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    # ####if our method can be consistent with standard method
    # annotation_standard <-
    #   urine_metabolomics_data@variable_info %>%
    #   dplyr::filter(!is.na(SS)) %>%
    #   dplyr::filter(!is.na(KEGG.ID)) %>%
    #   dplyr::filter(KEGG.ID %in% node_data$KEGG_ID) %>%
    #   dplyr::filter(KEGG.ID %in% metabolite_database@spectra.info$KEGG.ID)
    #
    # compared_result <-
    #   annotation_standard %>%
    #   dplyr::select(variable_id, Compound.name, KEGG.ID, SS) %>%
    #   dplyr::left_join(annotation_table_all[, c("variable_id",
    #                                             "Compound.name",
    #                                             "KEGG.ID",
    #                                             "compound_class",
    #                                             "score")], by = "variable_id")
    #
    # write.csv(compared_result, "compared_result.csv")
    #
    # temp <-
    #   unique(compared_result$variable_id) %>%
    #   purrr::map(function(x) {
    #     temp <-
    #       compared_result %>%
    #       dplyr::filter(variable_id == x)
    #
    #     accurate_score <-
    #       temp %>%
    #       dplyr::filter(KEGG.ID.x == KEGG.ID.y) %>%
    #       pull(score)
    #     max_score <- max(temp$score)
    #     c(accurate_score = accurate_score,
    #       max_score = max_score,
    #       diff_score = accurate_score - max_score)
    #   }) %>%
    #   do.call(rbind, .) %>%
    #   as.data.frame() %>%
    #   dplyr::filter(diff_score != 0)
    #
    # temp %>%
    #   dplyr::filter(max_score >= 80) %>%
    #   pull(accurate_score)
    #
    #
    # idx = temp %>% lapply(length) %>% unlist() %>% `==`(0) %>% which()
    #
    # temp2 = urine_metabolomics_data@variable_info %>%
    #   dplyr::filter(variable_id %in% unique(compared_result$variable_id)[idx]) %>%
    #   dplyr::left_join(metabolite_database@spectra.info[,c("KEGG.ID", "Formula")],
    #                    by = "KEGG.ID")
    # write.csv(temp2, "temp2.csv")
    #
    # x = masstools::sum_formula(formula = "C9H18NO4", adduct = "M+H")
    # Rdisop::getMass(molecule =  Rdisop::getMolecule("C9H18NO4"))
    # (156.040183 - 156.0514)*10^6/400
    
    
    # Remove redundant metabolite annotations
    # calculate_redundance(annotation_table = annotation_table_all)
    annotation_table_all <-
      remove_redundancy(annotation_table = annotation_table_all)
    
    # calculate_redundance(annotation_table = annotation_table_all2)
    
    # Filter the marker feature table to include only annotated markers
    feature_table_marker <-
      feature_table_marker %>%
      dplyr::filter(variable_id %in% annotation_table_all$variable_id)
    
    annotation_table_marker <-
      annotation_table_all %>%
      dplyr::filter(variable_id %in% feature_table_marker$variable_id)
    
    message("Detecting metabolic modules from metabolic reaction network.\n")
    
    ###detected_metabolites are IDs of metabolites that matched with the metabolic features
    detected_metabolites <-
      unique(annotation_table_marker$KEGG.ID)
    
    ###hidden metabolites are metabolites that are not matched with the metabolic features,
    ###and are connected to the detected_metabolites in the metabolic network within three steps
    message("Detecting hidden metabolites...\n")
    
    if (include_hidden_metabolites) {
      hidden_metabolites <-
        get_hidden_metabolites(
          metabolic_network = metabolic_network,
          detected_metabolites = detected_metabolites,
          threads = threads,
          max.reaction = 3
        )
      
      hidden_metabolites <-
        hidden_metabolites %>%
        unlist() %>%
        unique()
    } else{
      hidden_metabolites <- NULL
    }
    
    total_metabolites <-
      c(detected_metabolites, hidden_metabolites)
    
    ####extract subnetwork
    sub_metabolic_network <-
      metabolic_network %>%
      tidygraph::activate("nodes") %>%
      dplyr::filter(name %in% total_metabolites) %>%
      dplyr::mutate(node_class = ifelse(name %in% detected_metabolites, "Detected", "Hidden")) %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(degree2 = tidygraph::centrality_degree())
    
    # sub_metabolic_network2 <-
    # igraph::subgraph(graph = metabolic_network,
    #                  vids = total_metabolites) %>%
    #   tidygraph::as_tbl_graph() %>%
    #   dplyr::mutate(node_class = ifelse(name %in% detected_metabolites, "Detected", "Hidden")) %>%
    #   tidygraph::activate(what = "nodes") %>%
    #   dplyr::mutate(degree2 = tidygraph::centrality_degree())
    
    ###metabolic_module detection
    #####detecting metabolic modules
    message("Detecting metabolic modules...\n")
    
    metabolic_modules <-
      identify_metabolic_modules(
        sub_metabolic_network = sub_metabolic_network,
        detected_metabolites = detected_metabolites,
        hidden_metabolites = hidden_metabolites
      )
    
    ###calculate activity scores and impacts of metabolic modules
    message("Calculating quality scores of metabolic modules.\n")
    
    module_quality_score <-
      calculate_activity_socre(
        metabolic_modules = metabolic_modules,
        detected_metabolites = detected_metabolites,
        hidden_metabolites = hidden_metabolites,
        sub_metabolic_network = sub_metabolic_network,
        threads = threads
      )
    
    # module_quality_score <-
    #   metabolic_modules %>%
    #   purrr::map(function(x) {
    #     calculate_module_quality(graph = sub_metabolic_network, nodes = x)$conductance
    #   }) %>%
    #   unlist()
    
    ##Module impact
    message("Calculating impacts of metabolic modules...\n")
    
    module_impact <-
      lapply(metabolic_modules, function(temp_group) {
        temp_graph <-
          igraph::subgraph(graph = sub_metabolic_network, v = temp_group)
        centrality <- calculate_centrality(graph = temp_graph, type = "d")
        temp_idx <- which(temp_group %in% detected_metabolites)
        module_impact <- sum(centrality[temp_idx]) / sum(centrality)
        module_impact
      })
    
    module_impact <- unlist(module_impact)
    
    ####Hub Dominance Score (HDS) or Dominant Edge Ratio
    ##DER = \frac{\sum_{\text{top } K} e_i}{|E|}
    # where:
    #  K  is a small fraction (e.g., top 10-20%) of nodes sorted by degree.
    # 	 e_i  is the number of edges connected to those nodes.
    # 	 |E|  is the total number of edges.
    #
    # Rationale: If the top nodes account for almost all edges, DER → 1.
    message("Calculating Dominant Edge Ratio (DER) of metabolic modules...\n")
    
    module_dominant_edge_rate <-
      lapply(metabolic_modules, function(temp_group) {
        temp_graph <-
          igraph::subgraph(graph = sub_metabolic_network, v = temp_group)
        
        all_degree <- igraph::degree(temp_graph)
        
        # Sort nodes by degree (descending)
        sorted_degrees <- sort(all_degree, decreasing = TRUE)
        
        # Define top K nodes (e.g., top 20% most connected nodes)
        # K <- ceiling(0.2 * length(V(temp_graph)))  # 20% of nodes
        # top_nodes <- names(sorted_degrees)[1:K]
        
        top_nodes <-
          names(all_degree[all_degree > quantile(sorted_degrees, probs = 0.9)])
        if (length(top_nodes) == 0) {
          top_nodes <- names(sorted_degrees)[1]
        }
        
        
        all_edges <-
          igraph::as_edgelist(temp_graph)
        
        # Count edges where at least one node is in the top nodes
        dominant_edges <-
          sum(all_edges[, 1] %in% top_nodes |
                all_edges[, 2] %in% top_nodes)
        
        # Compute Dominant Edge Ratio (DER)
        total_edges <- igraph::gsize(temp_graph)
        dominant_edges / total_edges
      }) %>%
      unlist()
    
    ##get NULL distribution
    message("Identify null distributation of metabolic module quality scores...\n")
    
    null_quality_score <-
      generate_null_activity_score_distribution(
        annotation_table_all = annotation_table_all,
        feature_table_marker = feature_table_marker,
        metabolite_database = metabolite_database,
        metabolic_network = metabolic_network,
        permutation_times = 5,
        threads = threads,
        adduct.table = adduct.table,
        column = column,
        ms1.match.ppm = ms1.match.ppm,
        mz.ppm.thr = mz.ppm.thr,
        include_hidden_metabolites = include_hidden_metabolites
      )
    
    # null_quality_score <-
    #   generate_null_conductance_distribution(
    #     annotation_table_all = annotation_table_all,
    #     feature_table_marker = feature_table_marker,
    #     metabolite_database = metabolite_database,
    #     metabolic_network = metabolic_network,
    #     permutation_times = 5,
    #     threads = threads,
    #     adduct.table = adduct.table,
    #     column = column,
    #     ms1.match.ppm = ms1.match.ppm,
    #     mz.ppm.thr = mz.ppm.thr,
    #     include_hidden_metabolites = include_hidden_metabolites
    #   )
    
    null_quality_score <-
      unlist(null_quality_score)
    
    # mean(module_quality_score)
    # mean(null_quality_score)
    # median(mean(module_quality_score))
    # median(mean(null_quality_score))
    
    ####calculate the p values for each metabolic module
    module_quality_score <- module_quality_score + 0.00000001
    
    ###
    para <- MASS::fitdistr(x = null_quality_score, densfun = "gamma")[[1]]
    # standard.distribution <- rgamma(n = length(null_quality_score), shape = para[1], rate = para[2])
    standard.distribution <- rgamma(n = 100000,
                                    shape = para[1],
                                    rate = para[2])
    # ks.test(x = null_quality_score, y = standard.distribution)
    
    cdf <- ecdf(x = standard.distribution)
    
    metabolic_module_p <- 1 - cdf(module_quality_score)
    
    #order module according to p value
    # metabolic_modules <- metabolic_modules[order(metabolic_module_p)]
    # module_impact <- module_impact[order(metabolic_module_p)]
    # metabolic_module_p <- metabolic_module_p[order(metabolic_module_p)]
    # module_dominant_edge_rate <- module_dominant_edge_rate[order(metabolic_module_p)]
    
    ## get information of each metabolic module
    metabolic_module_result <-
      lapply(metabolic_modules, function(x) {
        temp_graph <-
          igraph::subgraph(graph = metabolic_network, v = x)
        
        detected_metabolite_number <- sum(x %in% detected_metabolites)
        hidden_metabolite_number <- sum(x %in% hidden_metabolites)
        total_metabolite_number <- length(x)
        
        i.id <- x[which(x %in% detected_metabolites)]
        h.id <- x[which(x %in% hidden_metabolites)]
        
        detected_metabolite_name <-
          data.frame(KEGG.ID = i.id) %>%
          dplyr::left_join(
            annotation_table_all %>%
              dplyr::select(KEGG.ID, Compound.name) %>%
              dplyr::distinct(KEGG.ID, .keep_all = TRUE),
            by = "KEGG.ID"
          ) %>%
          pull(Compound.name) %>%
          paste(collapse = "{}")
        
        hidden_metabolite_name <-
          data.frame(KEGG.ID = h.id) %>%
          dplyr::left_join(
            annotation_table_all %>%
              dplyr::select(KEGG.ID, Compound.name) %>%
              dplyr::distinct(KEGG.ID, .keep_all = TRUE),
            by = "KEGG.ID"
          ) %>%
          pull(Compound.name) %>%
          paste(collapse = "{}")
        
        total_metabolite_name <-
          data.frame(KEGG.ID = x) %>%
          dplyr::left_join(
            annotation_table_all %>%
              dplyr::select(KEGG.ID, Compound.name) %>%
              dplyr::distinct(KEGG.ID, .keep_all = TRUE),
            by = "KEGG.ID"
          ) %>%
          pull(Compound.name) %>%
          paste(collapse = "{}")
        
        i.id <- paste(i.id, collapse = "{}")
        h.id <- paste(h.id, collapse = "{}")
        t.id <- paste(x, collapse = "{}")
        
        temp_result <-
          c(
            total_metabolite_number,
            detected_metabolite_number,
            hidden_metabolite_number,
            t.id,
            i.id,
            h.id,
            total_metabolite_name,
            detected_metabolite_name,
            hidden_metabolite_name
          )
        return(temp_result)
      })
    
    metabolic_module_result <-
      do.call(rbind, metabolic_module_result)
    
    metabolic_module_result <-
      data.frame(
        module_impact,
        module_dominant_edge_rate,
        metabolic_module_p,
        metabolic_module_result,
        stringsAsFactors = FALSE
      ) %>%
      dplyr::arrange(metabolic_module_p) %>%
      dplyr::mutate(metabolic_module_name = paste("metabolic_module", 1:nrow(.), sep = "_")) %>%
      dplyr::select(metabolic_module_name, dplyr::everything())
    
    colnames(metabolic_module_result) <- c(
      "Name",
      "Impact",
      "Dominant_edge_rate",
      "p_value",
      "Total_metabolite_number",
      "Detected_metabolite_number",
      "Hidden_metabolite_number",
      "Total_metabolite_id",
      "Detected_metabolite_id",
      "hidden_metabolite_id",
      "Total_metabolite_name",
      "Detected_metabolite_name",
      "hidden_metabolite_name"
    )
    
    metabolic_module_result$Total_metabolite_number <-
      as.numeric(metabolic_module_result$Total_metabolite_number)
    
    metabolic_module_result$Detected_metabolite_number <-
      as.numeric(metabolic_module_result$Detected_metabolite_number)
    
    metabolic_module_result$Hidden_metabolite_number <-
      as.numeric(metabolic_module_result$Hidden_metabolite_number)
    
    #####extract the final dysregulated networks
    final_metabolite_id <-
      metabolic_module_result %>%
      dplyr::filter(p_value < 0.05) %>%
      dplyr::pull(Total_metabolite_id) %>%
      stringr::str_split(pattern = "\\{\\}") %>%
      unlist() %>%
      unique()
    
    annotation_table_final <-
      annotation_table_all %>%
      dplyr::filter(
        variable_id %in% feature_table_marker$variable_id &
          KEGG.ID %in% final_metabolite_id
      )
    
    # ####not in the workflow
    # calculate_redundance(annotation_table_final)
    #
    # annotation_standard <-
    #   urine_metabolomics_data@variable_info %>%
    #   dplyr::filter(!is.na(SS)) %>%
    #   dplyr::filter(!is.na(KEGG.ID)) %>%
    #   dplyr::filter(KEGG.ID %in% node_data$KEGG_ID) %>%
    #   dplyr::filter(KEGG.ID %in% metabolite_database@spectra.info$KEGG.ID)
    #
    # dim(annotation_standard)
    #
    # compared_result2 <-
    #   annotation_standard %>%
    #   dplyr::select(variable_id, Compound.name, KEGG.ID, SS) %>%
    #   dplyr::left_join(annotation_table_final[, c("variable_id",
    #                                               "Compound.name",
    #                                               "KEGG.ID",
    #                                               "compound_class",
    #                                               "score")], by = "variable_id")
    #
    # write.csv(compared_result2, file = "compared_result2.csv")
    #
    # unique(compared_result2$variable_id) %>%
    #   purrr::map(function(x) {
    #     temp <-
    #       compared_result2 %>%
    #       dplyr::filter(variable_id == x)
    #
    #     if (all(is.na(temp$KEGG.ID.y))) {
    #       return("No annotation")
    #     }
    #
    #     if (unique(temp$KEGG.ID.x) %in% temp$KEGG.ID.y) {
    #       return("Correct annotation")
    #     } else{
    #       return("Incorrect annotation")
    #     }
    #
    #   }) %>%
    #   unlist() %>%
    #   table()
    #
    # compared_result3 <-
    #   compared_result2 %>%
    #   dplyr::filter(!is.na(KEGG.ID.y))
    #
    # write.csv(compared_result3, file = "compared_result3.csv")
    
    final_metabolic_network <-
      metabolic_network %>%
      tidygraph::activate("nodes") %>%
      dplyr::filter(name %in% final_metabolite_id) %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(node_class = ifelse(name %in% detected_metabolites, "Detected", "Hidden")) %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(degree2 = tidygraph::centrality_degree())
    
    final_node_data <-
      tidygraph::as_tibble(final_metabolic_network, active = "nodes") %>%
      dplyr::select(name, KEGG_ID, HMDB_ID, Compound_name, node_class)
    
    final_edge_data <-
      tidygraph::as_tibble(final_metabolic_network, active = "edges") %>%
      dplyr::mutate(from = from_compound_KEGG_ID, to = to_compound_KEGG_ID) %>%
      dplyr::select(from, to) %>%
      dplyr::mutate(edge_class = "metabolic_reaction")
    
    final_edge_data2 <-
      annotation_table_final %>%
      dplyr::filter(isotope == "[M]") %>%
      dplyr::select(variable_id, KEGG.ID) %>%
      dplyr::rename(from = variable_id, to = KEGG.ID) %>%
      dplyr::mutate(edge_class = "feature_metabolite")
    
    final_node_data2 <-
      data.frame(name = unique(final_edge_data2$from)) %>%
      dplyr::mutate(
        KEGG_ID = NA,
        HMDB_ID = NA,
        Compound_name = NA,
        node_class = "Feature"
      )
    
    ###add feature information to the final metabolic network
    final_node_data <-
      rbind(final_node_data, final_node_data2)
    
    final_edge_data <-
      rbind(final_edge_data, final_edge_data2)
    
    final_node_data <-
      final_node_data %>%
      dplyr::filter(name %in% unique(c(
        final_edge_data$from, final_edge_data$to
      ))) %>%
      dplyr::distinct(name, .keep_all = TRUE)
    
    final_edge_data <-
      final_edge_data %>%
      dplyr::filter(from %in% final_node_data$name &
                      to %in% final_node_data$name)
    
    final_metabolic_network <-
      tidygraph::tbl_graph(nodes = final_node_data,
                           edges = final_edge_data,
                           directed = FALSE)
    
    ####Pathway enrichment for all the metabolites
    kegg_id <-
      annotation_table_final$KEGG.ID
    
    message("Enriching pathways for dysregulated metabolic network...\n")
    
    enriched_pathways <-
      enrich_pathways(
        query_id = kegg_id,
        query_type = "compound",
        id_type = "KEGG",
        pathway_database = pathway_database,
        p_cutoff = 0.05,
        p_adjust_method = "fdr",
        threads = threads
      )
    
    message("Enriching pathways for each metabolic module...\n")
    
    enriched_pathways_list <-
      vector("list", length = nrow(metabolic_module_result))
    
    names(enriched_pathways_list) <- metabolic_module_result$Name
    
    for (i in 1:sum(metabolic_module_result$p_value < 0.05)) {
      cat(i, " ")
      x <-
        metabolic_module_result$Total_metabolite_id[i] %>%
        stringr::str_split(pattern = "\\{\\}") %>%
        unlist() %>%
        unique()
      enriched_pathways_list[[i]] <-
        enrich_pathways(
          query_id = x,
          query_type = "compound",
          id_type = "KEGG",
          pathway_database = pathway_database,
          p_cutoff = 0.05,
          p_adjust_method = "fdr",
          threads = threads
        )
    }
    
    return(
      list(
        dysregulated_metabolic_module = metabolic_module_result,
        dysregulated_metabolic_network = final_metabolic_network,
        annotation_table = annotation_table_final,
        enriched_pathways = enriched_pathways,
        enriched_pathways_list = enriched_pathways_list
      )
    )
    
  }






#' Identify Metabolic Modules
#'
#' This function identifies metabolic modules from a metabolic reaction network.
#' It applies a clustering algorithm to detect modules based on network structure.
#' Large modules are recursively split, and weakly connected hidden metabolites are removed.
#'
#' @param sub_metabolic_network A metabolic network represented as a `tbl_graph` object.
#' @param detected_metabolites A character vector of detected metabolite IDs.
#' @param hidden_metabolites A character vector of hidden metabolite IDs.
#'
#' @return A list of metabolic modules, where each module is a character vector
#'   containing the metabolite IDs within the module.
#'
#' @export

identify_metabolic_modules <-
  function(sub_metabolic_network,
           detected_metabolites,
           hidden_metabolites) {
    # Apply Walktrap clustering to detect metabolic modules
    metabolic_modules <- igraph::cluster_walktrap(graph = sub_metabolic_network)
    
    # Extract module membership information
    node_name <- metabolic_modules$names
    membership <- metabolic_modules$membership
    membership_counts <- table(membership)
    
    # Group nodes by module membership
    group <- lapply(names(membership_counts), function(x) {
      node_name[which(membership == as.numeric(x))]
    })
    
    # Remove modules with fewer than 3 nodes
    group <- group[unlist(lapply(group, length)) > 3]
    
    # Recursively split modules containing more than 100 nodes
    large_module_idx <- which(unlist(lapply(group, length)) > 100)
    
    while (length(large_module_idx) > 0) {
      large_group <- group[large_module_idx]
      group <- group[-large_module_idx]
      
      new_group <- lapply(large_group, function(x) {
        # Extract subnetwork for the large module
        temp_graph <- sub_metabolic_network %>%
          tidygraph::activate("nodes") %>%
          dplyr::filter(name %in% x)
        
        # Apply clustering again within the subnetwork
        temp_metabolic_modules <- igraph::cluster_walktrap(graph = temp_graph)
        temp_node_name <- temp_metabolic_modules$names
        temp_membership <- temp_metabolic_modules$membership
        temp_membership_counts <- table(temp_membership)
        
        # Split into smaller groups
        temp_group <- lapply(names(temp_membership_counts), function(y) {
          temp_node_name[which(temp_membership == as.numeric(y))]
        })
        
        # Remove small groups with fewer than 3 nodes
        temp_group <- temp_group[unlist(lapply(temp_group, length)) > 3]
        temp_group
      })
      
      new_group <- unlist(new_group, recursive = FALSE)
      group <- c(group, new_group)
      
      # Identify any remaining large modules
      large_module_idx <- which(unlist(lapply(new_group, length)) > 100)
    }
    
    # Remove weakly connected groups where the number of edges is less than the number of nodes
    remove_idx <- which(unlist(lapply(group, function(x) {
      temp_graph <- sub_metabolic_network %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(name %in% x)
      
      node_count <- igraph::gorder(temp_graph)  # Number of nodes
      edge_count <- igraph::gsize(temp_graph)  # Number of edges
      degree_values <- igraph::degree(graph = temp_graph, v = x)
      
      # Remove groups where edges are fewer than nodes and all nodes have degree < 3
      if (node_count - edge_count > 0 & all(degree_values < 3)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    })))
    
    if (length(remove_idx) > 0) {
      group <- group[-remove_idx]
    }
    
    # Remove hidden metabolites with degree < 2
    group <- lapply(group, function(x) {
      temp_graph <- sub_metabolic_network %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(name %in% x)
      
      degree_values <- igraph::degree(graph = temp_graph, v = x)
      
      metabolite_class <- sapply(names(degree_values), function(y) {
        ifelse(y %in% detected_metabolites, "Detected", "Hidden")
      })
      
      remove_idx <- which(degree_values == 1 &
                            metabolite_class == "Hidden")
      
      while (length(remove_idx) > 0) {
        x <- x[-remove_idx]
        temp_graph <- sub_metabolic_network %>%
          tidygraph::activate("nodes") %>%
          dplyr::filter(name %in% x)
        
        degree_values <- igraph::degree(graph = temp_graph, v = x)
        metabolite_class <- sapply(names(degree_values), function(y) {
          ifelse(y %in% detected_metabolites, "Detected", "Hidden")
        })
        
        remove_idx <- which(degree_values == 1 &
                              metabolite_class == "Hidden")
      }
      
      return(x)
    })
    
    # Remove groups that still have fewer than 3 nodes
    group <- group[unlist(lapply(group, length)) > 3]
    
    return(group)
  }





#' Calculate Node Centrality in a Graph
#'
#' This function computes the centrality of nodes in a graph using either degree centrality or betweenness centrality.
#'
#' @param graph An `igraph` object representing the network.
#' @param type A character string specifying the type of centrality to calculate. Options are "degree" (default) or "betweenness".
#'
#' @return A named vector of centrality values for each node in the graph.
#'
#' @export
#'
calculate_centrality <-
  function(graph, type = c("degree", "betweenness")) {
    # Ensure that the specified centrality type is valid
    type <- match.arg(type)
    
    # Extract node names from the graph
    node <- igraph::V(graph)$name
    
    # Compute degree centrality if selected
    if (type == "degree") {
      node_degree <- igraph::degree(graph = graph, v = node)
      return(node_degree)
    }
    
    # Compute betweenness centrality
    node_betweenness <- lapply(node, function(x) {
      i <- match(x, node)
      if (i == length(node))
        return(NULL)
      
      from <- x
      to <- node[(i + 1):length(node)]
      
      # Compute shortest paths between the node and others
      temp <- igraph::shortest_paths(graph = graph,
                                     from = from,
                                     to = to)[[1]]
      
      # Extract nodes that lie on shortest paths
      temp <- lapply(temp, function(path) {
        if (length(path) <= 2)
          return(NULL)
        names(path)[-c(1, length(path))]
      })
      
      unlist(temp)
    })
    
    node_betweenness <- unlist(node_betweenness)
    node_betweenness <- table(node_betweenness)
    
    # Ensure that all nodes are included in the result
    missing_nodes <- setdiff(node, names(node_betweenness))
    if (length(missing_nodes) > 0) {
      add <- rep(0, length(missing_nodes))
      names(add) <- missing_nodes
      node_betweenness <- c(node_betweenness, add)
    }
    
    # Ensure correct ordering of results
    node_betweenness <- node_betweenness[match(node, names(node_betweenness))]
    
    return(node_betweenness)
  }



#' Identify Hidden Metabolites in a Metabolic Network
#'
#' This function identifies hidden metabolites that are connected to detected metabolites
#' within a specified number of reaction steps in a metabolic network.
#'
#' @param metabolic_network An `igraph` object representing the metabolic reaction network.
#' @param detected_metabolites A character vector of detected metabolite IDs.
#' @param threads Integer; number of parallel threads to use. Default is `3`.
#' @param max.reaction Integer; maximum number of reaction steps allowed to define a hidden metabolite. Default is `3`.
#'
#' @return A list of hidden metabolites connected to detected metabolites within the specified reaction steps.
#'
#' @export
#'
get_hidden_metabolites <-
  function(metabolic_network,
           detected_metabolites,
           threads = 8,
           max.reaction = 3) {
    # Helper function to find hidden metabolites for each detected metabolite
    temp_fun_hidden_metabolites <- function(name,
                                            metabolic_network,
                                            detected_metabolites,
                                            max.reaction) {
      options(warn = -1)  # Suppress warnings
      
      # Identify remaining detected metabolites for distance calculation
      idx <- match(name, detected_metabolites)
      to <- detected_metabolites[-c(1:idx)]
      
      # Compute shortest path distances from the current metabolite
      distance <- igraph::distances(
        graph = metabolic_network,
        v = name,
        to = to,
        mode = "all"
      )[1, ]
      
      # Set infinite distances to exceed max.reaction threshold
      distance[is.infinite(distance)] <- max.reaction + 10
      to <- to[which(distance <= max.reaction)]
      
      if (length(to) == 0) {
        return(NULL)  # No valid hidden metabolites found
      }
      
      # Extract shortest paths between the metabolite and its connected nodes
      result <- igraph::shortest_paths(
        graph = metabolic_network,
        from = name,
        to = to,
        mode = "all"
      )[[1]]
      
      # Extract node names from shortest paths and remove already detected metabolites
      result <- unique(unlist(lapply(result, names)))
      result <- setdiff(result, detected_metabolites)
      
      return(result)
    }
    
    # Select parallel processing backend based on OS
    if (masstools::get_os() == "windows") {
      bpparam <- BiocParallel::SnowParam(workers = threads, progressbar = TRUE)
    } else {
      bpparam <- BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)
    }
    
    # Identify hidden metabolites using parallel processing
    hidden_metabolites <- BiocParallel::bplapply(
      detected_metabolites[-length(detected_metabolites)],
      FUN = temp_fun_hidden_metabolites,
      BPPARAM = bpparam,
      metabolic_network = metabolic_network,
      detected_metabolites = detected_metabolites,
      max.reaction = max.reaction
    )
    
    return(hidden_metabolites)
  }



#' Calculate Activity Scores for Metabolic Modules
#'
#' This function calculates activity scores for metabolic modules based on their network properties,
#' including modularity, detected metabolite contributions, and overall network connectivity.
#'
#' @param metabolic_modules A list of metabolic modules, each containing metabolite IDs.
#' @param detected_metabolites A character vector of detected metabolite IDs.
#' @param hidden_metabolites A character vector of hidden metabolite IDs.
#' @param sub_metabolic_network An `igraph` object representing the metabolic subnetwork.
#' @param threads Integer; number of parallel threads to use. Default is `3`.
#'
#' @return A numeric vector of activity scores for each metabolic module.
#'
#' @export
#'
calculate_activity_score <-
  function(metabolic_modules,
           detected_metabolites,
           hidden_metabolites,
           sub_metabolic_network,
           threads = 8) {
    # Helper function to compute activity scores for a given metabolic module
    temp_fun_activity_scores <- function(x,
                                         detected_metabolites,
                                         hidden_metabolites,
                                         total_metabolite_number,
                                         sub_metabolic_network) {
      options(warn = -1)  # Suppress warnings
      
      # Extract the subgraph for the given metabolic module
      temp_graph <- igraph::subgraph(graph = sub_metabolic_network, v = x)
      detected_number <- sum(x %in% detected_metabolites)
      hidden_number <- sum(x %in% hidden_metabolites)
      
      total_edges <- length(igraph::E(sub_metabolic_network))  # Total edges in network
      module_edges <- length(igraph::E(temp_graph))  # Edges within the module
      
      node_names <- names(igraph::V(temp_graph))
      
      # Compute modularity-based contribution
      value <- sapply(node_names, function(node) {
        temp_value <- sapply(node_names, function(x) {
          degree_x <- igraph::degree(graph = sub_metabolic_network, v = x)
          degree_node <- igraph::degree(graph = sub_metabolic_network, v = node)
          
          (degree_x * degree_node) / (4 * total_edges^2)
        })
        sum(temp_value)
      })
      
      modularity <- module_edges / total_edges - sum(value)
      
      # Compute activity score
      q <- modularity * sqrt(total_metabolite_number / length(x))
      activity <- detected_number * q / length(x)
      
      return(activity)
    }
    
    # Choose appropriate parallelization method based on OS
    if (masstools::get_os() == "windows") {
      bpparam <- BiocParallel::SnowParam(workers = threads, progressbar = TRUE)
    } else {
      bpparam <- BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)
    }
    
    # Compute activity scores for all metabolic modules in parallel
    activity_scores <- BiocParallel::bplapply(
      metabolic_modules,
      temp_fun_activity_scores,
      BPPARAM = bpparam,
      detected_metabolites = detected_metabolites,
      hidden_metabolites = hidden_metabolites,
      total_metabolite_number = length(detected_metabolites),
      sub_metabolic_network = sub_metabolic_network
    )
    
    return(unlist(activity_scores))
  }



#' Annotate Metabolites for Feature Pathway Analysis (FPA)
#'
#' This function performs metabolite annotation for a given feature table using an MS1-based approach.
#' It processes both positive and negative ionization modes separately and extracts annotation results.
#'
#' @param feature_table A data frame containing feature data, including `variable_id`, `mz`, `rt`, and `polarity` columns.
#' @param metabolite_database A metabolite database object used for annotation.
#' @param column A character string specifying the chromatography column type. Options are `"rp"` (reverse phase) or `"hilic"` (hydrophilic interaction chromatography). Default is `"rp"`.
#' @param adduct.table Optional data frame specifying possible adducts. Default is `NULL`.
#' @param ms1.match.ppm Numeric; mass tolerance for MS1 matching in parts per million (ppm). Default is `10`.
#' @param rt.match.tol Numeric; retention time tolerance threshold for matching in seconds. Default is `5`.
#' @param mz.ppm.thr Numeric; M/Z tolerance threshold for filtering in ppm. Default is `400`.
#' @param threads Integer; number of parallel threads to use for annotation. Default is `3`.
#'
#' @return A data frame containing the annotated metabolites with polarity information.
#'
#' @export
#'
annotate_metabolites_fpa <-
  function(feature_table,
           metabolite_database,
           column = c("rp", "hilic"),
           adduct.table = NULL,
           ms1.match.ppm = 10,
           rt.match.tol = 5,
           mz.ppm.thr = 400,
           threads = 3) {
    column <- match.arg(column)  # Ensure valid column type
    
    ##### Separate features by polarity #####
    feature_table_pos <- feature_table %>%
      dplyr::filter(polarity == "positive")
    feature_table_neg <-
      feature_table %>% dplyr::filter(polarity == "negative")
    
    ##### Create mass dataset objects for positive and negative modes #####
    # Positive mode processing
    if (nrow(feature_table_pos) == 0) {
      object_pos <- NULL
    } else {
      expression_data_pos <- data.frame(sample1 = 1:nrow(feature_table_pos))
      rownames(expression_data_pos) <- feature_table_pos$variable_id
      sample_info_pos <- data.frame(sample_id = "sample1", class = "Subject")
      
      object_pos <-
        massdataset::create_mass_dataset(
          expression_data = expression_data_pos,
          sample_info = sample_info_pos,
          variable_info = feature_table_pos
        )
    }
    
    # Negative mode processing
    if (nrow(feature_table_neg) == 0) {
      object_neg <- NULL
    } else {
      expression_data_neg <- data.frame(sample1 = 1:nrow(feature_table_neg))
      rownames(expression_data_neg) <- feature_table_neg$variable_id
      sample_info_neg <- data.frame(sample_id = "sample1", class = "Subject")
      
      object_neg <- massdataset::create_mass_dataset(
        expression_data = expression_data_neg,
        sample_info = sample_info_neg,
        variable_info = feature_table_neg
      )
    }
    
    ##### Perform metabolite annotation #####
    message("Annotating metabolites in positive mode...\n")
    if (!is.null(object_pos)) {
      object_pos <- metid::annotate_metabolites(
        object = object_pos,
        database = metabolite_database,
        based_on = "ms1",
        polarity = "positive",
        column = column,
        adduct.table = adduct.table,
        ms1.match.ppm = ms1.match.ppm,
        candidate.num = 100,
        threads = threads,
        mz.ppm.thr = mz.ppm.thr
      )
      annotation_table_pos <-
        massdataset::extract_annotation_table(object_pos) %>%
        dplyr::mutate(polarity = "positive")
      
      ####add feature mz, rt and annotation formula to it
      annotation_table_pos <-
        annotation_table_pos %>%
        dplyr::left_join(feature_table_pos[, c("variable_id", "mz", "rt")], by = c("variable_id" = "variable_id")) %>%
        dplyr::left_join(metabolite_database@spectra.info[, c("Lab.ID", "Formula")],
                         by = c("Lab.ID" = "Lab.ID")) %>%
        dplyr::select(variable_id, mz, rt, Formula, Adduct, dplyr::everything()) %>%
        dplyr::mutate(isotope = "[M]")
      
      ####isotope annotation for annotations
      message("Annotating isotopes in positive mode...\n")
      future::plan(future::multisession, workers = threads)
      system.time(
        isotope_pos <-
          furrr::future_map(
            .x = as.data.frame(t(annotation_table_pos)),
            .f = function(x) {
              adduct <-
                stringr::str_extract(x[5], "\\(.+\\)") %>%
                stringr::str_replace("\\(", "") %>%
                stringr::str_replace("\\)", "")
              
              mz = as.numeric(stringr::str_trim(x[2], side = "both"))
              rt = as.numeric(stringr::str_trim(x[3], side = "both"))
              
              temp_iso <- try(annotate_isotope(
                formula = stringr::str_trim(x[4], side = "both"),
                adduct = adduct,
                mz = mz,
                rt = rt,
                peak.mz = feature_table_pos$mz,
                peak.rt = feature_table_pos$rt,
                rt.tol = rt.match.tol,
                mz.tol = ms1.match.ppm,
                max.isotope = 3
              ),
              silent = TRUE)
              
              if (is(temp_iso, "try-error")) {
                return(NULL)
              }
              
              if (is.null(temp_iso)) {
                return(NULL)
              }
              
              # return(temp_iso)
              
              temp_iso <-
                cbind(feature_table_pos[temp_iso$peakIndex, ], temp_iso) %>%
                dplyr::select(-c(peakIndex))
              
              colnames(temp_iso) <- c("variable_id",
                                      "mz",
                                      "rt",
                                      "polarity",
                                      "mz.error",
                                      "isotope",
                                      "RT.error")
              
              x <- matrix(x, nrow = 1) %>% as.data.frame()
              colnames(x) <- colnames(annotation_table_pos)
              
              temp_iso$Compound.name <- x$Compound.name
              temp_iso$Lab.ID <- x$Lab.ID
              temp_iso$Adduct <- x$Adduct
              temp_iso$Formula <- x$Formula
              # temp_iso$rt <- x$rt
              temp_iso$CAS.ID <- x$CAS.ID
              temp_iso$HMDB.ID <- x$HMDB.ID
              temp_iso$KEGG.ID <- x$KEGG.ID
              temp_iso$Database <- x$Database
              temp_iso$ms2_files_id <- NA
              temp_iso$ms2_spectrum_id <- NA
              temp_iso$CE <- NA
              temp_iso$SS <- NA
              temp_iso$Total.score <- NA
              temp_iso$Level <- NA
              
              temp_iso$mz.match.score <-
                (ms1.match.ppm - temp_iso$mz.error) / ms1.match.ppm
              temp_iso$RT.match.score <-
                (rt.match.tol - temp_iso$RT.error) / rt.match.tol
              temp_iso %>%
                dplyr::select(colnames(x))
            },
            .progress = TRUE
          )
      )
      
      isotope_pos <- isotope_pos %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      annotation_table_pos <-
        rbind(annotation_table_pos, isotope_pos) %>%
        dplyr::arrange(Compound.name, isotope, rt)
      
    } else {
      annotation_table_pos <- NULL
    }
    
    message("Annotating metabolites in negative mode...\n")
    if (!is.null(object_neg)) {
      object_neg <- metid::annotate_metabolites(
        object = object_neg,
        database = metabolite_database,
        based_on = "ms1",
        polarity = "negative",
        column = column,
        adduct.table = adduct.table,
        ms1.match.ppm = ms1.match.ppm,
        candidate.num = 100,
        threads = threads,
        mz.ppm.thr = mz.ppm.thr
      )
      
      annotation_table_neg <-
        massdataset::extract_annotation_table(object_neg) %>%
        dplyr::mutate(polarity = "negative")
      
      ####add feature mz, rt and annotation formula to it
      annotation_table_neg <-
        annotation_table_neg %>%
        dplyr::left_join(feature_table_neg[, c("variable_id", "mz", "rt")], by = c("variable_id" = "variable_id")) %>%
        dplyr::left_join(metabolite_database@spectra.info[, c("Lab.ID", "Formula")],
                         by = c("Lab.ID" = "Lab.ID")) %>%
        dplyr::select(variable_id, mz, rt, Formula, Adduct, dplyr::everything()) %>%
        dplyr::mutate(isotope = "[M]")
      
      ####isotope annotation for annotations
      message("Annotating isotopes in negative mode...\n")
      future::plan(future::multisession, workers = threads)
      system.time(
        isotope_neg <-
          furrr::future_map(
            .x = as.data.frame(t(annotation_table_neg)),
            .f = function(x) {
              adduct <-
                stringr::str_extract(x[5], "\\(.+\\)") %>%
                stringr::str_replace("\\(", "") %>%
                stringr::str_replace("\\)", "")
              
              mz = as.numeric(stringr::str_trim(x[2], side = "both"))
              rt = as.numeric(stringr::str_trim(x[3], side = "both"))
              
              temp_iso <- try(annotate_isotope(
                formula = stringr::str_trim(x[4], side = "both"),
                adduct = adduct,
                mz = mz,
                rt = rt,
                peak.mz = feature_table_neg$mz,
                peak.rt = feature_table_neg$rt,
                rt.tol = rt.match.tol,
                mz.tol = ms1.match.ppm,
                max.isotope = 3
              ),
              silent = TRUE)
              
              if (is(temp_iso, "try-error")) {
                return(NULL)
              }
              
              if (is.null(temp_iso)) {
                return(NULL)
              }
              
              # return(temp_iso)
              
              temp_iso <-
                cbind(feature_table_neg[temp_iso$peakIndex, ], temp_iso) %>%
                dplyr::select(-c(peakIndex))
              
              colnames(temp_iso) <- c("variable_id",
                                      "mz",
                                      "rt",
                                      "polarity",
                                      "mz.error",
                                      "isotope",
                                      "RT.error")
              
              x <- matrix(x, nrow = 1) %>% as.data.frame()
              colnames(x) <- colnames(annotation_table_neg)
              
              temp_iso$Compound.name <- x$Compound.name
              temp_iso$Lab.ID <- x$Lab.ID
              temp_iso$Adduct <- x$Adduct
              temp_iso$Formula <- x$Formula
              # temp_iso$rt <- x$rt
              temp_iso$CAS.ID <- x$CAS.ID
              temp_iso$HMDB.ID <- x$HMDB.ID
              temp_iso$KEGG.ID <- x$KEGG.ID
              temp_iso$Database <- x$Database
              temp_iso$ms2_files_id <- NA
              temp_iso$ms2_spectrum_id <- NA
              temp_iso$CE <- NA
              temp_iso$SS <- NA
              temp_iso$Total.score <- NA
              temp_iso$Level <- NA
              
              temp_iso$mz.match.score <-
                (ms1.match.ppm - temp_iso$mz.error) / ms1.match.ppm
              temp_iso$RT.match.score <-
                (rt.match.tol - temp_iso$RT.error) / rt.match.tol
              temp_iso %>%
                dplyr::select(colnames(x))
            },
            .progress = TRUE
          )
      )
      
      isotope_neg <- isotope_neg %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      annotation_table_neg <-
        rbind(annotation_table_neg, isotope_neg) %>%
        dplyr::arrange(Compound.name, isotope, rt)
      
    } else {
      annotation_table_neg <- NULL
    }
    
    ##### Combine annotation tables #####
    annotation_table <-
      rbind(annotation_table_pos, annotation_table_neg)
    
    return(annotation_table)
  }



#' Generate Null Distribution of Activity Scores
#'
#' This function generates a null distribution of activity scores by performing random sampling
#' of marker metabolites and calculating activity scores through network analysis.
#'
#' @param annotation_table_all A data frame containing all annotated metabolites.
#' @param feature_table_marker A data frame containing marker features.
#' @param metabolite_database A metabolite database object used for annotation.
#' @param metabolic_network An `igraph` object representing the metabolic network.
#' @param permutation_times Integer; number of permutations to perform. Default is `20`.
#' @param threads Integer; number of parallel threads to use. Default is `8`.
#' @param adduct.table Optional data frame specifying possible adducts. Default is `NULL`.
#' @param column Character; chromatography column type. Options are `"rp"` or `"hilic"`. Default is `"rp"`.
#' @param ms1.match.ppm Numeric; MS1 matching tolerance in ppm. Default is `25`.
#' @param mz.ppm.thr Numeric; M/Z tolerance threshold in ppm. Default is `400`.
#' @param include_hidden_metabolites Logical; whether to include hidden metabolites in the network analysis. Default is `FALSE`.
#'
#' @return A list of numeric vectors representing null activity scores for each permutation.
#'
#' @export
#'
generate_null_activity_score_distribution <-
  function(annotation_table_all,
           feature_table_marker,
           metabolite_database,
           metabolic_network,
           permutation_times = 20,
           threads = 8,
           adduct.table = NULL,
           column = c("rp", "hilic"),
           ms1.match.ppm = 25,
           mz.ppm.thr = 400,
           include_hidden_metabolites = FALSE) {
    column <- match.arg(column)  # Ensure valid column type
    
    # Initialize list to store activity scores for each permutation
    null_activity_score <- vector(mode = "list", length = permutation_times)
    marker_number <- nrow(feature_table_marker)  # Number of markers
    
    for (i in seq_len(permutation_times)) {
      cat(i, " ")
      
      # Randomly sample markers for permutation
      if (nrow(feature_table_marker) / length(unique(annotation_table_all$variable_id)) < 0.2) {
        temp_marker_variable_id <- sample(
          x = setdiff(
            unique(annotation_table_all$variable_id),
            feature_table_marker$variable_id
          ),
          size = marker_number,
          replace = FALSE
        )
      } else{
        temp_marker_variable_id <- sample(
          x = unique(annotation_table_all$variable_id),
          size = marker_number,
          replace = FALSE
        )
      }
      
      
      # Extract annotation information for sampled markers
      temp_annotation_table_marker <- annotation_table_all %>%
        dplyr::filter(variable_id %in% temp_marker_variable_id)
      
      # Extract detected metabolites from annotation table
      temp_detected_metabolites <- unique(temp_annotation_table_marker$KEGG.ID)
      
      # Identify hidden metabolites connected to detected metabolites
      if (include_hidden_metabolites) {
        temp_hidden_metabolites <- get_hidden_metabolites(
          metabolic_network = metabolic_network,
          detected_metabolites = temp_detected_metabolites,
          threads = threads,
          max.reaction = 3
        ) %>% unlist() %>% unique()
      } else {
        temp_hidden_metabolites <- NULL
      }
      
      temp_total_metabolites <- c(temp_detected_metabolites, temp_hidden_metabolites)
      
      # Extract subnetwork containing detected and hidden metabolites
      temp_sub_metabolic_network <- metabolic_network %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(name %in% temp_total_metabolites) %>%
        dplyr::mutate(node_class = ifelse(name %in% temp_detected_metabolites, "Detected", "Hidden")) %>%
        tidygraph::mutate(degree2 = tidygraph::centrality_degree())
      
      # Identify metabolic modules in the subnetwork
      temp_metabolic_modules <- identify_metabolic_modules(
        sub_metabolic_network = temp_sub_metabolic_network,
        detected_metabolites = temp_detected_metabolites,
        hidden_metabolites = temp_hidden_metabolites
      )
      
      # Calculate activity scores for the identified metabolic modules
      temp_activity_score <- calculate_activity_score(
        metabolic_modules = temp_metabolic_modules,
        detected_metabolites = temp_detected_metabolites,
        hidden_metabolites = temp_hidden_metabolites,
        sub_metabolic_network = temp_sub_metabolic_network
      )
      
      # Store the computed activity score in the list
      null_activity_score[[i]] <- temp_activity_score
    }
    
    return(null_activity_score)
  }



#####g
#' Group Retention Time (RT) Peaks
#'
#' This function groups retention time (RT) values based on a specified tolerance,
#' clustering peaks that fall within the same RT range.
#'
#' @param rt A numeric vector of retention times.
#' @param rt.tol Numeric; retention time tolerance for grouping. Default is `10`.
#'
#' @return A data frame with two columns:
#'   - `rt`: The retention times of grouped peaks.
#'   - `class`: The assigned group ID.
#'
#' @export
#'
group_peaks_rt <-
  function(rt, rt.tol = 10) {
    # Identify clusters of retention times within the specified tolerance
    rt_class <- lapply(rt, function(x) {
      which(rt >= x & rt < x + rt.tol)
    })
    
    # Ensure each peak is uniquely assigned to a single class
    rt_class <- lapply(seq_along(rt_class)[-1], function(i) {
      setdiff(rt_class[[i]], unique(unlist(rt_class[1:(i - 1)])))
    }) %>% c(rt_class[1], .)
    
    # Remove empty groups
    rt_class <- rt_class[which(lapply(rt_class, length) != 0)]
    
    # Assign group names
    names(rt_class) <- seq_along(rt_class)
    
    # Convert grouping information into a data frame
    rt_class <- purrr::map2(
      .x = rt_class,
      .y = names(rt_class),
      .f = function(x, y) {
        data.frame(rt = rt[x],
                   class = y,
                   stringsAsFactors = FALSE)
      }
    ) %>% do.call(rbind, .)
    
    rownames(rt_class) <- NULL
    
    return(rt_class)
  }


#' Score a Peak Group Based on Adduct and Polarity Information
#'
#' This function assigns a score to a peak group based on the detected adduct types,
#' polarity, and isotopic information. Higher scores indicate a higher likelihood
#' of correct identification.
#'
#' @param peak_group A data frame containing peak group information, including `Adduct`, `polarity`, and `isotope` columns.
#'
#' @return A numeric score representing the confidence of the peak group assignment.
#'
#' @details
#' The scoring system is based on the following criteria:
#' - `+50` if the positive adduct is `(M+H)+`.
#' - `+20` if `(M+H)+` is present and is not the `[M]` isotope.
#' - `+20` if the polarity is positive but the adduct is not `(M+H)+`.
#' - `+10` if the polarity is positive, the adduct is not `(M+H)+`, and the isotope is not `[M]`.
#' - `+50` if the negative adduct is `(M-H)-`.
#' - `+20` if `(M-H)-` is present and is not the `[M]` isotope.
#' - `+20` if the polarity is negative but the adduct is not `(M-H)-`.
#' - `+10` if the polarity is negative, the adduct is not `(M-H)-`, and the isotope is not `[M]`.
#'
#' @export
#'
score_peak_group <-
  function(peak_group) {
    score <- 0
    
    ## This should be optimized in the future
    
    # Positive mode scoring
    if (any(peak_group$Adduct == "(M+H)+")) {
      score <- score + 50
    }
    
    if (any(peak_group$Adduct == "(M+H)+" &
            peak_group$isotope != "[M]")) {
      score <- score + 20
    }
    
    if (any(peak_group$polarity == "positive" &
            peak_group$Adduct != "(M+H)+")) {
      score <- score + 20
    }
    
    if (any(
      peak_group$polarity == "positive" &
      peak_group$Adduct != "(M+H)+" &
      peak_group$isotope != "[M]"
    )) {
      score <- score + 10
    }
    
    # Negative mode scoring
    if (any(peak_group$Adduct == "(M-H)-")) {
      score <- score + 50
    }
    
    if (any(peak_group$Adduct == "(M-H)-" &
            peak_group$isotope != "[M]")) {
      score <- score + 20
    }
    
    if (any(peak_group$polarity == "negative" &
            peak_group$Adduct != "(M-H)-")) {
      score <- score + 20
    }
    
    if (any(
      peak_group$polarity == "negative" &
      peak_group$Adduct != "(M-H)-" &
      peak_group$isotope != "[M]"
    )) {
      score <- score + 10
    }
    
    return(score)
  }





#' Calculate Redundancy in Metabolite Annotation
#'
#' This function calculates redundancy metrics in a metabolite annotation table.
#' It evaluates the number of distinct compound classes per compound and the number
#' of compounds assigned to each peak.
#'
#' @param annotation_table A data frame containing metabolite annotations, including `Lab.ID`, `compound_class`, and `variable_id`.
#'
#' @return A numeric vector of length two:
#'   - The first value represents the average number of compound classes per compound (`redundancy1`).
#'   - The second value represents the average number of compounds assigned to each peak (`redundancy2`).
#'
#' @export
#'
calculate_redundance <-
  function(annotation_table) {
    # Ensure the annotation table is in a consistent format
    annotation_table <- dplyr::bind_rows(annotation_table)
    
    # Redundancy 1: Average number of compound classes per compound
    redundancy1 <- annotation_table %>%
      dplyr::group_by(Lab.ID) %>%
      dplyr::distinct(compound_class) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::pull(n) %>%
      mean()
    
    # Redundancy 2: Average number of compounds assigned to each peak
    redundancy2 <- annotation_table %>%
      dplyr::group_by(variable_id) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::pull(n) %>%
      mean()
    
    return(c(redundancy1, redundancy2))
  }





#' Remove Redundant Annotations in Metabolite Data
#'
#' This function removes redundant metabolite annotations by filtering out low-confidence matches.
#' It iteratively refines the annotation table based on score thresholds.
#'
#' @param annotation_table A data frame containing metabolite annotations,
#' including `Lab.ID`, `variable_id`, and `score` columns.
#'
#' @return A filtered version of the annotation table with reduced redundancy.
#'
#' @details
#' The function follows these rules:
#' - If a compound has an annotation with a score > 80, it removes annotations with scores ≤ 20.
#' - If a peak has an annotation with a score > 80, it removes other annotations with scores ≤ 20.
#' - Iterates until redundancy stabilizes, ensuring that lower-quality annotations are removed.
#'
#' @export
#'
remove_redundancy <-
  function(annotation_table) {
    ## Track redundancy change across iterations
    redundancy_diff <- c(-1, -1)
    
    while (any(redundancy_diff < 0)) {
      # Calculate redundancy before filtering
      before_redundancy <- calculate_redundance(annotation_table = annotation_table)
      
      # Filter compound annotations: Remove low-scoring annotations if high-scoring ones exist
      annotation_table <- annotation_table %>%
        dplyr::group_by(Lab.ID) %>%
        dplyr::filter(if (any(score > 80)) {
          score > 20
        } else {
          score > 0
        }) %>%
        dplyr::ungroup()
      
      # Filter peak annotations: Remove lower-scoring annotations when a high-scoring annotation exists
      annotation_table <- annotation_table %>%
        dplyr::group_by(variable_id) %>%
        dplyr::filter(if (any(score > 80)) {
          score > 20
        } else {
          score > 0
        }) %>%
        dplyr::ungroup()
      
      # Calculate redundancy after filtering
      after_redundancy <- calculate_redundance(annotation_table = annotation_table)
      
      # Compute change in redundancy
      redundancy_diff <- after_redundancy - before_redundancy
    }
    
    return(annotation_table)
  }




#' Annotate Isotopes for a Given Molecular Formula
#'
#' This function identifies isotopic peaks for a given molecular formula and adduct,
#' matching them with observed peaks based on mass-to-charge ratio (m/z) and retention time (RT).
#'
#' @param formula Character. The molecular formula of the compound (e.g., `"C9H14N2O12P2"`).
#' @param adduct Character. The adduct ion type (e.g., `"M-H"`).
#' @param charge Numeric. The charge state of the ion. Default is `1`.
#' @param mz Numeric. The accurate mass-to-charge ratio (m/z) of the monoisotopic peak.
#' @param rt Numeric. The retention time (RT) of the compound in seconds.
#' @param peak.mz Numeric vector. Observed m/z values of peaks in the dataset.
#' @param peak.rt Numeric vector. Observed retention times (RT) of peaks in the dataset.
#' @param rt.tol Numeric. The retention time tolerance (in seconds) for matching peaks. Default is `5`.
#' @param mz.tol Numeric. The m/z tolerance (in ppm) for matching peaks. Default is `15`.
#' @param max.isotope Numeric. The maximum number of isotopes to consider. Default is `4`.
#'
#' @return A data frame containing:
#'   \item{peakIndex}{Index of the matched peak in the dataset.}
#'   \item{mzError.ppm}{Mass error (in ppm) between the observed and theoretical isotope peak.}
#'   \item{isotopes}{Label of the identified isotope (e.g., `[M+1]`, `[M+2]`).}
#'   \item{rtError.s}{Retention time error (in seconds).}
#'
#' If no isotopic peaks are matched within the given tolerances, the function returns `NULL`.
#'
#' @details
#' The function calculates the theoretical isotope distribution using `Rdisop::getMolecule()`
#' and then matches observed peaks based on m/z and RT tolerances. The function ensures that
#' at least the `[M+1]` isotope is present for a valid match.
#'
#' @importFrom masstools sum_formula
#' @importFrom Rdisop getMolecule getIsotope
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @export

annotate_isotope <-
  function(formula = "C9H14N2O12P2",
           adduct = "M-H",
           charge = 1,
           mz = 402.9998,
           rt = 823.1462,
           ## peak information
           peak.mz,
           peak.rt,
           ## other parameters
           rt.tol = 5,
           mz.tol = 15,
           max.isotope = 4) {
    formula1 <- masstools::sum_formula(formula = formula, adduct = adduct)
    ###should be fix latter
    if (is.na(formula1)) {
      formula1 <- formula
    }
    
    molecule <- Rdisop::getMolecule(formula = formula1, # z = charge,
                                    maxisotopes = max.isotope + 1)
    isotopes <- t(Rdisop::getIsotope(molecule = molecule))
    rownames(isotopes) <-
      c("[M]", paste("[M", "+", c(1:(nrow(
        isotopes
      ) - 1)), "]", sep = ""))
    isotopes <- data.frame(isotopes, rownames(isotopes), stringsAsFactors = FALSE)
    colnames(isotopes) <- c("mz", "intensity", "isotope")
    accurate.mz <- mz
    
    isotopes <- isotopes[-1, , drop = FALSE]
    
    ###rt filtering
    rt.error <- abs(rt - peak.rt)
    index1 <- which(rt.error <= rt.tol)
    if (length(index1) == 0)
      return(NULL)
    
    iso.info <-
      purrr::map(
        as.data.frame(t(isotopes)),
        .f = function(x) {
          temp.mz <- as.numeric(x[1])
          ## calculate error
          peak.mz.error <-
            abs(temp.mz - peak.mz) * 10^6 / ifelse(temp.mz >= 400, temp.mz, 400)
          ##has peak matched the mz tolerance
          idx <- which(peak.mz.error <= mz.tol & rt.error < rt.tol)
          ###idx=0, no peaks matched
          if (length(idx) == 0) {
            # return(c(NA, NA, NA))
            return(NULL)
          }
          ## more than one matched, see the intensity ratio error
          if (length(idx) > 0) {
            idx <- idx[which.min(peak.mz.error[idx])]
            return(c(idx, peak.mz.error[idx], x[3]))
          }
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    if (nrow(iso.info) == 0) {
      return(NULL)
    }
    
    if (!"[M+1]" %in% iso.info$V3) {
      return(NULL)
    }
    
    colnames(iso.info) <- c("peakIndex", "mzError.ppm", "isotopes")
    if (!"[M+2]" %in% iso.info$isotopes) {
      iso.info <- iso.info[1, , drop = FALSE]
    }
    
    iso.info$`rtError.s` <-
      rt.error[as.numeric(iso.info$peakIndex)]
    
    iso.info <-
      iso.info %>%
      dplyr::mutate(
        peakIndex = as.numeric(peakIndex),
        mzError.ppm = as.numeric(mzError.ppm),
        rtError.s = as.numeric(rtError.s)
      )
    
    if (!"[M+2]" %in% iso.info$isotopes) {
      iso.info <- iso.info[1, ]
    }
    
    # rm(list = c("peak.mz", "peak.rt", "peak.int", "cor"))
    # gc()
    return(iso.info)
  }



#' Plot a Metabolic Module from FPA Results
#'
#' This function visualizes a selected metabolic module from the results of a
#' Functional Pathway Analysis (FPA), displaying interactions between metabolites
#' and features in the metabolic network.
#'
#' @param fpa_result List. The output of the Functional Pathway Analysis (FPA),
#'   containing `dysregulated_metabolic_module`, `dysregulated_metabolic_network`,
#'   and `annotation_table`.
#' @param feature_table_marker Data frame. The feature table containing metabolic
#'   markers. Default is `feature_table_marker`.
#' @param include_feature Logical. Whether to include detected metabolic features
#'   in the plot. Default is `FALSE`.
#' @param include_hidden_metabolites Logical. Whether to include hidden metabolites
#'   in the plot. Default is `FALSE`.
#' @param add_compound_name Logical. Whether to add compound names as labels in
#'   the visualization. Default is `TRUE`.
#' @param metabolic_module_index Numeric. The index of the metabolic module to
#'   visualize. Default is `1`.
#' @param layout The layout of the network, such as `kk` or `fr`.
#' @param add_pathways Add pathways beside of the network or not. Default is `FALSE`.
#'
#' @return A `ggplot2` object representing the metabolic module network.
#'
#' @details
#' The function extracts and filters a metabolic subnetwork based on the selected
#' module index. It then uses `ggraph` to visualize the metabolic interactions
#' while allowing customization, such as including metabolic features or hidden
#' metabolites.
#'
#' The network layout follows the Fruchterman-Reingold force-directed algorithm
#' (`layout = "fr"`). Nodes represent metabolites and features, while edges
#' represent different types of metabolic interactions.
#'
#' @importFrom stringr str_split
#' @importFrom dplyr filter pull mutate
#' @importFrom tidygraph activate centrality_degree
#' @importFrom ggraph ggraph geom_edge_arc geom_node_point geom_node_text theme_graph
#' @importFrom ggraph scale_edge_color_manual
#' @export


plot_metabolic_module_fpa <-
  function(fpa_result,
           feature_table_marker = feature_table_marker,
           include_feature = FALSE,
           include_hidden_metabolites = FALSE,
           add_compound_name = TRUE,
           metabolic_module_index = 1,
           layout = "fr",
           add_pathways = FALSE) {
    dysregulated_metabolic_module <-
      fpa_result$dysregulated_metabolic_module
    
    dysregulated_metabolic_network <-
      fpa_result$dysregulated_metabolic_network
    
    annotation_table <-
      fpa_result$annotation_table
    
    remain_compound <-
      dysregulated_metabolic_module$Total_metabolite_id[metabolic_module_index] %>%
      stringr::str_split(pattern = "\\{\\}") %>%
      unlist() %>%
      unique()
    
    remain_feature <-
      annotation_table %>%
      dplyr::filter(KEGG.ID %in% remain_compound) %>%
      dplyr::pull(variable_id) %>%
      unique()
    
    remain_feature <-
      remain_feature[remain_feature %in% igraph::V(dysregulated_metabolic_network)$name]
    
    temp_graph <-
      dysregulated_metabolic_network %>%
      dplyr::filter(name %in% c(remain_compound, remain_feature)) %>%
      tidygraph::activate("nodes") %>%
      dplyr::mutate(degree = tidygraph::centrality_degree()) %>%
      tidygraph::activate("nodes") %>%
      dplyr::filter(degree > 0)
    
    if (!include_feature) {
      temp_graph <-
        temp_graph %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(node_class != "Feature")
    }
    
    if (!include_hidden_metabolites) {
      temp_graph <-
        temp_graph %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(node_class != "Hidden")
    }
    
    plot <-
      temp_graph %>%
      ggraph::ggraph(layout = layout) +
      ggraph::geom_edge_arc(strength = 0.1, aes(color = edge_class)) +
      ggraph::scale_edge_color_manual(values = c(
        "feature_metabolite" =  "#eaeaea",
        "metabolic_reaction" = "#6e8da9"
      )) +
      ggraph::geom_node_point(aes(
        color = node_class,
        shape = node_class,
        size = degree
      )) +
      scale_color_manual(values = c(
        "Feature" = "#ffca3a",
        "Detected" = "#ff595e",
        "Hidden" = "#6e8da9"
      )) +
      ggraph::theme_graph()
    
    if (add_compound_name) {
      plot <- plot +
        ggraph::geom_node_text(aes(label = Compound_name), repel = TRUE)
    }
    
    Total_metabolite_number <-
      fpa_result$dysregulated_metabolic_module$Total_metabolite_number[metabolic_module_index]
    Impact <-
      fpa_result$dysregulated_metabolic_module$Impact[metabolic_module_index]
    Dominant_edge_rate <-
      fpa_result$dysregulated_metabolic_module$Dominant_edge_rate[metabolic_module_index]
    p_value <-
      fpa_result$dysregulated_metabolic_module$p_value[metabolic_module_index]
    Detected_metabolite_number <-
      fpa_result$dysregulated_metabolic_module$Detected_metabolite_number[metabolic_module_index]
    
    plot <-
      plot +
      labs(
        title = fpa_result$dysregulated_metabolic_module$Name[metabolic_module_index],
        subtitle = paste0(
          "Total Metabolite Number: ",
          Total_metabolite_number,
          "\n",
          "Detected Metabolite Number: ",
          Detected_metabolite_number,
          "\n",
          "Impact: ",
          Impact,
          "\n",
          "Dominant Edge Rate: ",
          round(Dominant_edge_rate, 3),
          "\n",
          "P-value: ",
          p_value
        )
      )
    
    if (add_pathways) {
      pathway <-
        fpa_result$enriched_pathways_list[[metabolic_module_index]]
      if (is.null(pathway)) {
        plot_pathway <-
          ggplot() +
          theme_bw() +
          labs(title = "No enriched pathways")
      } else{
        plot_pathway <-
          enrich_scatter_plot(object = pathway, label = TRUE)
      }
      
      plot <-
        plot + patchwork::wrap_plots(plot_pathway, nrow = 1)
      
    }
    
    plot
  }

#' Plot the Dysregulated Metabolic Network from FPA Results
#'
#' This function visualizes the dysregulated metabolic network obtained from Functional Pathway Analysis (FPA),
#' highlighting interactions between metabolites and features in the network.
#'
#' @param fpa_result List. The output of the Functional Pathway Analysis (FPA),
#'   containing `dysregulated_metabolic_module`, `dysregulated_metabolic_network`,
#'   and `annotation_table`.
#' @param feature_table_marker Data frame. The feature table containing metabolic markers.
#' @param include_feature Logical. Whether to include detected metabolic features in the plot.
#'   Default is `FALSE`.
#' @param node_color_by_module Logical. When doesn't include features (include_feature = FALSE),
#' whether set the colors of the nodes by metabolics modules.
#' @param include_hidden_metabolites Logical. Whether to include hidden metabolites in the plot.
#'   Default is `FALSE`.
#' @param add_compound_name Logical. Whether to add compound names as labels in
#'   the visualization. Default is `TRUE`.
#' @param layout The layout of the network, such as `kk` or `fr`.
#'
#' @return A `ggplot2` object representing the metabolic network.
#'
#' @details
#' The function processes the metabolic network by filtering and selecting relevant nodes and edges.
#' It allows the inclusion of metabolic features and hidden metabolites based on user preference.
#' The network is visualized using the Fruchterman-Reingold force-directed algorithm (`layout = "fr"`)
#' to clearly illustrate metabolic interactions.
#'
#' @importFrom tidygraph activate centrality_degree
#' @importFrom dplyr mutate filter
#' @importFrom ggraph ggraph geom_edge_arc geom_node_point geom_node_text theme_graph
#' @importFrom ggplot2 scale_color_manual
#' @importFrom grDevices colorRampPalette
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assume `fpa_result` is a precomputed list from an FPA analysis
#'   plot_metabolic_network_fpa(
#'     fpa_result = fpa_result,
#'     feature_table_marker = feature_table_marker,
#'     include_feature = TRUE,
#'     include_hidden_metabolites = FALSE,
#'     add_compound_name = TRUE
#'   )
#' }


plot_metabolic_network_fpa <-
  function(fpa_result,
           feature_table_marker,
           include_feature = FALSE,
           node_color_by_module = FALSE,
           include_hidden_metabolites = FALSE,
           add_compound_name = TRUE,
           layout = "fr") {
    dysregulated_metabolic_module <-
      fpa_result$dysregulated_metabolic_module
    
    dysregulated_metabolic_network <-
      fpa_result$dysregulated_metabolic_network
    
    annotation_table <-
      fpa_result$annotation_table
    
    temp_graph <-
      dysregulated_metabolic_network %>%
      tidygraph::activate("nodes") %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())
    
    if (!include_feature) {
      temp_graph <-
        temp_graph %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(node_class != "Feature")
      if (node_color_by_module) {
        temp_data <-
          seq_len(nrow(dysregulated_metabolic_module)) %>%
          purrr::map(function(i) {
            data.frame(
              module = dysregulated_metabolic_module$Name[i],
              KEGG.ID = stringr::str_split(
                dysregulated_metabolic_module$Total_metabolite_id[i],
                "\\{\\}"
              )[[1]]
            ) %>%
              dplyr::distinct(KEGG.ID, .keep_all = TRUE)
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::mutate(module = factor(module, levels = stringr::str_sort(unique(module), numeric = TRUE)))
        
        temp_graph <-
          temp_graph %>%
          tidygraph::activate("nodes") %>%
          dplyr::left_join(temp_data, by = c("name" = "KEGG.ID"))
        
      }
    }
    
    if (!include_hidden_metabolites) {
      temp_graph <-
        temp_graph %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(node_class != "Hidden")
    }
    
    if (!include_feature & node_color_by_module) {
      plot <-
        temp_graph %>%
        ggraph::ggraph(layout = layout) +
        ggraph::geom_edge_arc(strength = 0.1, aes(color = edge_class)) +
        ggraph::scale_edge_color_manual(values = c(
          "feature_metabolite" =  "#eaeaea",
          "metabolic_reaction" = "#6e8da9"
        )) +
        ggraph::geom_node_point(aes(
          color = module,
          shape = node_class,
          size = degree
        )) +
        scale_color_manual(values =
                             colorRampPalette(ggsci::pal_aaas("default")(10))(length(unique(
                               as.character(
                                 tidygraph::as_tibble(temp_graph, active = "nodes") %>% pull(module)
                               )
                             )))) +
        ggraph::theme_graph()
    } else{
      plot <-
        temp_graph %>%
        ggraph::ggraph(layout = layout) +
        ggraph::geom_edge_arc(strength = 0.1, aes(color = edge_class)) +
        ggraph::scale_edge_color_manual(values = c(
          "feature_metabolite" =  "#eaeaea",
          "metabolic_reaction" = "#6e8da9"
        )) +
        ggraph::geom_node_point(aes(
          color = node_class,
          shape = node_class,
          size = degree
        )) +
        scale_color_manual(values = c(
          "Feature" = "#ffca3a",
          "Detected" = "#ff595e",
          "Hidden" = "#6e8da9"
        )) +
        ggraph::theme_graph()
    }
    
    
    if (add_compound_name) {
      plot <- plot +
        ggraph::geom_node_text(aes(label = Compound_name), repel = TRUE)
    }
    
    plot
  }




#' Calculate Node Centrality in a Graph
#'
#' This function calculates the centrality of nodes in a given network graph based on either **degree centrality** or **betweenness centrality**.
#'
#' @param graph An `igraph` object representing the network.
#' @param type A character string specifying the type of centrality to compute. Options are `"degree"` (default) or `"betweenness"`.
#'
#' @return A named numeric vector where names represent node IDs and values represent centrality scores.
#'
#' @details
#' - **Degree centrality** measures the number of direct connections a node has.
#' - **Betweenness centrality** measures the number of shortest paths passing through a node, indicating its role as a bridge in the network.
#' - The function uses `igraph::degree()` for degree centrality and a custom approach for betweenness centrality using `igraph::shortest_paths()`.
#'
#' @importFrom igraph V degree shortest_paths
#'
#' @examples
#' \dontrun{
#' library(igraph)
#'
#' # Create a sample graph
#' g <- make_ring(10)
#' V(g)$name <- paste0("Node", 1:10)
#'
#' # Calculate degree centrality
#' degree_centrality <- calculate_centrality(g, type = "degree")
#' print(degree_centrality)
#'
#' # Calculate betweenness centrality
#' betweenness_centrality <- calculate_centrality(g, type = "betweenness")
#' print(betweenness_centrality)
#' }
#'
#' @export

calculate_centrality <-
  function(graph, type = c("degree", "betweenness")) {
    type = match.arg(type)
    node <- igraph::V(graph)$name
    if (type == "degree") {
      node.degree <- igraph::degree(graph = graph, v = node)
      return(node.degree)
    }
    
    node.betweenness <- lapply(node, function(x) {
      i <- match(x, node)
      if (i == length(node))
        return(NULL)
      from <- x
      to <- node[(i + 1):length(node)]
      temp  <- igraph::shortest_paths(graph = graph,
                                      from = from,
                                      to = to)[[1]]
      temp <- lapply(temp, function(x) {
        if (length(x) <= 2)
          return(NULL)
        names(x)[-c(1, length(x))]
      })
      unlist(temp)
    })
    
    node.betweenness <- unlist(node.betweenness)
    node.betweenness <- table(node.betweenness)
    node0 <- setdiff(node, names(node.betweenness))
    if (length(node0) == 0)
      return(node.betweenness)
    add <- rep(0, length(node0))
    names(add) <- node0
    node.betweenness <- c(node.betweenness, add)
    node.betweenness <- node.betweenness[match(node, names(node.betweenness))]
    return(node.betweenness)
  }




#' Calculate Activity Score for Metabolic Modules
#'
#' This function calculates an activity score for each metabolic module based on its detected and hidden metabolites,
#' as well as its structure within a sub-metabolic network. The calculation incorporates modularity and degree-based
#' weighting to assess the significance of the module.
#'
#' @param metabolic_modules A list of metabolic modules, where each module is represented as a vector of metabolite IDs.
#' @param detected_metabolites A character vector of detected metabolite IDs.
#' @param hidden_metabolites A character vector of hidden metabolite IDs.
#' @param sub_metabolic_network An `igraph` object representing the sub-metabolic network.
#' @param threads An integer specifying the number of threads to use for parallel computation. Default is 3.
#'
#' @return A numeric vector where each value represents the activity score of a corresponding metabolic module.
#'
#' @details
#' The activity score is computed by assessing the modularity of each metabolic module within the sub-metabolic network.
#' The function extracts the subgraph corresponding to each module and computes a modularity-based score adjusted for the
#' total number of metabolites in the dataset.
#'
#' The computation is performed in parallel using `BiocParallel`, selecting either `MulticoreParam` (for Unix-based systems)
#' or `SnowParam` (for Windows) to handle multi-threaded execution efficiently.
#'
#' @import igraph
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam
#' @importFrom masstools get_os
#'
#' @examples
#' \dontrun{
#' library(igraph)
#'
#' # Example metabolic network
#' g <- make_ring(10)
#' V(g)$name <- paste0("M", 1:10)
#'
#' # Example input data
#' metabolic_modules <- list(
#'   c("M1", "M2", "M3"),
#'   c("M4", "M5", "M6", "M7")
#' )
#' detected_metabolites <- c("M1", "M3", "M5", "M7")
#' hidden_metabolites <- c("M2", "M6")
#'
#' # Compute activity scores
#' scores <- calculate_activity_score(
#'   metabolic_modules = metabolic_modules,
#'   detected_metabolites = detected_metabolites,
#'   hidden_metabolites = hidden_metabolites,
#'   sub_metabolic_network = g,
#'   threads = 2
#' )
#'
#' print(scores)
#' }
#'
#' @export

calculate_activity_socre <-
  function(metabolic_modules,
           detected_metabolites,
           hidden_metabolites,
           sub_metabolic_network,
           threads = 3) {
    temp_fun_activity_scores <-
      function(x,
               detected_metabolites,
               hidden_metabolites,
               total_metabolite_number,
               sub_metabolic_network = sub_metabolic_network) {
        requireNamespace("igraph", quietly = TRUE)
        # options(warn = -1)
        
        temp_graph <- igraph::subgraph(graph = sub_metabolic_network, v = x)
        detected_number <- sum(x %in% detected_metabolites)
        hidden_number <- sum(x %in% hidden_metabolites)
        
        m <- length(igraph::E(sub_metabolic_network))
        e <- length(igraph::E(temp_graph))
        
        node_name <- names(igraph::V(temp_graph))
        
        value <- sapply(node_name, function(node) {
          temp_value <- sapply(node_name, function(x) {
            temp1 <- igraph::degree(graph = sub_metabolic_network, v = x)
            temp2 <- igraph::degree(graph = sub_metabolic_network, v = node)
            temp1 * temp2 / (4 * m^2)
          })
          sum(temp_value)
        })
        
        modularity <- e / m - sum(value)
        # modularity
        q <- modularity * sqrt(total_metabolite_number / length(x))
        a <- detected_number * q / length(x)
        a
      }
    
    if (masstools::get_os() == "windows") {
      bpparam = BiocParallel::SnowParam(workers = threads, progressbar = TRUE)
    } else{
      bpparam = BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)
    }
    
    activity_score <-
      BiocParallel::bplapply(
        metabolic_modules,
        temp_fun_activity_scores,
        BPPARAM = bpparam,
        detected_metabolites = detected_metabolites,
        hidden_metabolites = hidden_metabolites,
        total_metabolite_number = length(detected_metabolites),
        sub_metabolic_network = sub_metabolic_network
      )
    
    activity_score <- unlist(activity_score)
    return(activity_score)
    
  }

















# Function to compute module quality metrics
calculate_module_quality <-
  function(graph, nodes) {
    temp_graph <-
      igraph::induced_subgraph(graph, vids = nodes)
    internal_edges <- igraph::gsize(temp_graph)
    external_edges <- sum(igraph::degree(graph, v = nodes)) - (2 * internal_edges)
    
    # Conductance
    conductance <- ifelse((internal_edges + external_edges) == 0,
                          0,
                          external_edges / (internal_edges + external_edges)
    )
    
    # Internal Density
    num_nodes <- length(nodes)
    internal_density <-
      ifelse(num_nodes < 2, 0, (2 * internal_edges) / (num_nodes * (num_nodes - 1)))
    
    # Edge Density Ratio (EDR)
    edge_density_ratio <-
      ifelse(external_edges == 0, 1, internal_edges / external_edges)
    
    # Return metrics
    return(
      data.frame(
        module_size = num_nodes,
        internal_edges = internal_edges,
        external_edges = external_edges,
        conductance = conductance,
        internal_density = internal_density,
        edge_density_ratio = edge_density_ratio
      )
    )
  }












generate_null_conductance_distribution <-
  function(annotation_table_all,
           feature_table_marker,
           metabolite_database,
           metabolic_network,
           permutation_times = 5,
           threads = 8,
           adduct.table = NULL,
           column = c("rp", "hilic"),
           ms1.match.ppm = 25,
           mz.ppm.thr = 400,
           include_hidden_metabolites = FALSE) {
    column <- match.arg(column)  # Ensure valid column type
    
    # Initialize list to store activity scores for each permutation
    null_conductance <- vector(mode = "list", length = permutation_times)
    marker_number <- nrow(feature_table_marker)  # Number of markers
    
    for (i in seq_len(permutation_times)) {
      cat(i, " ")
      
      # Randomly sample markers for permutation
      if (nrow(feature_table_marker) / length(unique(annotation_table_all$variable_id)) < 0.2) {
        temp_marker_variable_id <- sample(
          x = setdiff(
            unique(annotation_table_all$variable_id),
            feature_table_marker$variable_id
          ),
          size = marker_number,
          replace = FALSE
        )
      } else{
        temp_marker_variable_id <- sample(
          x = unique(annotation_table_all$variable_id),
          size = marker_number,
          replace = FALSE
        )
      }
      
      
      # Extract annotation information for sampled markers
      temp_annotation_table_marker <- annotation_table_all %>%
        dplyr::filter(variable_id %in% temp_marker_variable_id)
      
      # Extract detected metabolites from annotation table
      temp_detected_metabolites <- unique(temp_annotation_table_marker$KEGG.ID)
      
      # Identify hidden metabolites connected to detected metabolites
      if (include_hidden_metabolites) {
        temp_hidden_metabolites <- get_hidden_metabolites(
          metabolic_network = metabolic_network,
          detected_metabolites = temp_detected_metabolites,
          threads = threads,
          max.reaction = 3
        ) %>% unlist() %>% unique()
      } else {
        temp_hidden_metabolites <- NULL
      }
      
      temp_total_metabolites <-
        c(temp_detected_metabolites, temp_hidden_metabolites)
      
      # Extract subnetwork containing detected and hidden metabolites
      temp_sub_metabolic_network <- metabolic_network %>%
        tidygraph::activate("nodes") %>%
        dplyr::filter(name %in% temp_total_metabolites) %>%
        dplyr::mutate(node_class = ifelse(name %in% temp_detected_metabolites, "Detected", "Hidden")) %>%
        tidygraph::mutate(degree2 = tidygraph::centrality_degree())
      
      # Identify metabolic modules in the subnetwork
      temp_metabolic_modules <- identify_metabolic_modules(
        sub_metabolic_network = temp_sub_metabolic_network,
        detected_metabolites = temp_detected_metabolites,
        hidden_metabolites = temp_hidden_metabolites
      )
      
      # Calculate activity scores for the identified metabolic modules
      temp_conductance <-
        temp_metabolic_modules %>%
        purrr::map(function(x) {
          calculate_module_quality(graph = temp_sub_metabolic_network, nodes = x)$conductance
        }) %>%
        unlist()
      
      # Store the computed activity score in the list
      null_conductance[[i]] <- temp_conductance
    }
    
    return(null_conductance)
  }







check_feature_table_all <-
  function(feature_table_all) {
    ####feature_table_all
    # 1. Check if it's a data frame
    if (!is.data.frame(feature_table_all)) {
      stop("Error: feature_table_all must be a data frame.")
    }
    
    # 2. Check column names
    expected_colnames <- c("variable_id", "mz", "rt", "polarity")
    if (!identical(colnames(feature_table_all), expected_colnames)) {
      stop("Error: Column names must be exactly: 'variable_id', 'mz', 'rt', 'polarity'.")
    }
    
    # 3. Check if variable_id is unique
    if (anyDuplicated(feature_table_all$variable_id) > 0) {
      stop("Error: 'variable_id' column must have unique values.")
    }
    
    # 4. Check if mz is numeric and has no NA
    if (!is.numeric(feature_table_all$mz) ||
        any(is.na(feature_table_all$mz))) {
      stop("Error: 'mz' column must be numeric and contain no missing values.")
    }
    
    # 5. Check if rt is numeric and has no NA
    if (!is.numeric(feature_table_all$rt) ||
        any(is.na(feature_table_all$rt))) {
      stop("Error: 'rt' column must be numeric and contain no missing values.")
    }
    
    # 6. Check if rt is in seconds (not in minutes)
    max_rt <- max(feature_table_all$rt, na.rm = TRUE)  # Get the max retention time
    if (max_rt < 100) {
      # If max RT is too small, it might be in minutes
      stop("Error: 'rt' values seem to be in minutes, but they should be in seconds.")
    }
    
    # 7. Check polarity column
    valid_polarity <- c("positive", "negative")
    if (any(is.na(feature_table_all$polarity))) {
      stop("Error: 'polarity' column must not contain NA values.")
    }
    if (any(grepl("\\s", feature_table_all$polarity))) {
      stop("Error: 'polarity' column must not contain spaces.")
    }
    if (!all(feature_table_all$polarity %in% valid_polarity)) {
      stop("Error: 'polarity' must only contain 'positive' or 'negative'.")
    }
    
    # If all checks pass
    message("feature_table_all is correct.\n")
    
  }


check_feature_table_marker <-
  function(feature_table_marker) {
    # 1. Check if it's a data frame
    if (!is.data.frame(feature_table_marker)) {
      stop("Error: feature_table_marker must be a data frame.")
    }
    
    # 2. Check column names
    expected_colnames <- c("variable_id", "mz", "rt", "degree", "p_value", "polarity")
    if (!identical(colnames(feature_table_marker), expected_colnames)) {
      stop(
        "Error: Column names must be exactly: 'variable_id', 'mz', 'rt', 'degree', 'p_value', 'polarity'."
      )
    }
    
    # 3. Check if variable_id is unique
    if (anyDuplicated(feature_table_marker$variable_id) > 0) {
      stop("Error: 'variable_id' column must have unique values.")
    }
    
    # 4. Check if mz is numeric and has no NA
    if (!is.numeric(feature_table_marker$mz) ||
        any(is.na(feature_table_marker$mz))) {
      stop("Error: 'mz' column must be numeric and contain no missing values.")
    }
    
    # 5. Check if rt is numeric, has no NA, and represents time in seconds
    if (!is.numeric(feature_table_marker$rt) ||
        any(is.na(feature_table_marker$rt))) {
      stop("Error: 'rt' column must be numeric and contain no missing values.")
    }
    
    # 6. Ensure rt is in seconds (not minutes)
    max_rt <- max(feature_table_marker$rt, na.rm = TRUE)  # Get the max retention time
    if (max_rt < 100) {
      # If max RT is too small, it might be in minutes
      stop("Error: 'rt' values seem to be in minutes, but they should be in seconds.")
    }
    
    # 7. Check if degree is numeric and has no NA
    if (!is.numeric(feature_table_marker$degree) ||
        any(is.na(feature_table_marker$degree))) {
      stop("Error: 'degree' column must be numeric and contain no missing values.")
    }
    
    # 8. Check if p_value is numeric and has no NA
    if (!is.numeric(feature_table_marker$p_value) ||
        any(is.na(feature_table_marker$p_value))) {
      stop("Error: 'p_value' column must be numeric and contain no missing values.")
    }
    
    # 9. Check polarity column
    valid_polarity <- c("positive", "negative")
    if (any(is.na(feature_table_marker$polarity))) {
      stop("Error: 'polarity' column must not contain NA values.")
    }
    if (any(grepl("\\s", feature_table_marker$polarity))) {
      stop("Error: 'polarity' column must not contain spaces.")
    }
    if (!all(feature_table_marker$polarity %in% valid_polarity)) {
      stop("Error: 'polarity' must only contain 'positive' or 'negative'.")
    }
    
    # If all checks pass
    message("feature_table_marker is correct.")
  }

# Example usage:
# check_feature_table_marker(feature_table_marker)