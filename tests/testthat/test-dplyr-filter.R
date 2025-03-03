library(metpath)
library(stringr)
data("kegg_hsa_pathway")
test_that("filter.pathway_database", {
  kegg_hsa_pathway <-
    kegg_hsa_pathway %>% 
    filter(stringr::str_detect(pathway_class, "Human Diseases"))
  testthat::expect_match(object = unlist(kegg_hsa_pathway@pathway_class), "Human Diseases")
})



data("kegg_hsa_pathway")

remain_idx =
  kegg_hsa_pathway@pathway_class %>%
  unlist() %>%
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

pathway_database =
  filter_pathway(object = kegg_hsa_pathway, 
                 remain_idx = remain_idx[1:50])

data("query_id_kegg")

kegg_enrichment =
  enrich_kegg(
    query_id = query_id_kegg,
    query_type = "compound",
    id_type = "KEGG",
    pathway_database = pathway_database,
    p_cutoff = 0.05,
    p_adjust_method = "BH",
    method = "hypergeometric",
    threads = 5
  )

test_that("filter.enrich_result", {
  kegg_enrichment <-
    kegg_enrichment %>% 
    filter(p_value < 0.05)
  testthat::expect_true(all(kegg_enrichment@result$p_value < 0.05))
})
