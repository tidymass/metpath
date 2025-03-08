% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_feature_based_pathway_enrichment.R
\name{generate_null_activity_score_distribution}
\alias{generate_null_activity_score_distribution}
\title{Generate Null Distribution of Activity Scores}
\usage{
generate_null_activity_score_distribution(
  annotation_table_all,
  feature_table_marker,
  metabolite_database,
  metabolic_network,
  permutation_times = 20,
  threads = 8,
  adduct.table = NULL,
  column = c("rp", "hilic"),
  ms1.match.ppm = 25,
  mz.ppm.thr = 400,
  include_hidden_metabolites = FALSE
)
}
\arguments{
\item{annotation_table_all}{A data frame containing all annotated metabolites.}

\item{feature_table_marker}{A data frame containing marker features.}

\item{metabolite_database}{A metabolite database object used for annotation.}

\item{metabolic_network}{An `igraph` object representing the metabolic network.}

\item{permutation_times}{Integer; number of permutations to perform. Default is `20`.}

\item{threads}{Integer; number of parallel threads to use. Default is `8`.}

\item{adduct.table}{Optional data frame specifying possible adducts. Default is `NULL`.}

\item{column}{Character; chromatography column type. Options are `"rp"` or `"hilic"`. Default is `"rp"`.}

\item{ms1.match.ppm}{Numeric; MS1 matching tolerance in ppm. Default is `25`.}

\item{mz.ppm.thr}{Numeric; M/Z tolerance threshold in ppm. Default is `400`.}

\item{include_hidden_metabolites}{Logical; whether to include hidden metabolites in the network analysis. Default is `FALSE`.}
}
\value{
A list of numeric vectors representing null activity scores for each permutation.
}
\description{
This function generates a null distribution of activity scores by performing random sampling
of marker metabolites and calculating activity scores through network analysis.
}
