% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_feature_based_pathway_enrichment.R
\name{calculate_activity_socre}
\alias{calculate_activity_socre}
\title{Calculate Activity Score for Metabolic Modules}
\usage{
calculate_activity_socre(
  metabolic_modules,
  detected_metabolites,
  hidden_metabolites,
  sub_metabolic_network,
  threads = 3
)
}
\arguments{
\item{metabolic_modules}{A list of metabolic modules, where each module is represented as a vector of metabolite IDs.}

\item{detected_metabolites}{A character vector of detected metabolite IDs.}

\item{hidden_metabolites}{A character vector of hidden metabolite IDs.}

\item{sub_metabolic_network}{An `igraph` object representing the sub-metabolic network.}

\item{threads}{An integer specifying the number of threads to use for parallel computation. Default is 3.}
}
\value{
A numeric vector where each value represents the activity score of a corresponding metabolic module.
}
\description{
This function calculates an activity score for each metabolic module based on its detected and hidden metabolites,
as well as its structure within a sub-metabolic network. The calculation incorporates modularity and degree-based
weighting to assess the significance of the module.
}
\details{
The activity score is computed by assessing the modularity of each metabolic module within the sub-metabolic network.
The function extracts the subgraph corresponding to each module and computes a modularity-based score adjusted for the
total number of metabolites in the dataset.

The computation is performed in parallel using `BiocParallel`, selecting either `MulticoreParam` (for Unix-based systems)
or `SnowParam` (for Windows) to handle multi-threaded execution efficiently.
}
\examples{
\dontrun{
library(igraph)

# Example metabolic network
g <- make_ring(10)
V(g)$name <- paste0("M", 1:10)

# Example input data
metabolic_modules <- list(
  c("M1", "M2", "M3"),
  c("M4", "M5", "M6", "M7")
)
detected_metabolites <- c("M1", "M3", "M5", "M7")
hidden_metabolites <- c("M2", "M6")

# Compute activity scores
scores <- calculate_activity_score(
  metabolic_modules = metabolic_modules,
  detected_metabolites = detected_metabolites,
  hidden_metabolites = hidden_metabolites,
  sub_metabolic_network = g,
  threads = 2
)

print(scores)
}

}
