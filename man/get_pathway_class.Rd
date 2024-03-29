% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway_database-methods.R
\name{get_pathway_class}
\alias{get_pathway_class}
\title{get_pathway_class}
\usage{
get_pathway_class(object)
}
\arguments{
\item{object}{pathway_class object.}
}
\description{
Extract the class of pathways.
}
\examples{
library(metpath)
kegg_hsa_pathway = get_kegg_pathway(local = TRUE, organism = "hsa")
kegg_hsa_pathway

get_pathway_class(kegg_hsa_pathway)
}
\author{
Xiaotao Shen
\email{shenxt1990@outlook.com}
}
