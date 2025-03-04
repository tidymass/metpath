% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway_database-methods.R
\name{[.pathway_database}
\alias{[.pathway_database}
\alias{names.pathway_database}
\alias{length.pathway_database}
\alias{database_info,pathway_database-method}
\alias{pathway_id,pathway_database-method}
\alias{pathway_name,pathway_database-method}
\alias{describtion,pathway_database-method}
\alias{pathway_class,pathway_database-method}
\alias{gene_list,pathway_database-method}
\alias{compound_list,pathway_database-method}
\alias{reference_list,pathway_database-method}
\alias{related_disease,pathway_database-method}
\alias{related_module,pathway_database-method}
\title{[ method}
\usage{
\method{[}{pathway_database}(x, i, ...)

\method{names}{pathway_database}(x)

\method{length}{pathway_database}(x)

\S4method{database_info}{pathway_database}(object)

\S4method{pathway_id}{pathway_database}(object)

\S4method{pathway_name}{pathway_database}(object)

\S4method{describtion}{pathway_database}(object)

\S4method{pathway_class}{pathway_database}(object)

\S4method{gene_list}{pathway_database}(object)

\S4method{compound_list}{pathway_database}(object)

\S4method{reference_list}{pathway_database}(object)

\S4method{related_disease}{pathway_database}(object)

\S4method{related_module}{pathway_database}(object)
}
\arguments{
\item{x}{A `pathway_database` class object.}

\item{i}{i}

\item{...}{other parameters}

\item{object}{A `pathway_database` class object.}
}
\value{
pathway_database

A vector of pathway IDs.

An integer indicating the number of pathways.

A list containing database metadata.

A character vector of pathway IDs.

A character vector of pathway names.

A character vector of pathway descriptions.

A character vector of pathway classes.

A list of gene sets for each pathway.

A list of compound sets for each pathway.

A list of references for each pathway.

A list of related diseases for each pathway.

A list of related modules for each pathway.
}
\description{
Returns pathway IDs from a `pathway_database` object.

Returns the number of pathways in a `pathway_database` object.

Extracts database metadata from a `pathway_database` object.

Extracts pathway IDs from a `pathway_database` object.

Extracts pathway names from a `pathway_database` object.

Extracts pathway descriptions from a `pathway_database` object.

Extracts pathway classification from a `pathway_database` object.

Extracts lists of genes involved in pathways from a `pathway_database` object.

Extracts lists of compounds involved in pathways from a `pathway_database` object.

Extracts literature references for pathways in a `pathway_database` object.

Extracts disease associations from a `pathway_database` object.

Extracts related modules from a `pathway_database` object.
}
\author{
Xiaotao Shen
\email{shenxt1990@outlook.com}

Xiaotao Shen (\email{shenxt1990@outlook.com})
}
