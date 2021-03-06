% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureFiltration.R
\name{featureFiltration}
\alias{featureFiltration}
\title{Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes}
\usage{
featureFiltration(
  seurat.obj,
  mito.low = 0.05,
  mito.high = 0.975,
  feature.cut = 0.975
)
}
\arguments{
\item{seurat.obj}{A seurat object.}

\item{mito.low}{Filter cells that have total mitochondrial genes low than this quantile}

\item{mito.high}{Filter cells that have total mitochondrial genes higher than this quantile}

\item{feature.cut}{Filter cells that have total features higher than this quantile}
}
\description{
Filters cells from a Seurat object based upon the amount of features and percentage of mitochondrial genes
}
\examples{
featureFiltration(seurat.obj, mito.low = .5, mito.high = .975, feature.cut = .975)
}
