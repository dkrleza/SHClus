\name{SigmaIndex-SyntheticTests}
\alias{SHC_TestCase_ClustersAndOutliers_SigmaIndex}
\alias{SHC_TestCase_ClustersAndOutliers_SigmaIndex_Theorems}
\title{Sigma-Index Synthetic Testing Procedures}
\usage{
SHC_TestCase_ClustersAndOutliers_SigmaIndex(n=100000)
SHC_TestCase_ClustersAndOutliers_SigmaIndex_Theorems()
}
\arguments{
\item{n}{(integer) - A number of data points generated for the test.}
}
\description{
Sigma-Index is a statistical organizing structure aiming to improve query and processing times for statistical algorithms, when searching space comprising statistical populations. These synthetic tests compare sequential scan approach (iterating through all populations) with Sigma-Index approach. Generated diagrams are placed in \code{./inst} folder.
\itemize{
  \item \code{SHC_TestCase_ClustersAndOutliers_SigmaIndex()} - Generates comparison diagrams for synthetically generated set of clusters and outliers. This is the basic testing of the Sigma-Index improvements and usability.
  \item \code{SHC_TestCase_ClustersAndOutliers_SigmaIndex_Theorems()} - Generates computational cost reduction histograms that experimentally proves theorems about balacing the Sigma-Index DAG.
}
}
\author{
Dalibor Krleža
}
\references{
[1] Krleža D, Vrdoljak B, and Brčić M, Statistical hierarchical clustering algorithm for outlier detection in evolving data streams, \emph{Machine Learning}, Sep. 2020
}