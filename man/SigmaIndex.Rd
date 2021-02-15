\docType{class}
\name{SigmaIndex-class}
\alias{SigmaIndex}
\alias{convertFromDSD}
\alias{addPopulation}
\alias{addDataPoints}
\alias{addDataPointsInc}
\alias{queryDataPoints}
\alias{getTotalPopulations}
\alias{print}
\alias{getPopulations}
\alias{removePopulation}
\alias{resetStatistics}
\alias{getHistogram}
\alias{getStatistics}
\title{Sigma-Index: a statistical indexing organizing structure}
\usage{
SigmaIndex(theta=3, neighborhood=9, precision_switch=TRUE) 
convertFromDSD(x, total_elements=500, theta=3, neighborhood=9, precision_switch=TRUE) 
}
\arguments{
\item{theta}{(numeric) - A threshold that defines the bound of the statistical population. Used in statistical classification and inference.}
\item{neighborhood}{(numeric) - A neighborhood statistical threshold, used to form the sigma-index DAG. Must be >\code{theta}.}
\item{precision_switch}{(logical) - Keep the precision over the speed.}
\item{x}{(DSD_Gaussians) - Gaussians data stream definition from the \pkg{stream} package. Gaussian definitions and outliers are used to covert to the Sigma-Index DAG. The same DSD can be used later to test the created Sigma-Index.}
\item{total_elements}{(integer) - When converting a Sigma-Index DAG from DSD, we need to supply a total number of data points, to calculate population elements and numbers, since the Sigma-Index is probability structured.}
}
\section{Methods}{
\describe{
All methods here are detailed in the stream package.
\item{\code{addPopulation(x, id, mean, covariance, elements=1, ...)}}{
  Adds a new population to the SigmaIndex object \code{x}. Population is characterized by its identifier (\code{id} - (string)), a centroid (\code{mean} - (numeric vector)), covariance matrix (\code{covariance} - (matrix)), and a number of elements (\code{elements} - (integer)). The covariance matrix must be non-singular and invertible.}
\item{\code{addDataPoints(x, id, data_points, neighborhood, ...)}}{
  Updates a statistical population with one or more data points in the SigmaIndex object \code{x}. If adding one value, a caller must supply the population identifier (\code{id} - (string)) and a new data point (\code{data_points} - (numeric vector)). When supplying multiple data points, \code{id} must be a string vector containing population identifiers for all added data points, and \code{data_points} must be a \emph{data frame} containing all data points.
  }
\item{\code{addDataPointsInc(x, data_points, query_results, ...)}}{
  Incrementally updates a statistical population with one or more data points in the SigmaIndex object \code{x}. Prior to calling this method, a query must be done on the same set of data points. Results of this query are then used to incrementally update the SigmaIndex DAG. See the example.
  }
\item{\code{queryDataPoints(x, data_points, ...)}}{
Performs a query on the SigmaIndex object \code{x}. \code{data_points} parameter can be one value or a data frame comprising a set of values. Returns a list having the following members:
  \itemize{
    \item \code{classified} - a set of population identifiers for which a tested data point is closer (<=) than the \code{theta} threshold to the population centroid
    \item \code{neighborhood} - a set of popluations for which each tested data point is between \code{theta} and \code{neighborhood} thresholds. During this query call, SigmaIndex object is collecting computation cost reduction statistics.}
  }
\item{\code{getTotalPopulations(x, ...)}}{
Returns the number of populations in the SigmaIndex object \code{x}.}
\item{\code{print(x, ...)}}{
Prints details and structure of the SigmaIndex object \code{x}.}
\item{\code{getPopulations(x, ...)}}{
Returns details of all populations in the SigmaIndex object \code{x}. Returns a named list having the following members:
  \itemize{
    \item the name of the element is the population identifier
    \item \code{mean} (numeric vector) - the population centroid
    \item \code{icovariance} (matrix) - the inverted covariance matrix of the population
    \item \code{elements} (integer) - the number of the population elements}
  }
\item{\code{removePopulation(x, id, ...)}}{
Removes the population having identifier \code{id} from the SigmaIndex object \code{x}.}
\item{\code{resetStatistics(x, ...)}}{
Resets the statistical counters of the SigmaIndex object \code{x}.}
\item{\code{getHistogram(x, ...)}}{
Returns a numeric vector comprising 101 integers, which represents a histogram of the SigmaIndex object \code{x}. Such histogram contains the number of population elements related to the specific computation cost reduction.}
\item{\code{getStatistics(x, ...)}}{
Returns the statistical counters for the SigmaIndex object \code{x}. These counters can be resetted by the \code{resetStatistics} call. This method returns a list comprising the following members:
  \itemize{
    \item \code{totalCount} - Total populations searched in all queries
    \item \code{classifiedNodes} - Total number of populations that had queried data points in range <=\code{neighborhood} from their centroids
    \item \code{missedNodes} - Total number of populations that had queried data points in range >\code{neihgborhood} from their centroids
    \item \code{sequentialNodes} - Total number of populations that would be searched in the \emph{sequential scan} approach
    \item \code{computationCostReduction} - The overall computational cost reduction comparing to the \emph{sequential scan} approacj
  }
}
}}
\description{
The Sigma-Index statistical organizing structure.
}
\details{
A Sigma-Index is a statistical organizing structure aiming to improve query and processing times for statistical algorithms, when searching space comprising statistical populations. Those interested in inspecting the Sigma-Index code, must start from the following units:
\itemize{
  \item \code{SigmaIndex.hpp}, \code{SigmaIndex.cpp}, \code{SigmaIndex_Inline.cpp}, \code{SigmaIndexProxy.hpp}, and \code{SigmaIndexProxy.cpp} - C++ units comprising the Sigma-Index implementation suitable for the Statistical Hierachical Clusterer.
  \item \code{SHC.cpp} - C++ unit that uses the Sigma-Index. After construction of the SHC object [1], we create and add a Sigma-Index using the \code{void useSigmaIndex(int neighborhood_mutiplier=2, bool precision_switch=true);} method. Each processing is done through the SHC method \code{shared_ptr<ClassificationResult> process(VectorXd *newElement, bool classifyOnly=false);}, which is the primary place for the Sigma-Index invocations.
}

The Sigma-Index computational cost reduction is covered by three (3) distinct tests:
\enumerate{
\item A generic synthetic test - \link{SHC_TestCase_ClustersAndOutliers_SigmaIndex}
\item A synthetic test covering balancing theorems - \link{SHC_TestCase_ClustersAndOutliers_SigmaIndex_Theorems}
\item A real-life sensor dataset test - \link{SHC_TestCase_Sensors_SigmaIndex}
}
}
\examples{ % must be 150
d <- DSD_Gaussians(k=40,outliers=40,separation_type="Mahalanobis",separation=4,
                   space_limit=c(0,150),variance_limit=4,
                   outlier_options=list(outlier_horizon=200000))
si <- convertFromDSD(d,200000,theta=3.2,neighborhood=12.9)
res <- queryDataPoints(si, get_points(d, 200000))
hist <- getHistogram(si)
stat <- getStatistics(si)

# incremental update
points <- get_points(d, 200)
query_res <- queryDataPoints(si, points)
addDataPointsInc(si, points, query_res)
}
\references{
[1] Krleža D, Vrdoljak B, and Brčić M, Statistical hierarchical clustering algorithm for outlier detection in evolving data streams, \emph{Machine Learning}, Sep. 2020
}