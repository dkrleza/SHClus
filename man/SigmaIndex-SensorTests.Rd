\name{SigmaIndex-SensorTests}
\alias{SHC_TestCase_Sensors_SigmaIndex}
\alias{SHC_TestCase_Sensors_SigmaIndex_Display}
\alias{prepareSensorDataset}
\title{Sigma-Index Sensor Dataset Testing Procedures}
\usage{
SHC_TestCase_Sensors_SigmaIndex(endingRound = 2000)
SHC_TestCase_Sensors_SigmaIndex_Display()
prepareSensorDataset()
}
\arguments{
\item{endingRound}{(integer) - A number of data points processed from the sensor dataset in thousands, i.e., 2000 rounds = 2 million data points.}
}
\description{
Sigma-Index is a statistical organizing structure aiming to improve query and processing times for statistical algorithms, when searching space comprising statistical populations. These sensor dataset tests compare sequential scan approach (iterating through all populations) with Sigma-Index approach. Generated diagrams are placed in the \code{./inst} folder.
\itemize{
  \item \code{SHC_TestCase_Sensors_SigmaIndex()} - Processes the sensor dataset and places the results in the RDS files in the \code{./inst} folder. These RDS files can be further used to generate diagrams.
  \item \code{SHC_TestCase_Sensors_SigmaIndex_Display()} - Generates diagrams from RDS files in the \code{./inst} folder.
  \item \code{prepareSensorDataset()} - A function implicitly called by the \code{SHC_TestCase_Sensors_SigmaIndex()} function to download and prepare the sensor dataset. Since CRAN does not allow packages comprising big files, the sensor dataset is downloaded directly from the site (\url{http://db.csail.mit.edu/labdata/labdata.html}) and prepared (cleaned) by this function.
}
}
\author{
Dalibor Krleža
}
\references{
[1] Krleža D, Vrdoljak B, and Brčić M, Statistical hierarchical clustering algorithm for outlier detection in evolving data streams, \emph{Machine Learning}, Sep. 2020
}