SigmaIndex <- function(theta=3, neighborhood=9, balance=FALSE) {
  e1 <- new.env() # we need this to maintain the state of the stream generator
  e1$RObj <- new(SigmaIndex_R,theta,neighborhood,balance)
  
  si <- list(description = "Sigma-Index", env = e1)
  class(si) <- c("SigmaIndex")
  si
}

convertFromDSD <- function(x, total_elements=500, theta=3, neighborhood=9, balance=FALSE) {
  if(!"DSD_Gaussians" %in% class(x)) stop("Conversion from: [DSD_Gaussians] is currently supported")
  e1 <- new.env() # we need this to maintain the state of the stream generator
  e1$RObj <- new(SigmaIndex_R,theta,neighborhood,balance)
  if(is(x, "DSD_Gaussians")) {
    els <- as.integer(total_elements / x$k)
    ix <- 1
    for(i in 1:nrow(x$mu)) {
      e1$RObj$addPopulation(paste0(ix),unlist(x$mu[i,]),x$sigma[[i]],els)
      els <- els + 1
      ix <- ix + 1
    }
    if(x$o>0) {
      out_cov <- matrix(data=rep(0,ncol(x$mu)^2),nrow=ncol(x$mu),ncol=ncol(x$mu))
      diag(out_cov) <- x$outs_vv
      for(i in 1:nrow(x$outs)) {
        e1$RObj$addPopulation(paste0(ix),unname(unlist(x$outs[i,])),out_cov,1)
        ix <- ix + 1
      }
    }
  }
  
  si <- list(description = "Sigma-Index", env = e1)
  class(si) <- c("SigmaIndex")
  si
}

addPopulation <- function(x, id, mean, covariance, elements=1, ...) UseMethod("addPopulation");
addPopulation.SigmaIndex <- function(x, id, mean, covariance, elements=1, ...) {
  if(!is.character(id)) stop("Identifier must be a string")
  if(!is.matrix(covariance)) stop("Supplied covariance must be a matrix")
  if(!is.vector(mean)) stop("Supplied mean value must be a vector")
  if(length(mean)!=nrow(covariance) || length(mean)!=ncol(covariance))
    stop("Mean value and covariance dimensionalities do not match")
  x$env$RObj$addPopulation(id, mean, covariance, elements)
}

addDataPoints <- function(x, id, data_points, ...) UseMethod("addDataPoints");
addDataPoints.SigmaIndex <- function(x, id, data_points, ...) {
  if((is.matrix(data_points) || is.data.frame(data_points)) &&
     (is.list(id) || is.vector(id))) {
    if(is.list(id)) id <- unlist(id)
    for(i in 1:nrow(data_points)) {
      t1 <- data_points[i,,drop=T]
      if(is.list(t1)) t1 <- unlist(t1)
      x$env$RObj$addDataPoint(paste0(id[[i]]), t1, list())
    }
  } else if(is.list(data_points)) x$env$RObj$addDataPoint(paste0(id), unlist(data_points), list())
  else if(is.vector(data_points)) x$env$RObj$addDataPoint(paste0(id), data_points, list())
  else stop("Data points type missmatch")
}

addDataPointsInc <- function(x, data_points, query_results, ...) UseMethod("addDataPointsInc");
addDataPointsInc.SigmaIndex <- function(x, data_points, query_results, ...) {
  if(!is.list(query_results)) stop("Query results must be a list returned by the queryDataPoints method")
  if(!is.data.frame(data_points)) stop("Data points must be a data frame")
  if(length(query_results)!=nrow(data_points)) stop("Query results must have the same number of results as data points rows")
  for(i in 1:length(query_results)) {
    item <- query_results[[i]]
    if(!is.null(item$classified) && length(item$classified)>0) {
      id <- names(item$classified)[[1]]
      neighborhood <- list()
      if(!is.null(item$classified) && length(item$classified)>1)
        for(i in 2:length(item$classified))
          append(neighborhood, list(id=names(item$classified)[[i]],distance=item$classified[[i]]))
      for(id_neigh in names(item$neighborhood))
        append(neighborhood, list(id=id_neigh,distance=item$neighborhood[[id_neigh]]))
      t1 <- data_points[i,,drop=T]
      if(is.list(t1)) t1 <- unlist(t1)
      x$env$RObj$addDataPoint(id, t1, neighborhood)
    }
  }
}

queryDataPoints <- function(x, data_points, ...) UseMethod("queryDataPoints");
queryDataPoints.SigmaIndex <- function(x, data_points, ...) {
  res <- list()  
  if(is.matrix(data_points) || is.data.frame(data_points)) {
    for(i in 1:nrow(data_points)) {
      t1 <- data_points[i,,drop=T]
      if(is.list(t1)) t1 <- unlist(t1)
      res[[length(res)+1]] <- x$env$RObj$queryDataPoint(t1)
    }
  } else if(is.list(data_points)) res[[length(res)+1]] <- x$env$RObj$queryDataPoint(unlist(data_points))
  else if(is.vector(data_points)) res[[length(res)+1]] <- x$env$RObj$queryDataPoint(data_points)
  if(length(res)==0) warning("No classification returned from the sigma-index, means it is a new outlier")
  return(res)
}

getTotalPopulationNumber <- function(x, ...) UseMethod("getTotalPopulationNumber");
getTotalPopulationNumber.SigmaIndex <- function(x, ...) {
  x$env$RObj$getTotalPopulationNumber()
}

print <- function(x, ...) UseMethod("print");
print.SigmaIndex <- function(x, ...) {
  x$env$RObj$print()
}

getPopulations <- function(x, ...) UseMethod("getPopulations");
getPopulations.SigmaIndex <- function(x, ...) {
  x$env$RObj$getPopulations()
}

removePopulation <- function(x, id, ...) UseMethod("removePopulation");
removePopulation.SigmaIndex <- function(x, id, ...) {
  x$env$RObj$removePopulation(id)
}

resetStatistics <- function(x, ...) UseMethod("resetStatistics");
resetStatistics.SigmaIndex <- function(x, ...) {
  x$env$RObj$resetStatistics()
}

getHistogram <- function(x, ...) UseMethod("getHistogram");
getHistogram.SigmaIndex <- function(x, ...) {
  x$env$RObj$getHistogram()
}

getStatistics <- function(x, ...) UseMethod("getStatistics");
getStatistics.SigmaIndex <- function(x, ...) {
  x$env$RObj$getStatistics()
}