DSD_SHCGaussianGenerator <- function(clusters = 2, outliers = 0, movingClusters = 0, initMovingClusterSteps = 5,
                                     minDimensionValues = c(0,0), maxDimensionValues = c(100,100), initPopulation = 500,
                                     cacheAndRepeat = FALSE, minVariance = 0.5, maxVariance = 8, theta = 3.5,
                                     virtualVariance = 1.0, dimensions = 2, minCorrelation = 0.0, maxCorrelation = 0.2,
                                     generateVariableClusterPopulations = FALSE, checkOverlapping = TRUE,
                                     checkOverlapping_MD = FALSE, regenerateWithOverflow = TRUE) {
  if(!requireNamespace("Matrix",quietly=T)) stop("Obtain Matrix package to use this DSD")
  if(!requireNamespace("uuid",quietly=T)) stop("Obtain uuid package to use this DSD")
  if(!requireNamespace("MASS",quietly=T)) stop("Obtain MASS package to use this DSD")
  if(!requireNamespace("parallel",quietly=T)) stop("Obtain parallel package to use this DSD")
  if(!requireNamespace("rlist",quietly=T)) stop("Obtain rlist package to use this DSD")
  
  p1 <- initPopulation - outliers
  p2 <- trunc(p1 / (clusters + movingClusters * initMovingClusterSteps))
  env <- environment()
  generator <- SHC_MVRGenerator(dimension=2, minDimensionValues=minDimensionValues, maxDimensionValues=maxDimensionValues,
                                minVariance=minVariance, maxVariance=maxVariance, theta=theta, virtualVariance=virtualVariance,
                                minCorrelation=minCorrelation, maxCorrelation=maxCorrelation, clusters = clusters, 
                                clusterPopulation = p2, outliers = outliers, movingClusters = movingClusters, 
                                initMovingClusterSteps = initMovingClusterSteps, totalPopulation = initPopulation,
                                dimensionNamePrefix = "X", variableClusterPopulation = generateVariableClusterPopulations,
                                checkOverlapping = checkOverlapping, checkOverlapping_MD = checkOverlapping_MD)
  pointer <- 1
  rm(envir = env, list = c("p1", "p2", "clusters", "outliers", "movingClusters", "minDimensionValues",
                           "maxDimensionValues", "initPopulation", "minVariance", "maxVariance",
                           "theta", "virtualVariance", "dimensions", "minCorrelation", "maxCorrelation",
                           "initMovingClusterSteps", "checkOverlapping", "checkOverlapping_MD"))
  this <- list(
    description = "SHC multivariate synthetic stream",
    env = env
  )
  assign("this", this, env)
  class(this) <- c("DSD_SHCGaussianGenerator", "DSD_OutlierGenerator", "DSD_R", "DSD_data.frame", "DSD")
  return(this)
}

get_points.DSD_SHCGaussianGenerator <- function(x, n=1, outofpoints=c("stop", "warn", "ignore"), cluster = FALSE,
                                                class = FALSE, outlier = FALSE, ...) {
  if(class) cluster <- FALSE
  print(paste("*** Getting points from the stream (", n, ", class=", class, ", cluster=", cluster, ")", sep = ""))
  df <- data.frame()
  pointer <- x$env$pointer
  generator <- x$env$generator
  data <- generator$population
  points <- data[,!(colnames(data) %in% .reserved),drop=FALSE]
  czes <- generator$classes
  car <- x$env$cacheAndRepeat
  pnew <- pointer + (n - 1)
  if(!car && pnew>nrow(points)) {
    diff <- pnew - nrow(points) + 1
    generator$generate(n=diff,regenerate=x$env$regenerateWithOverflow)
    czes <- generator$classes
    data <- generator$population
    points <- data[,!(colnames(data) %in% .reserved),drop=FALSE]
  } else if(car && pnew > nrow(points))
    pnew <- nrow(points)
  clusters <- c()
  outliers <- c()
  cz <- c()
  df <- points[pointer:pnew,]
  for(i in pointer:pnew) {
    element <- data[i,]
    if(class) {
      cindex <- as.integer(which(czes==element$class))
      if(length(cindex) > 0)
        cz[[length(cz)+1]] <- cindex
      else
        cz[[length(cz)+1]] <- 0
    }
    if(cluster) clusters[[length(clusters)+1]] <- which(czes==element$class)
    if(outlier) {
      if(element$isOutlier) outliers[[length(outliers)+1]] <- TRUE
      else outliers[[length(outliers)+1]] <- FALSE
    }
  }
  if(class) df <- cbind(df, class=cz)
  if(!car) pointer <- pnew + 1
  else pointer <- 1
  x$env$pointer <- pointer
  if(cluster) attr(df, "cluster") <- unlist(clusters)
  if(outlier) attr(df, "outlier") <- unlist(outliers)
  return(df)
}

reset_stream.DSD_SHCGaussianGenerator <- function(dsd, pos=1) {
  dsd$env$pointer <- pos
}

plot.DSD_SHCGaussianGenerator <- function(x, n = 500, ...) {
  loadNamespace("ggplot2")
  d <- get_points(x, n, cluster = TRUE)
  assignment <- attr(d, "cluster")
  df <- cbind(d, c = assignment)
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x="X1", y="X2", color=c)) +
    ggplot2::geom_point(shape = 1, size = 2.5, show.legend = FALSE, stroke = 1.5) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                                         panel.grid.minor = ggplot2::element_blank(),
                                         axis.ticks = ggplot2::element_blank(), 
                                         axis.text = ggplot2::element_blank(), 
                                         axis.title = ggplot2::element_blank())
  return(plot)
}