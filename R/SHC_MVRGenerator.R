SHC_MVRGenerator <- setRefClass("SHC_MVRGenerator",
                                fields = list(
                                  clusters = "numeric",
                                  outliers = "numeric",
                                  movingClusters = "numeric",
                                  clusterPopulation = "numeric",
                                  movingClusterSteps = "numeric",
                                  totalPopulation = "numeric",
                                  dimension = "numeric",
                                  minDimensionValues = "numeric",
                                  maxDimensionValues = "numeric",
                                  checkOverlapping = "logical",
                                  checkOverlapping_MD = "logical",
                                  minVariance = "numeric",
                                  maxVariance = "numeric",
                                  minCorrelation = "numeric",
                                  maxCorrelation = "numeric",
                                  theta = "numeric",
                                  virtualVariance = "numeric",
                                  shuffle = "logical",
                                  population = "ANY",
                                  definitions = "ANY",
                                  classes = "ANY",
                                  oldOutliers = "ANY",
                                  dimNamePrefix = "character",
                                  variableClusterPopulation = "logical"
                                ),
                                methods = list(
                                  initialize = function(dimension = 2, minDimensionValues = c(0,0), maxDimensionValues = c(100,100),
                                                        checkOverlapping = TRUE, minVariance = 0.5, maxVariance = 16,
                                                        minCorrelation = 0.0, maxCorrelation = 0.2, theta = 3.5,
                                                        virtualVariance = 1.0, clusters = 1, movingClusters = 0, 
                                                        clusterPopulation = 100, outliers = 0, shuffle = TRUE, 
                                                        initMovingClusterSteps = 5, predefinedMovingClusterPaths = list(), 
                                                        totalPopulation = -1, dimensionNamePrefix = "", variableClusterPopulation = FALSE,
                                                        checkOverlapping_MD = FALSE) {
                                    clusters <<- clusters
                                    outliers <<- outliers
                                    movingClusters <<- movingClusters
                                    clusterPopulation <<- clusterPopulation
                                    movingClusterSteps <<- initMovingClusterSteps
                                    totalPopulation <<- totalPopulation
                                    dimension <<- dimension
                                    minDimensionValues <<- minDimensionValues
                                    maxDimensionValues <<- maxDimensionValues
                                    checkOverlapping <<- checkOverlapping
                                    checkOverlapping_MD <<- checkOverlapping_MD
                                    minVariance <<- minVariance
                                    maxVariance <<- maxVariance
                                    minCorrelation <<- minCorrelation
                                    maxCorrelation <<- maxCorrelation
                                    theta <<- theta
                                    virtualVariance <<- virtualVariance
                                    shuffle <<- shuffle
                                    population <<- c()
                                    definitions <<- c()
                                    classes <<- c()
                                    oldOutliers <<- c()
                                    dimNamePrefix <<- dimensionNamePrefix
                                    variableClusterPopulation <<- variableClusterPopulation
                                    
                                    mcv1 <- c()
                                    if(movingClusters > 0 && length(predefinedMovingClusterPaths) > 0 && length(predefinedMovingClusterPaths) != movingClusters)
                                      stop(paste("list of the predefined moving cluster paths (", length(predefinedMovingClusterPaths),") must match the number of moving clusters (", movingClusters, ")", sep = ""))
                                    if(movingClusters > 0) {
                                      print("*** Generating moving clusters")
                                      for(i in 1:movingClusters) {
                                        done <- FALSE
                                        iterCounter <- 0
                                        while(!done) {
                                          if(length(predefinedMovingClusterPaths) == 0)
                                            mcv2 <- generateMovingCluster(clusterPopulation, steps = movingClusterSteps)
                                          else {
                                            path <- predefinedMovingClusterPaths[[i]]
                                            mcv2 <- generateMovingCluster(clusterPopulation, steps = movingClusterSteps, predefinedPoints = path)
                                          }
                                          done <- TRUE
                                          if(checkOverlapping && length(mcv1) > 0)
                                            for(z in 1:length(mcv2)) {
                                              mcv2sp <- mcv2[[z]]
                                              for(i in 1:length(mcv1)) {
                                                mcvd <- mcv1[[i]]
                                                for(j in 1:length(mcvd)) {
                                                  mcvsp <- mcvd[[j]]
                                                  if(checkOverlapping_MD) {
                                                    if(checkOverlappingClusters_MD(mcv2sp$mean, mcv2sp$covariance, clusterPopulation,
                                                                                   mcvsp$mean, mcvsp$covariance, clusterPopulation)) {
                                                      done <- FALSE
                                                      break()
                                                    }
                                                  } else {
                                                    if(checkOverlappingClusters(mcv2sp$mean, mcv2sp$covariance, clusterPopulation,
                                                                                mcvsp$mean, mcvsp$covariance, clusterPopulation)) {
                                                      done <- FALSE
                                                      break()
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          iterCounter <- iterCounter + 1
                                          if(iterCounter > 200)
                                            stop("problem unsolvable, increase dimension limits, or decrease number of clusters and outliers")
                                        }
                                        mcv1[[length(mcv1)+1]] <- mcv2
                                        lmcv2 <- mcv2[[length(mcv2)]]
                                        classes <<- c(classes, lmcv2$id)
                                        definitions[[length(definitions)+1]] <<- lmcv2
                                        definitions[[length(definitions)]]$population <<- NULL
                                      }
                                    }
                                    cv1 <- c()
                                    if(clusters > 0) {
                                      print("*** Generating stationary clusters")
                                      for(i in 1:clusters) {
                                        done <- FALSE
                                        iterCounter <- 0
                                        while(!done) {
                                          if(variableClusterPopulation)
                                            cv2 <- generateCluster(as.integer(runif(1,clusterPopulation,3*clusterPopulation)))
                                          else
                                            cv2 <- generateCluster(clusterPopulation)
                                          done <- TRUE
                                          if(checkOverlapping && length(cv1) > 0)
                                            for(j in 1:length(cv1)) {
                                              cv3 <- cv1[[j]]
                                              if(checkOverlapping_MD) {
                                                if(checkOverlappingClusters_MD(cv2$mean, cv2$covariance, cv2$n, 
                                                                               cv3$mean, cv3$covariance, cv2$n)) {
                                                  done <- FALSE
                                                  break()
                                                }
                                              } else {
                                                if(checkOverlappingClusters(cv2$mean, cv2$covariance, cv2$n, 
                                                                            cv3$mean, cv3$covariance, cv2$n)) {
                                                  done <- FALSE
                                                  break()
                                                }
                                              }
                                            }
                                          if(checkOverlapping && length(mcv1) > 0)
                                            for(i in 1:length(mcv1)) {
                                              mcvd <- mcv1[[i]]
                                              for(j in 1:length(mcvd)) {
                                                mcvsp <- mcvd[[j]]
                                                if(checkOverlapping_MD) {
                                                  if(checkOverlappingClusters_MD(cv2$mean, cv2$covariance, cv2$n, 
                                                                                 mcvsp$mean, mcvsp$covariance, clusterPopulation)) {
                                                    done <- FALSE
                                                    break()
                                                  }
                                                } else {
                                                  if(checkOverlappingClusters(cv2$mean, cv2$covariance, cv2$n, 
                                                                              mcvsp$mean, mcvsp$covariance, clusterPopulation)) {
                                                    done <- FALSE
                                                    break()
                                                  }
                                                }
                                              }
                                            }
                                          iterCounter <- iterCounter + 1
                                          if(iterCounter > 200)
                                            stop("problem unsolvable, increase dimension limits, or decrease number of clusters and outliers")
                                        }
                                        classes <<- c(classes, cv2$id)
                                        cv1[[length(cv1)+1]] <- cv2
                                        definitions[[length(definitions)+1]] <<- cv2
                                        definitions[[length(definitions)]]$population <<- NULL
                                      }
                                    }
                                    diff <- totalPopulation - ((movingClusters * movingClusterSteps + clusters) * clusterPopulation + outliers)
                                    if(diff > 0) {
                                      if(movingClusters > 0) {
                                        mcv2 <- mcv1[[length(mcv1)]]
                                        lmcv2 <- mcv2[[length(mcv2)]]
                                        pop <- MASS::mvrnorm(diff, lmcv2$mean, lmcv2$covariance)
                                        if(diff==1) pop <- matrix(pop, ncol=dimension)
                                        pop <- data.frame(pop)
                                        colnames(pop) <- paste0(dimNamePrefix,as.character(seq(1,ncol(pop))))
                                        lmcv2$population <- rbind(lmcv2$population, pop)
                                        mcv2[[length(mcv2)]] <- lmcv2
                                        mcv1[[length(mcv1)]] <- mcv2
                                      } else if(clusters > 0) {
                                        cv2 <- cv1[[length(cv1)]]
                                        pop <- MASS::mvrnorm(diff, cv2$mean, cv2$covariance)
                                        if(diff==1) pop <- matrix(pop, ncol=dimension)
                                        pop <- data.frame(pop)
                                        colnames(pop) <- paste0(dimNamePrefix,as.character(seq(1,ncol(pop))))
                                        cv2$population <- rbind(cv2$population, pop)
                                        cv1[[length(cv1)]] <- cv2
                                      } else
                                        stop(paste0("total population (", totalPopulation, ") is bigger than data generator can generate"))
                                    }
                                    ov1 <- c()
                                    if(outliers > 0)
                                      for(i in 1:outliers) {
                                        print(paste("*** Generating outlier ", i, sep = ""))
                                        done <- FALSE
                                        iterCounter <- 0
                                        while(!done) {
                                          ov2 <- generateOutlier()
                                          done <- TRUE
                                          if(checkOverlapping && length(cv1) > 0)
                                            for(j in 1:length(cv1)) {
                                              cv3 <- cv1[[j]]
                                              if(checkOverlappingClusterToPoint(cv3$mean, cv3$covariance, cv3$n, ov2$outlier)) {
                                                done <- FALSE
                                                break()
                                              }
                                            }
                                          if(done && checkOverlapping && length(mcv1) > 0)
                                            for(i in 1:length(mcv1)) {
                                              mcvd <- mcv1[[i]]
                                              for(j in 1:length(mcvd)) {
                                                mcvsp <- mcvd[[j]]
                                                if(checkOverlappingClusterToPoint(mcvsp$mean, mcvsp$covariance, clusterPopulation, ov2$outlier)) {
                                                  done <- FALSE
                                                  break()
                                                }
                                              }
                                            }
                                          if(done && checkOverlapping && length(ov1) > 0)
                                            for(j in 1:length(ov1)) {
                                              ov3 <- ov1[[j]]
                                              if(checkOverlappingPointToPoint(ov3$outlier, ov2$outlier)) {
                                                done <- FALSE
                                                break()
                                              }
                                            }
                                          iterCounter <- iterCounter + 1
                                          if(iterCounter > 200)
                                            stop("problem unsolvable, increase dimension limits, or decrease number of clusters and outliers")
                                        }
                                        classes <<- c(classes, ov2$id)
                                        ov1[[length(ov1)+1]] <- ov2
                                        definitions[[length(definitions)+1]] <<- ov2
                                      }
                                    pop <- data.frame()
                                    if(clusters > 0)
                                      for(x in 1:clusters) {
                                        clus_def <- cv1[[x]]
                                        clus_def_n <- nrow(clus_def$population)
                                        temp_pop <- cbind(clus_def$population, class=rep(clus_def$id, clus_def_n), stopHere=rep(FALSE, clus_def_n), isOutlier=rep(FALSE, clus_def_n))
                                        pop <- rbind(pop, temp_pop)
                                      }
                                    if(outliers > 0)
                                      for(x in 1:outliers) {
                                        out_df <- data.frame(t(ov1[[x]]$outlier))
                                        colnames(out_df) <- paste0(dimNamePrefix,as.character(seq(1, ncol(out_df))))
                                        out_df <- cbind(out_df, class=ov1[[x]]$id, stopHere=FALSE, isOutlier=TRUE)
                                        pop <- rbind(pop, out_df)
                                      }
                                    if(shuffle)
                                      pop <- pop[sample(1:nrow(pop)),]
                                    if(movingClusters > 0 && length(mcv1) > 0) {
                                      pop2 <- data.frame()
                                      popslice <- trunc(nrow(pop) / movingClusterSteps)
                                      for(j in 1:movingClusterSteps) {
                                        pop1 <- data.frame()
                                        for(i in 1:movingClusters) {
                                          mcvsp <- mcv1[[i]][[j]]
                                          mcvsp_n <- nrow(mcvsp$population)
                                          pop1 <- rbind(pop1, cbind(mcvsp$population, class=rep(mcvsp$id, mcvsp_n), stopHere=rep(FALSE, mcvsp_n), isOutlier=rep(FALSE, mcvsp_n)))
                                        }
                                        li <- (j - 1) * popslice + 1
                                        if(j < movingClusterSteps) ui <- j * popslice 
                                        else ui <- nrow(pop)
                                        pop1 <- rbind(pop1, pop[li:ui,])
                                        if(shuffle)
                                          pop1 <- pop1[sample(1:nrow(pop1)),]
                                        pop1[nrow(pop1),]$stopHere <- TRUE
                                        pop2 <- rbind(pop2, pop1)
                                      }
                                      population <<- pop2
                                    } else population <<- pop
                                  }
                                ))

SHC_MVRGenerator$methods(list(
  generateCluster = function(n, mu = NA, covariance = NA) {
    cdef <- generateClusterDefinition(mu = mu, covariance = covariance)
    pop <- MASS::mvrnorm(n, cdef$mean, cdef$covariance)
    if(n==1) pop <- matrix(pop, ncol=dimension)
    pop <- data.frame(pop)
    colnames(pop) <- paste0(dimNamePrefix,as.character(seq(1,ncol(pop))))
    cdef$population <- pop
    cdef$n <- n
    return(cdef)
  },
  
  generateClusterDefinition = function(mu = NA, covariance = NA, id = NA) {
    if(all(is.na(mu)))
      mu <- runif(dimension, min = minDimensionValues, max = maxDimensionValues)
    if(all(is.na(covariance)))
      covariance <- generateRandomCovariance()
    if(all(is.na(id)))
      id <- paste("CLUS(", uuid::UUIDgenerate(), ")", sep = "")
    list(id = id, mean = mu, covariance = covariance, type = "stationary")
  },
  
  generateRandomCovariance = function() {
    covariance <- matrix(data = c(0, length = (dimension * dimension)), ncol = dimension, nrow = dimension)
    for(i in 1:dimension) {
      for(j in i:dimension) {
        if(i==j) covariance[i, j] <- runif(1, min = minVariance, max = maxVariance)
        else {
          covariance[i, j] <- runif(1, min = minCorrelation, max = maxCorrelation) * sqrt(covariance[i, i]) * sqrt(covariance[j, j])
          covariance[j, i] <- covariance[i, j]
        }
      }
    }
    covDPO <- Matrix::nearPD(covariance)$mat
    covariance <- matrix(data = covDPO@x, ncol = covDPO@Dim[[1]], nrow = covDPO@Dim[[2]])
    covariance
  },
  
  generateMovingCluster = function(n, steps = 5, predefinedPoints = NA) {
    mcdefs <- generateMovingClusterDefinitions(steps = steps, predefinedPoints = predefinedPoints)
    pops <- c()
    for(i in 1:length(mcdefs)) {
      mcdef <- mcdefs[[i]]
      pop <- MASS::mvrnorm(n, mcdef$mean, mcdef$covariance)
      if(n==1) pop <- matrix(pop, ncol=dimension)
      pop <- data.frame(pop)
      colnames(pop) <- paste0(dimNamePrefix,as.character(seq(1,ncol(pop))))
      mcdef$population <- pop
      mcdef$n <- n
      pops[[i]] <- mcdef
    }
    pops
  },
  
  generateMovingClusterDefinitions = function(steps = 5, predefinedPoints = NA) {
    if(!all(is.na(predefinedPoints))) {
      startPoint <- predefinedPoints$startPoint
      endPoint <- predefinedPoints$endPoint
    } else {
      startPoint <- runif(dimension, min = minDimensionValues, max = maxDimensionValues)
      allsteps <- rep(as.double(steps)*as.double(maxVariance),length(startPoint))
      tminDim <- startPoint - allsteps
      for(i in 1:dimension)
        if(tminDim[[i]] < minDimensionValues[[i]]) tminDim[[i]] <- minDimensionValues[[i]]
      tmaxDim <- startPoint + allsteps
      for(i in 1:dimension)
        if(tmaxDim[[i]] > maxDimensionValues[[i]]) tmaxDim[[i]] <- maxDimensionValues[[i]]
      endPoint <- runif(dimension, min = tminDim, max = tmaxDim)
    }
    deltaPoint <- (endPoint - startPoint) / as.double(steps)
    pops <- list()
    covariance <- generateRandomCovariance()
    mclusid <- paste0("MCLUS(", uuid::UUIDgenerate(), ")")
    if(all(is.na(predefinedPoints))) {
      print(paste("*** Adjusting moving cluster ", mclusid, sep = ""))
      closeEnough <- FALSE
      d <- 1.0
      doExpand <- FALSE
      while(!closeEnough) {
        p2 <- startPoint + d * deltaPoint
        if(!checkOverlappingClusters(startPoint, covariance, 1000, p2, covariance, 1000)) {
          d <- d - 0.1
          if(doExpand) 
            closeEnough <- TRUE
        } else {
          if(d >= 1.0) {
            doExpand <- TRUE
            d <- d + 0.1
          } else closeEnough <- TRUE
        }
      }
      d <- 0.45 * d
      deltaPoint <- d * deltaPoint
    }
    for(i in 0:(steps-1)) {
      point <- (startPoint + i * deltaPoint)
      pops[[i+1]] <- list(id = mclusid, mean = point, covariance = covariance, type = "moving", movement = deltaPoint)
    }
    pops
  },
  
  generateOutlier = function() {
    return(list(id = paste0("OUTL(", uuid::UUIDgenerate(), ")"),
                outlier = runif(dimension, min = minDimensionValues, max = maxDimensionValues),
                type = "outlier"))
  },
  
  generateClusterPopulation = function(n, def) {
    if(def$type == "moving")
      def$mean = def$mean + def$movement
    pop <- MASS::mvrnorm(n, def$mean, def$covariance)
    if(n==1) pop <- matrix(pop, ncol=dimension)
    pop <- data.frame(pop)
    colnames(pop) <- paste0(dimNamePrefix,as.character(seq(1,ncol(pop))))
    outpop <- cbind(pop, class=rep(def$id, n), stopHere=rep(FALSE, n), isOutlier=rep(FALSE, n))
    #outpop <- c()
    #for(i in 1:nrow(pop)) {
    #  item <- list(value = pop[i,], original_clazz = def$id, stopHere = FALSE, isOutlier = FALSE)
    #  outpop[[length(outpop)+1]] <- item
    #}
    if(shuffle)
      outpop <- outpop[sample(1:nrow(outpop)),]
    list(population = outpop, definition = def)
  },
  
  generate = function(n = NA, regenerate = TRUE) {
    if(all(is.na(n)))
      if(totalPopulation>0)
        n <- totalPopulation
      else 
        n <- (movingClusters * movingClusterSteps + clusters) * clusterPopulation + outliers
    clusts <- length(definitions[which(sapply(definitions, function(x) x$type) != "outlier")])
    outs <- length(definitions) - clusts
    if(clusts > 0) {
      v1 <- trunc((n - outs) / clusts)
      remainder <- (n - outs) - (v1 * clusts)
    } else {
      v1 <- 0
      remainder <- 0
    }
    if(!regenerate) {
      sclusts <- length(definitions[which(sapply(definitions, function(x) x$type) == "stationary")])
      mclustDefs <- definitions[which(sapply(definitions, function(x) x$type) == "moving")]
      mclusts <- length(mclustDefs)
      clusters <<- sclusts
      movingClusters <<- mclusts
      outliers <<- outs
      clusterPopulation <<- v1+remainder
    }
    outpop <- data.frame()
    newdefs <- c()
    if(regenerate)
      for(def in definitions) 
        if(def$type == "outlier") 
          oldOutliers[[length(oldOutliers)+1]] <<- def
    for(def in definitions) {
      if(regenerate && def$type == "outlier") {
        done <- FALSE
        iterCounter <- 0
        while(!done) {
          out <- generateOutlier()
          done <- TRUE
          if(checkOverlapping && clusts>0)
            for(def2 in definitions) {
              if(def2$type %in% c("stationary", "moving") &&
                   checkOverlappingClusterToPoint(def2$mean, def2$covariance, clusterPopulation, out$outlier)) 
                  done <- FALSE
            }
          if(done && checkOverlapping && length(oldOutliers)>0)
            for(oo in oldOutliers) {
              if(checkOverlappingPointToPoint(oo$outlier, out$outlier)) 
                done <- FALSE
            }
          if(done && checkOverlapping && length(newdefs)>0)
            for(def2 in newdefs) {
              if(def2$type %in% c("outlier") &&
                 checkOverlappingPointToPoint(def2$outlier, out$outlier)) 
                done <- FALSE
            }
          iterCounter <- iterCounter + 1
          if(iterCounter > 200)
            stop("problem unsolvable, increase dimension limits, or decrease number of clusters and outliers")
        }
        out_df <- data.frame(t(out$outlier))
        colnames(out_df) <- paste0(dimNamePrefix,as.character(seq(1,ncol(out_df))))
        outpop <- rbind(outpop, cbind(out_df, class=out$id, stopHere=FALSE, isOutlier=TRUE))
        #outpop[[length(outpop)+1]] <- list(value = out$outlier, original_clazz = out$id, stopHere = FALSE, isOutlier = TRUE)
        if(!out$id %in% classes)
          classes <<- c(classes, out$id)
        newdefs[[length(newdefs)+1]] <- out
      } else if(!regenerate && def$type == "outlier") {
        out_df <- data.frame(t(def$outlier))
        colnames(out_df) <- paste0(dimNamePrefix,as.character(seq(1,ncol(out_df))))
        outpop <- rbind(outpop, cbind(out_df, class=def$id, stopHere=FALSE, isOutlier=TRUE))
        #outpop[[length(outpop)+1]] <- list(value = def$outlier, original_clazz = def$id, stopHere = FALSE, isOutlier = TRUE)
        if(!def$id %in% classes)
          classes <<- c(classes, def$id)
        newdefs[[length(newdefs)+1]] <- def
      } else {
        if(variableClusterPopulation)
          popres <- generateClusterPopulation(as.integer(runif(1,v1 + remainder,3*(v1 + remainder))), def)
        else
          popres <- generateClusterPopulation(v1 + remainder, def)
        outpop <- rbind(outpop, popres$population)
        remainder <- 0
        #outpop <- append(outpop, popres$population)
        if(!def$id %in% classes)
          classes <<- c(classes, def$id)
        newdefs[[length(newdefs)+1]] <- popres$definition
      }
    }
    if(shuffle)
      outpop <- outpop[sample(1:nrow(outpop)),]
    
    population <<- rbind(population, outpop)
    definitions <<- newdefs
  },
  
  plot = function(fromPos = -1, toPos = -1) {
    loadNamespace("ggplot2")
    if(fromPos>0 && toPos>0 && toPos<fromPos) stop("must be: fromPos<=toPos")
    if(nrow(population)<1) stop("No elements in data list")
    if(fromPos>0)
      if(toPos>0)
        p <- population[fromPos:toPos,]
      else
        p <- population[fromPos:nrow(population),]
    else
      if(toPos>0)
        p <- population[1:toPos,]
      else
        p <- population
    tmp_data <- p[,!(colnames(population) %in% .reserved),drop=FALSE]
    if(ncol(tmp_data)!=2) stop("Data elements are not 2-dimensional")
    x <- c()
    y <- c()
    clazz <- c()
    for(i in 1:nrow(tmp_data)) {
      x <- c(x, tmp_data[i,1])
      y <- c(y, tmp_data[i,2])
      clazz <- c(clazz, p[i,]$class)
    }
    df <- data.frame(x, y, clazz)
    plot <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = clazz)) +
      ggplot2::geom_point(shape = 1, size = 2, show.legend = FALSE) +
      ggplot2::theme_minimal()
    print(plot)
  },
  
  calculateClusterFromPoint = function(newElement) {
    K <- length(newElement)
    N <- 0
    mean0 <- rep(0, K)
    mean1 <- shc_CalculateNewMean(mean0, newElement, N)
    covar0 <- matrix(0.0, nrow = K, ncol = K)
    covar1 <- shc_CalculateCovariance(mean0, mean1, covar0, N, newElement, FALSE)
    list(mean = mean1, covariance = covar1)
  },
  
  generateVirtualVariance = function(K) {
    vvar0 <- array(0.0, dim = K)
    for(i in 1:K)
      vvar0[i] <- virtualVariance
    return(vvar0)
  },
  
  checkOverlappingClusters = function(mu1, cov1, N1, mu2, cov2, N2) {
    item <- list(key = 0, 
                 c1Component = list(
                   mean = mu1, covariance = cov1, virtualVariance = generateVirtualVariance(length(mu1)), 
                   N = N1, isInversion = FALSE
                 ),
                 c2Component = list(
                   mean = mu2, covariance = cov2, virtualVariance = generateVirtualVariance(length(mu2)), 
                   N = N2, isInversion = FALSE
                 ),
                 c1th = theta,
                 c2th = theta)
    iter_mdi <- c()
    iter_mdi[[1]] <- item
    eval_1 <- parallel::mclapply(iter_mdi, calcCM)
    eval_2 <- rlist::list.filter(eval_1, measure <= 1.0)
    if(length(eval_2) > 0) return(TRUE)
    else return(FALSE)
  },
  
  checkOverlappingClusters_MD = function(mu1, cov1, N1, mu2, cov2, N2) {
    item <- list(key = 0, 
                 c1Component = list(
                   mean = mu1, covariance = cov1, virtualVariance = generateVirtualVariance(length(mu1)), 
                   N = N1, isInversion = FALSE
                 ),
                 c2Component = list(
                   mean = mu2, covariance = cov2, virtualVariance = generateVirtualVariance(length(mu2)), 
                   N = N2, isInversion = FALSE
                 ))
    iter_mdi <- c()
    iter_mdi[[1]] <- item
    eval_1 <- parallel::mclapply(iter_mdi, calcMMM)
    eval_2 <- rlist::list.filter(eval_1, measure <= theta)
    if(length(eval_2) > 0) return(TRUE)
    else return(FALSE)
  },
  
  checkOverlappingClusterToPoint = function(mu1, cov1, N1, point) {
    res1 <- calculateClusterFromPoint(point)
    item <- list(key = 0, 
                 c1Component = list(
                   mean = mu1, covariance = cov1, virtualVariance = generateVirtualVariance(length(mu1)), 
                   N = N1, isInversion = FALSE
                 ),
                 c2Component = list(
                   mean = res1$mean, covariance = res1$covariance, virtualVariance = generateVirtualVariance(length(point)), 
                   N = 1, isInversion = FALSE
                 ),
                 c1th = theta,
                 c2th = theta)
    iter_mdi <- c()
    iter_mdi[[1]] <- item
    eval_1 <- parallel::mclapply(iter_mdi, calcCM)
    eval_2 <- rlist::list.filter(eval_1, measure <= 1.0)
    if(length(eval_2) > 0) return(TRUE)
    else return(FALSE)
  },
  
  checkOverlappingPointToPoint = function(point1, point2) {
    res1 <- calculateClusterFromPoint(point1)
    res2 <- calculateClusterFromPoint(point2)
    item <- list(key = 0, 
                 c1Component = list(
                   mean = res1$mean, covariance = res1$covariance, virtualVariance = generateVirtualVariance(length(point1)), 
                   N = 1, isInversion = FALSE
                 ),
                 c2Component = list(
                   mean = res2$mean, covariance = res2$covariance, virtualVariance = generateVirtualVariance(length(point2)), 
                   N = 1, isInversion = FALSE
                 ),
                 c1th = theta,
                 c2th = theta)
    iter_mdi <- c()
    iter_mdi[[1]] <- item
    eval_1 <- parallel::mclapply(iter_mdi, calcCM)
    eval_2 <- rlist::list.filter(eval_1, measure <= 1.0)
    if(length(eval_2) > 0) return(TRUE)
    else return(FALSE)
  }
))

SHC_Predefined_MVRGenerator <- setRefClass("SHC_Predefined_MVRGenerator",
                                           contains = "SHC_MVRGenerator",
                                           methods = list(
                                             initialize = function(definitions, totalPopulation, dimension = 2, 
                                                                   minDimensionValues = c(0,0), maxDimensionValues = c(100,100),
                                                                   checkOverlapping = TRUE, minVariance = 0.5, maxVariance = 16,
                                                                   minCorrelation = 0.0, maxCorrelation = 0.2, theta = 3.5,
                                                                   virtualVariance = 1.0, shuffle = TRUE, dimensionNamePrefix = "") {
                                               totalPopulation <<- totalPopulation
                                               dimension <<- dimension
                                               minDimensionValues <<- minDimensionValues
                                               maxDimensionValues <<- maxDimensionValues
                                               checkOverlapping <<- checkOverlapping
                                               minVariance <<- minVariance
                                               maxVariance <<- maxVariance
                                               minCorrelation <<- minCorrelation
                                               maxCorrelation <<- maxCorrelation
                                               theta <<- theta
                                               virtualVariance <<- virtualVariance
                                               shuffle <<- shuffle
                                               population <<- c()
                                               updateDefinitions(definitions)
                                               classes <<- c()
                                               oldOutliers <<- c()
                                               dimNamePrefix <<- dimensionNamePrefix
                                               
                                               generate()
                                             }
                                           ))

SHC_Predefined_MVRGenerator$methods(list(
  generate = function(n = NA, regenerate = FALSE) {
    callSuper(n, regenerate)
  },
  
  updateDefinitions = function(definitions) {
    definitions <<- definitions
    checkCovariances()
  },
  
  updateDefinition = function(i, attr, value) {
    definitions[[i]][[attr]] <<- value
    checkCovariances()
  },
  
  addDefinition = function(def) {
    definitions[[length(definitions)+1]] <<- def
    checkCovariances()
  },
  
  checkCovariances = function() {
    for(i in 1:length(definitions)) {
      def <- definitions[[i]]
      if(def$type %in% c("moving", "stationary") && all(is.na(def$covariance))) {
        cov <- generateRandomCovariance()
        definitions[[i]]$covariance <<- cov
      }
    }
  }
))