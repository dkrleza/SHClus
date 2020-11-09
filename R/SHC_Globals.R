AgglomerationType <- list(NormalAgglomeration=0,AggresiveAgglomeration=1,RelaxedAgglomeration=2)
DriftType <- list(NormalDrift=0,FastDrift=1,SlowDrift=2,NoDrift=3,UltraFastDrift=4)
.reserved <- c("class","assigned","clazz","component","cids","isOutlier","stopHere","clazzAsOutlier")

calcCM <- function(item) {
  mean1 <- item$c1Component$mean
  covar1 <- item$c1Component$covariance
  vvar1 <- item$c1Component$virtualVariance
  N1 <- item$c1Component$N
  isInv1 <- item$c1Component$isInversion
  mean2 <- item$c2Component$mean
  covar2 <- item$c2Component$covariance
  vvar2 <- item$c2Component$virtualVariance
  N2 <- item$c2Component$N
  isInv2 <- item$c2Component$isInversion
  mdsep <- shc_MDSeparation(mean1, covar1, vvar1, N1, isInv1, item$c1th,
                            mean2, covar2, vvar2, N2, isInv2, item$c2th)
  list(key=item$key, measure=mdsep)
}

calcMMM <- function(item) {
  mean1 <- item$c1Component$mean
  covar1 <- item$c1Component$covariance
  vvar1 <- item$c1Component$virtualVariance
  N1 <- item$c1Component$N
  isInv1 <- item$c1Component$isInversion
  mean2 <- item$c2Component$mean
  covar2 <- item$c2Component$covariance
  vvar2 <- item$c2Component$virtualVariance
  N2 <- item$c2Component$N
  isInv2 <- item$c2Component$isInversion
  md_min <- shc_MutualMinMahalanobis(mean1, covar1, vvar1, N1, isInv1,
                                     mean2, covar2, vvar2, N2, isInv2)
  list(key=item$key, measure=md_min)
}

normalize.data.frame <- function(x, take_cols, class_col=NULL, outlier_col=NULL, multiplier=1) {
  spec_cols <- c()
  if(!is.null(class_col)) {
    if(is.character(class_col)) spec_cols <- c(spec_cols,which(colnames(x)==class_col))
    else spec_cols <- c(spec_cols,class_col)
  }
  if(!is.null(outlier_col)) {
    if(is.character(outlier_col)) spec_cols <- c(spec_cols,which(colnames(x)==outlier_col))
    else spec_cols <- c(spec_cols,outlier_col)
  }
  take_cols_i <- c()
  for(v1 in take_cols) {
    if(is.character(v1)) take_cols_i <- c(take_cols_i,which(colnames(x)==v1))
    else take_cols_i <- c(take_cols_i,v1)
  }
  
  tdf <- x[,setdiff(take_cols_i,spec_cols)]
  maxdiff <- 0
  for(col in colnames(tdf)) {
    cdiff <- max(tdf[,col])-min(tdf[,col])
    if(cdiff>maxdiff) maxdiff <- cdiff
  }
  for(col in colnames(tdf)) {
    m_low <- min(tdf[,col])
    m_upper <- max(tdf[,col])
    normalized <- (tdf[,col]-m_low)/maxdiff
    tdf[col] <- normalized * multiplier
    #tdf[,col] <- 
  }
  tdf <- cbind(tdf,x[,spec_cols])
  tdf
}