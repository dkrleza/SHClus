SHC_TestCase_ClustersAndOutliers_SigmaIndex <- function(n=100000) {
  set.seed(1000)
#  ds <- DSD_SHCGaussianGenerator(clusters = 50, outliers = 50, initPopulation = n, 
#                                 cacheAndRepeat = TRUE, maxDimensionValues = c(150,150),
#                                 generateVariableClusterPopulations = F,regenerateWithOverflow = FALSE)
#
  ds <- stream::DSD_Gaussians(k=50,o=50,separation_type="Mahalanobis",lim=c(0,180),
                      outlier_horizon=30000,variance_lim=8)

  print("#SHC-sequential")
  c1 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = FALSE)
  setPseudoOfflineCounter(c1,500)
  res1 <- evaluate_with_callbacks(c1, ds, n=n, measure = c("cRand","queryTime","updateTime",
                                                           "processTime","nodeCount",
                                                           "computationCostReduction"), 
                                  type="macro", callbacks=list(shc=SHCEvalCallback()), 
                                  single_pass_update=T, use_outliers=T)
  
  reset_stream(ds)
  
  print("#SHC-index1")
  c2 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 2)
  setPseudoOfflineCounter(c2,500)
  res2 <- evaluate_with_callbacks(c2, ds, n=n, measure = c("cRand","queryTime","updateTime",
                                                           "processTime","nodeCount",
                                                           "computationCostReduction"), 
                                  type="macro", callbacks=list(shc=SHCEvalCallback()), 
                                  single_pass_update=T, use_outliers=T)
  hist2 <- getHistogram(c2)
  pdf("./inst/shc_idx1_histogram.pdf",  width=3.2, height=4)
  par(mar=c(5, 4, 1, 1) + 0.1)
  barplot(unname(hist2[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)", 
          space=0, beside=T)
  axis(side=1,at=c(0,20,40,60,80,100))
  dev.off()
  
  reset_stream(ds)
  
  print("#SHC-index2")
  c3 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 3)
  setPseudoOfflineCounter(c3,500)
  res3 <- evaluate_with_callbacks(c3, ds, n=n, measure = c("cRand","queryTime","updateTime",
                                                           "processTime","nodeCount",
                                                           "computationCostReduction"), 
                                  type="macro", callbacks=list(shc=SHCEvalCallback()), 
                                  single_pass_update=T, use_outliers=T)
  hist3 <- getHistogram(c3)
  pdf("./inst/shc_idx2_histogram.pdf",  width=3.2, height=4)
  par(mar=c(5, 4, 1, 1) + 0.1)
  barplot(unname(hist3[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)", 
          space=0, beside=T)
  axis(side=1,at=c(0,20,40,60,80,100))
  dev.off()
  
  reset_stream(ds)
  
  print("#SHC-index3")
  c4 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 4)
  setPseudoOfflineCounter(c4,500)
  res4 <- evaluate_with_callbacks(c4, ds, n=n, measure = c("cRand","queryTime","updateTime",
                                                           "processTime","nodeCount",
                                                           "computationCostReduction"), 
                                  type="macro", callbacks=list(shc=SHCEvalCallback()), 
                                  single_pass_update=T, use_outliers=T)
  hist4 <- getHistogram(c4)
  pdf("./inst/shc_idx3_histogram.pdf",  width=3.2, height=4)
  par(mar=c(5, 4, 1, 1) + 0.1)
  barplot(unname(hist4[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)", 
          space=0, beside=T)
  axis(side=1,at=c(0,20,40,60,80,100))
  dev.off()
  
  reset_stream(ds)
  
  print("#SHC-balanced")
  c5 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 3, balancedSigmaIndex = T)
  setPseudoOfflineCounter(c5,500)
  res5 <- evaluate_with_callbacks(c5, ds, n=n, measure = c("cRand","queryTime","updateTime",
                                                           "processTime","nodeCount",
                                                           "computationCostReduction"), 
                                  type="macro", callbacks=list(shc=SHCEvalCallback()), 
                                  single_pass_update=T, use_outliers=T)
  hist5 <- getHistogram(c5)
  pdf("./inst/shc_idx_bal_histogram.pdf",  width=3.2, height=4)
  par(mar=c(5, 4, 1, 1) + 0.1)
  barplot(unname(hist5[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)", 
          space=0, beside=T)
  axis(side=1,at=c(0,20,40,60,80,100))
  dev.off()
  
  pdf("./inst/st_clusout_crand.pdf",  width=3, height=4)
  par(mar=c(6.8, 4, 1, 1) + 0.1)
  df1 <- data.frame('SHC sequential'=res1$cRand,'SHC index 1'=res2$cRand,
                    'SHC index 2'=res3$cRand,'SHC index 3'=res4$cRand,
                    'SHC balanced'=res5$cRand,check.names=F)
  mns <- colMeans(df1)
  df1 <- df1[,order(mns,decreasing=T)]
  boxplot(df1, las=2, ylab="Corrected Rand", ylim=c(0,1))
  dev.off()
  
  pdf("./inst/st_clusout_querytimes.pdf",  width=3, height=4)
  par(mar=c(6, 4, 1, 1) + 0.1)
  df1 <- data.frame('SHC sequential'=res1$queryTime,'SHC index 1'=res2$queryTime,
                    'SHC index 2'=res3$queryTime,'SHC index 3'=res4$queryTime,
                    'SHC balanced'=res5$queryTime,check.names=F)
  mns <- colMeans(df1)
  df1 <- df1[,order(mns)]
  barplot(colMeans(df1), las=2, cex.names=.85, ylab="Query time (ms)", ylim=c(0,(max(mns)*1.1)))
  dev.off()
  
  pdf("./inst/st_clusout_nc.pdf",  width=3, height=4)
  par(mar=c(6, 4, 1, 1) + 0.1)
  df1 <- data.frame('SHC sequential'=res1$nodeCount,'SHC index 1'=res2$nodeCount,
                    'SHC index 2'=res3$nodeCount,'SHC index 3'=res4$nodeCount,
                    'SHC balanced'=res5$nodeCount,check.names=F)
  mns <- colMeans(df1)
  df1 <- df1[,order(mns)]
  barplot(colMeans(df1), las=2, cex.names=.85, cex.axis=0.6, ylab="Nodes visited/calculated", ylim=c(0,(max(mns)*1.1)))
  dev.off()
  
  pdf("./inst/st_clusout_ccr.pdf",  width=3, height=4)
  par(mar=c(6, 4, 1, 1) + 0.1)
  df1 <- data.frame('SHC sequential'=res1$computationCostReduction,'SHC index 1'=res2$computationCostReduction,
                    'SHC index 2'=res3$computationCostReduction,'SHC index 3'=res4$computationCostReduction,
                    'SHC balanced'=res5$computationCostReduction,check.names=F)
  mns <- colMeans(df1)
  df1 <- df1[,order(mns)]
  barplot(colMeans(df1), las=2, cex.names=.85, ylab="Comp. cost reduction (%)", ylim=c(0,(max(mns)*1.1)))
  dev.off()
  
  pdf("./inst/st_clusout_ut.pdf",  width=3, height=4)
  par(mar=c(6, 4, 1, 1) + 0.1)
  df1 <- data.frame('SHC sequential'=res1$updateTime,'SHC index 1'=res2$updateTime,
                    'SHC index 2'=res3$updateTime,'SHC index 3'=res4$updateTime,
                    'SHC balanced'=res5$updateTime,check.names=F)
  mns <- colMeans(df1)
  df1 <- df1[,order(mns,decreasing=T)]
  barplot(colMeans(df1), las=2, cex.names=.85, ylab="Update time (ms)", ylim=c(0,(max(mns)*1.1)))
  dev.off()
  
  pdf("./inst/st_clusout_pt.pdf",  width=3, height=4)
  par(mar=c(6, 4, 1, 1) + 0.1)
  df1 <- data.frame('SHC sequential'=res1$processTime,'SHC index 1'=res2$processTime,
                    'SHC index 2'=res3$processTime,'SHC index 3'=res4$processTime,
                    'SHC balanced'=res5$processTime,check.names=F)
  mns <- colMeans(df1)
  df1 <- df1[,order(mns,decreasing=T)]
  barplot(colMeans(df1), las=2, cex.names=.85, ylab="Processing time (ms)", ylim=c(0,(max(mns)*1.1)))
  dev.off()
}

SHC_TestCase_ClustersAndOutliers_SigmaIndex_Theorems <- function() {
  set.seed(0)
  ds <- DSD_SHCGaussianGenerator(clusters = 5, outliers = 0, initPopulation = 1000, 
                                 cacheAndRepeat = TRUE, maxDimensionValues = c(50,50),
                                 generateVariableClusterPopulations = TRUE, theta=10,
                                 maxVariance = 1.5, checkOverlapping_MD = T)
  
  print("#SHC, N=2")
  c1 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 2)
  reset_stream(ds)
  update(c1, ds, n=1000, use_outliers=T)
  
  hist <- getHistogram(c1)
  pdf("./inst/st_theorem_shc1_histogram.pdf",  width=3.2, height=4)
  par(mar=c(5, 4, 1, 1) + 0.1)
  barplot(unname(hist[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)", 
          space=0, beside=T)
  axis(side=1,at=c(0,20,40,60,80,100))
  dev.off()
  
  print("#SHC, N=3")
  c2 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 3)
  reset_stream(ds)
  update(c2, ds, n=1000, use_outliers=T)
  
  hist <- getHistogram(c2)
  pdf("./inst/st_theorem_shc2_histogram.pdf",  width=3.2, height=4)
  par(mar=c(5, 4, 1, 1) + 0.1)
  barplot(unname(hist[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)", 
          space=0, beside=T)
  axis(side=1,at=c(0,20,40,60,80,100))
  dev.off()
  
  print("#SHC, N=4")
  c3 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 4)
  reset_stream(ds)
  update(c3, ds, n=1000, use_outliers=T)
  
  hist <- getHistogram(c3)
  pdf("./inst/st_theorem_shc3_histogram.pdf",  width=3.2, height=4)
  par(mar=c(5, 4, 1, 1) + 0.1)
  barplot(unname(hist[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)",
          space=0, beside=T)
  axis(side=1,at=c(0,20,40,60,80,100))
  dev.off()
}

SHC_TestCase_Sensors_SigmaIndex <- function(endingRound = 2000) {
  prepareSensorDataset()
  ds <- DSD_ReadCSV("./inst/datasets/sensors_final.csv", take=c(2:6), class=2, header=T)
  
  print("Sensors SHC-sequential")
  c1 <- DSC_SHC.man(4, 0.789, 0.617, compAssimilationCheckCounter = 1000, cbNLimit = 20, 
                    driftRemoveCompSizeRatio = 0.166, driftCheckingSizeRatio = 0.732, 
                    driftMovementMDThetaRatio = 0.258, decaySpeed = 31, recStats = T, 
                    compFormingMinVVRatio = 0.476)
  reset_stream(ds)
  shc_res1 <- evaluate_cluster_with_callbacks(c1, ds, n = endingRound*1000, type="macro",
                                              measure = c("cRand","queryTime", "updateTime", "processTime",
                                                          "nodeCount", "computationCostReduction"),
                                              callbacks = list(shc=SHCEvalCallback()), horizon = 1000,
                                              verbose = T, use_outliers=T)
  saveRDS(shc_res1, file="./inst/sensors_shc1.RDS")
  
  print("Sensors SHC-index 1")
  c2 <- DSC_SHC.man(4, 0.789, 0.617, compAssimilationCheckCounter = 1000, cbNLimit = 20, 
                    driftRemoveCompSizeRatio = 0.166, driftCheckingSizeRatio = 0.732, 
                    driftMovementMDThetaRatio = 0.258, decaySpeed = 31, recStats = T, 
                    compFormingMinVVRatio = 0.476, sigmaIndex = T, sigmaIndexNeighborhood = 3)
  reset_stream(ds)
  shc_res2 <- evaluate_cluster_with_callbacks(c2, ds, n = endingRound*1000, type="macro",
                                              measure = c("cRand","queryTime", "updateTime", "processTime",
                                                          "nodeCount", "computationCostReduction"),
                                              callbacks = list(shc=SHCEvalCallback()), horizon = 1000,
                                              verbose = T, use_outliers=T)
  saveRDS(shc_res2, file="./inst/sensors_shc2.RDS")
  
  print("Sensors SHC-index 2")
  c3 <- DSC_SHC.man(4, 0.789, 0.617, compAssimilationCheckCounter = 1000, cbNLimit = 20, 
                    driftRemoveCompSizeRatio = 0.166, driftCheckingSizeRatio = 0.732, 
                    driftMovementMDThetaRatio = 0.258, decaySpeed = 31, recStats = T, 
                    compFormingMinVVRatio = 0.476, sigmaIndex = T, sigmaIndexNeighborhood = 6)
  reset_stream(ds)
  shc_res3 <- evaluate_cluster_with_callbacks(c3, ds, n = endingRound*1000, type="macro",
                                              measure = c("cRand","queryTime", "updateTime", "processTime",
                                                          "nodeCount", "computationCostReduction"),
                                              callbacks = list(shc=SHCEvalCallback()), horizon = 1000,
                                              verbose = T, use_outliers=T)
  saveRDS(shc_res3, file="./inst/sensors_shc3.RDS")
}

SHC_TestCase_Sensors_SigmaIndex_Display <- function() {
  if(!file.exists("./inst/sensors_shc1.RDS") || !file.exists("./inst/sensors_shc2.RDS") ||
     !file.exists("./inst/sensors_shc3.RDS"))
    stop("One of the sensor resulting files is missing. Run SHC_TestCase_Sensors_SigmaIndex() again.")
  shc_res1 <- readRDS("./inst/sensors_shc1.RDS")
  shc_res2 <- readRDS("./inst/sensors_shc2.RDS")
  shc_res3 <- readRDS("./inst/sensors_shc3.RDS")
  
  pts <- shc_res1[,"points"]/1000
  shc_res1_qt <- aggregate(shc_res1[,"queryTime"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res2_qt <- aggregate(shc_res2[,"queryTime"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res3_qt <- aggregate(shc_res3[,"queryTime"], by = list(pts %/% 10 *10), FUN=mean)
  
  pdf("./inst/st_sensors_qt.pdf", width=9, height=4)
  par(mar=c(6, 4, 2, 2) + 0.1)
  plot(shc_res1_qt, type="l", 
       xlab="Position in Stream (1000s)", ylab="Query time (ms)")
  lines(shc_res2_qt, type="l", lty=2)
  lines(shc_res3_qt, type="l", lty=3, lwd=2)
  legend(x="topleft", legend=c("SHC sequential","SHC index 1","SHC index 2"),
         lty=c(1,2,3), lwd=c(1,1,2), bty="n")
  dev.off()
  
  shc_res1_ut <- aggregate(shc_res1[,"updateTime"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res2_ut <- aggregate(shc_res2[,"updateTime"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res3_ut <- aggregate(shc_res3[,"updateTime"], by = list(pts %/% 10 *10), FUN=mean)
  
  pdf("./inst/st_sensors_ut.pdf", width=9, height=4)
  par(mar=c(6, 4, 2, 2) + 0.1)
  plot(shc_res2_ut, type="l", 
       xlab="Position in Stream (1000s)", ylab="Update time (ms)", ylim=c(0,max(max(shc_res2_ut[2]),max(shc_res3_ut[2]))))
  lines(shc_res3_ut, type="l", lty=2)
  legend(x="topleft", legend=c("SHC index 1","SHC index 2"),
         lty=c(1,2), lwd=c(1,1), bty="n")
  dev.off()
  
  shc_res1_pt <- aggregate(shc_res1[,"processTime"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res2_pt <- aggregate(shc_res2[,"processTime"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res3_pt <- aggregate(shc_res3[,"processTime"], by = list(pts %/% 10 *10), FUN=mean)
  
  pdf("./inst/st_sensors_pt.pdf", width=9, height=4)
  par(mar=c(6, 4, 2, 2) + 0.1)
  plot(shc_res1_pt, type="l", 
       xlab="Position in Stream (1000s)", ylab="Total processing time (ms)")
  lines(shc_res2_pt, type="l", lty=2)
  lines(shc_res3_pt, type="l", lty=3, lwd=2)
  legend(x="topleft", legend=c("SHC sequential","SHC index 1","SHC index 2"),
         lty=c(1,2,3), lwd=c(1,1,2), bty="n")
  dev.off()
  
  shc_res1_nc <- aggregate(shc_res1[,"nodeCount"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res2_nc <- aggregate(shc_res2[,"nodeCount"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res3_nc <- aggregate(shc_res3[,"nodeCount"], by = list(pts %/% 10 *10), FUN=mean)
  
  pdf("./inst/st_sensors_nc.pdf", width=9, height=4)
  par(mar=c(6, 4, 2, 2) + 0.1)
  plot(shc_res1_nc, type="l",
       xlab="Position in Stream (1000s)", ylab="Nodes visited/calculated")
  lines(shc_res2_nc, type="l", lty=2)
  lines(shc_res3_nc, type="l", lty=3, lwd=2)
  legend(x="topleft", legend=c("SHC sequential","SHC index 1","SHC index 2"),
         lty=c(1,2,3), lwd=c(1,1,2), bty="n")
  dev.off()
  
  shc_res2_nc <- aggregate(shc_res2[,"computationCostReduction"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res3_nc <- aggregate(shc_res3[,"computationCostReduction"], by = list(pts %/% 10 *10), FUN=mean)
  
  pdf("./inst/st_sensors_ccr.pdf", width=9, height=4)
  par(mar=c(6, 4, 2, 2) + 0.1)
  plot(shc_res2_nc, type="l",
       xlab="Position in Stream (1000s)", ylab="Comp. cost reduction (%)", ylim=c(0,100))
  lines(shc_res3_nc, type="l", lty=2)
  legend(x="bottomright", legend=c("SHC index 1","SHC index 2"),
         lty=c(1,2), lwd=c(1,1), bty="n")
  dev.off()
  
  shc_res1_a <- aggregate(shc_res1[,"cRand"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res2_a <- aggregate(shc_res2[,"cRand"], by = list(pts %/% 10 *10), FUN=mean)
  shc_res3_a <- aggregate(shc_res3[,"cRand"], by = list(pts %/% 10 *10), FUN=mean)
  
  pdf("./inst/st_sensors_crand.pdf",  width=4, height=4)
  par(mar=c(6, 4, 2, 2) + 0.1)
  df1 <- data.frame('SHC sequential'=shc_res1[,"cRand"],
                    'SHC index 1'=shc_res2[,"cRand"],
                    'SHC index 2'=shc_res3[,"cRand"],
                    check.names=F)
  mns <- colMeans(df1)
  df1 <- df1[,order(mns,decreasing=T)]
  barplot(colMeans(df1), las=2, cex.names=.85, ylab="Avg. Corr. Rand", ylim=c(0,1))
  dev.off()
}

prepareSensorDataset <- function() {
  # since we cannot pack big datasets into R packages, downloading it directly is the only
  # viable solution
  if(file.exists("./inst/datasets/sensors_final.csv")) return()
  if(!dir.exists("./inst/datasets")) dir.create("./inst/datasets")
  download.file("http://db.csail.mit.edu/labdata/data.txt.gz","./inst/datasets/sensors.csv.gz")
  loadNamespace("R.utils")
  message("unpacking...")
  R.utils::gunzip("./inst/datasets/sensors.csv.gz")
  message("reading from CSV...")
  r <- read.csv("./inst/datasets/sensors.csv",sep=" ",header=F)
  message("cleaning the dataset...")
  colnames(r) <- c("date","time","epoch","sensor","temperature","humidity","light","voltage")
  r <- transform(r, datetime=as.POSIXct(paste(date, time), format="%Y-%m-%d %H:%M:%S"))
  r <- r[,-which(colnames(r) %in% c("date","time","epoch"))]
  r <- r[order(r$datetime),]
  r <- r[!is.na(r$sensor),]
  r[is.na(r$light),"light"] <- 0
  r[is.na(r$temperature),"temperature"] <- 0
  r[is.na(r$humidity),"humidity"] <- 0
  r[is.na(r$voltage),"voltage"] <- 0
  message("deleting the original file...")
  file.remove("./inst/datasets/sensors.csv")
  message("writing the final dataset to CSV...")
  write.csv(r,"./inst/datasets/sensors_final.csv")
}