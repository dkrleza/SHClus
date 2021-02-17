### R code from vignette source 'SigmaIndex.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("SHClus")
library("stream")


###################################################
### code chunk number 2: SigmaIndex.Rnw:697-698
###################################################
library(SHClus)


###################################################
### code chunk number 3: SigmaIndex.Rnw:700-701
###################################################
si <- SigmaIndex(theta = 3.2, neighborhood = 9.6)


###################################################
### code chunk number 4: SigmaIndex.Rnw:705-712
###################################################
covariance1 <- matrix(data=c(0.8, 0, 0, 0.8), nrow=2, ncol=2)
covariance1
addPopulation(si, "1", c(0.5,0.5), covariance1, 800)
covariance2 <- matrix(data=c(0.5, 0, 0, 0.5), nrow=2, ncol=2)
addPopulation(si, "2", c(7,0.5), covariance2, 500)
covariance3 <- matrix(data=c(0.4, 0, 0, 0.4), nrow=2, ncol=2)
addPopulation(si, "3", c(3.5,10.0), covariance3, 400)


###################################################
### code chunk number 5: SigmaIndex.Rnw:715-718
###################################################
out_covariance <- matrix(data=c(0.3, 0, 0, 0.3), nrow=2, ncol=2)
addPopulation(si, "4", c(8.0,6.2), out_covariance, 1)
addPopulation(si, "5", c(0.5,15.0), out_covariance, 1)


###################################################
### code chunk number 6: SigmaIndex.Rnw:727-728
###################################################
getTotalPopulations(si)


###################################################
### code chunk number 7: SigmaIndex.Rnw:731-736
###################################################
pop <- getPopulations(si)
names(pop)
outlier <- pop$"4"
outlier$mean
outlier$icovariance


###################################################
### code chunk number 8: SigmaIndex.Rnw:741-743
###################################################
query_data <- data.frame(X=c(0.37,6.5,8.05),Y=c(0.505,0.4,6.3))
query_data


###################################################
### code chunk number 9: SigmaIndex.Rnw:746-747
###################################################
res <- queryDataPoints(si, query_data)


###################################################
### code chunk number 10: SigmaIndex.Rnw:750-751
###################################################
unlist(res[[1]]$outcome)


###################################################
### code chunk number 11: SigmaIndex.Rnw:754-755
###################################################
unlist(res[[1]]$neighborhood)


###################################################
### code chunk number 12: SigmaIndex.Rnw:759-762
###################################################
res <- queryDataPoints(si, c(3.75,0.5))
unlist(res[[1]]$outcome)
unlist(res[[1]]$neighborhood)


###################################################
### code chunk number 13: SigmaIndex.Rnw:767-768
###################################################
resetStatistics(si)


###################################################
### code chunk number 14: SigmaIndex.Rnw:771-773
###################################################
res <- queryDataPoints(si, c(0.9,0.7))
unlist(res)


###################################################
### code chunk number 15: SigmaIndex.Rnw:776-777
###################################################
unlist(getStatistics(si))


###################################################
### code chunk number 16: SigmaIndex.Rnw:780-792
###################################################
si <- SigmaIndex(theta = 3.2, neighborhood = 9.6)
covariance1 <- matrix(data=c(0.8, 0, 0, 0.8), nrow=2, ncol=2)
addPopulation(si, "1", c(0.5,0.5), covariance1, 800)
covariance2 <- matrix(data=c(0.5, 0, 0, 0.5), nrow=2, ncol=2)
addPopulation(si, "2", c(7,0.5), covariance2, 500)
covariance3 <- matrix(data=c(0.4, 0, 0, 0.4), nrow=2, ncol=2)
addPopulation(si, "3", c(3.5,10.0), covariance3, 400)
out_covariance <- matrix(data=c(0.3, 0, 0, 0.3), nrow=2, ncol=2)
addPopulation(si, "4", c(6.5,15.0), out_covariance, 1)
addPopulation(si, "5", c(0.5,15.0), out_covariance, 1)
res <- queryDataPoints(si, c(0.9,0.7))
unlist(getStatistics(si))


###################################################
### code chunk number 17: SigmaIndex.Rnw:798-800
###################################################
hist <- getHistogram(si)
hist


###################################################
### code chunk number 18: si1
###################################################
par(mar=c(5, 4, 1, 1) + 0.1)
barplot(unname(hist[1,]), cex.names=.75, ylab="Density", xlab="Comp. cost reduction (%)", space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100))


###################################################
### code chunk number 19: SigmaIndex.Rnw:814-815
###################################################
sum(hist)


###################################################
### code chunk number 20: SigmaIndex.Rnw:826-828
###################################################
new_observ <- data.frame(X=c(0.6,6.9),Y=c(0.4,0.55))
new_observ


###################################################
### code chunk number 21: SigmaIndex.Rnw:831-833
###################################################
ids <- c("1","3")
addDataPoints(si, ids, new_observ)


###################################################
### code chunk number 22: SigmaIndex.Rnw:843-845
###################################################
new_observ <- data.frame(X=c(0.62,6.91),Y=c(0.41,0.551))
new_observ


###################################################
### code chunk number 23: SigmaIndex.Rnw:848-852
###################################################
for(r in 1:nrow(new_observ)) {
query_res <- queryDataPoints(si, new_observ[r,])
addDataPointsInc(si, new_observ[r,], query_res)
}


###################################################
### code chunk number 24: SigmaIndex.Rnw:860-864
###################################################
library(stream)
library(SHClus)
set.seed(1000)
n <- 20000


###################################################
### code chunk number 25: SigmaIndex.Rnw:866-869
###################################################
dsd <- DSD_Gaussians(k=50,outliers=50,separation_type="Mahalanobis",
  space_limit=c(0,150),separation=4,variance_limit=8,
  outlier_options=list(outlier_horizon=20000))


###################################################
### code chunk number 26: si2
###################################################
par(mar=c(3, 3, 0.8, 0.8), mgp=c(1.8,0.7,0), cex.axis=1.5, cex.lab=1.5)
plot(dsd, n)


###################################################
### code chunk number 27: SigmaIndex.Rnw:877-879
###################################################
si <- convertFromDSD(dsd, total_elements = 20000, theta = 3.2, 
neighborhood = 9.6)


###################################################
### code chunk number 28: SigmaIndex.Rnw:882-884
###################################################
options("scipen"=100, "digits"=4)
reset_stream(dsd)


###################################################
### code chunk number 29: SigmaIndex.Rnw:886-890
###################################################
resetStatistics(si)
res1 <- queryDataPoints(si, get_points(dsd, 20000))
stats1 <- getStatistics(si)
unlist(stats1)


###################################################
### code chunk number 30: SigmaIndex.Rnw:893-906
###################################################
hist1 <- getHistogram(si)
reset_stream(dsd)
si <- convertFromDSD(dsd, total_elements = n, theta = 3.2, 
neighborhood = 6.4)
res2 <- queryDataPoints(si, get_points(dsd, n))
stats2 <- getStatistics(si)
hist2 <- getHistogram(si)
reset_stream(dsd)
si <- convertFromDSD(dsd, total_elements = n, theta = 3.2, 
neighborhood = 12.9)
res3 <- queryDataPoints(si, get_points(dsd, n))
stats3 <- getStatistics(si)
hist3 <- getHistogram(si)


###################################################
### code chunk number 31: si3
###################################################
par(mar=c(6, 4, 1, 1) + 0.1)
df1 <- data.frame('6.4'=stats2$computationCostReduction,
'9.6'=stats1$computationCostReduction,
'12.9'=stats3$computationCostReduction,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns)]
barplot(colMeans(df1), las=2, ylab="Comp. cost reduction (%)")


###################################################
### code chunk number 32: si4
###################################################
par(mar=c(5, 4, 1, 1) + 0.1)
barplot(unname(hist2[1,]), ylab="Density", xlab="Comp. cost reduction (%)", space=0, beside=T, ylim=c(0,(max(mns)*1.1)))
axis(side=1,at=c(0,20,40,60,80,100))


###################################################
### code chunk number 33: si5
###################################################
par(mar=c(5, 4, 1, 1) + 0.1)
barplot(unname(hist1[1,]), ylab="Density", xlab="Comp. cost reduction (%)", space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100))


###################################################
### code chunk number 34: si6
###################################################
par(mar=c(5, 4, 1, 1) + 0.1)
barplot(unname(hist3[1,]), ylab="Density", xlab="Comp. cost reduction (%)", space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100))


###################################################
### code chunk number 35: SigmaIndex.Rnw:970-971
###################################################
set.seed(1000)


###################################################
### code chunk number 36: SigmaIndex.Rnw:973-976
###################################################
dsd <- DSD_Gaussians(k=50,outliers=50,separation_type="Mahalanobis",
                     separation=4,space_limit=c(0,150),variance_limit=8,
                     outlier_options=list(outlier_horizon=20000))


###################################################
### code chunk number 37: SigmaIndex.Rnw:980-983
###################################################
shc_c <- DSC_SHC.man(2, 3.5, 0.9, cbVarianceLimit = 0, cbNLimit = 0, 
decaySpeed = 0, sharedAgglomerationThreshold = 100, sigmaIndex = T, 
sigmaIndexNeighborhood = 3)


###################################################
### code chunk number 38: SigmaIndex.Rnw:987-989
###################################################
evaluate(shc_c, dsd, n=20000, type="macro", measure=c("crand",
"outlierjaccard"), single_pass_update=T)


###################################################
### code chunk number 39: SigmaIndex.Rnw:993-994
###################################################
shc_c$RObj$getComputationCostReduction()


###################################################
### code chunk number 40: SigmaIndex.Rnw:997-1017
###################################################
n <- 20000
reset_stream(dsd)
c1 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = FALSE)
c1$RObj$setPseudoOfflineCounter(500)
res1 <- evaluate_with_callbacks(c1, dsd, n=n, measure = c("cRand", "queryTime", "updateTime", "processTime","nodeCount",  "computationCostReduction","OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()),single_pass_update=T, use_outliers=T)

reset_stream(dsd)
c2 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 2)
c2$RObj$setPseudoOfflineCounter(500)
res2 <- evaluate_with_callbacks(c2, dsd, n=n, measure = c("cRand", "queryTime", "updateTime", "processTime","nodeCount", "computationCostReduction", "OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()), single_pass_update=T, use_outliers=T)

reset_stream(dsd)
c3 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 3)
c3$RObj$setPseudoOfflineCounter(500)
res3 <- evaluate_with_callbacks(c3, dsd, n=n, measure = c("cRand", "queryTime","updateTime", "processTime","nodeCount", "computationCostReduction", "OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()), single_pass_update=T, use_outliers=T)

reset_stream(dsd)
c4 <- DSC_SHC.behavioral(2, AgglomerationType$NormalAgglomeration, DriftType$NoDrift, 0, sigmaIndex = TRUE, sigmaIndexNeighborhood = 4)
c4$RObj$setPseudoOfflineCounter(500)
res4 <- evaluate_with_callbacks(c4, dsd, n=n, measure = c("cRand", "queryTime","updateTime", "processTime","nodeCount", "computationCostReduction", "OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()), single_pass_update=T, use_outliers=T)


###################################################
### code chunk number 41: si7
###################################################
hist2 <- getHistogram(c2)
par(mar=c(4, 3, 0.5, 0.5), mgp=c(1.5,0.5,0))
barplot(unname(hist2[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 42: si8
###################################################
hist3 <- getHistogram(c3)
par(mar=c(4, 3, 0.5, 0.5), mgp=c(1.5,0.5,0))
barplot(unname(hist3[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 43: si9
###################################################
hist4 <- getHistogram(c4)
par(mar=c(4, 3, 0.5, 0.5), mgp=c(1.5,0.5,0))
barplot(unname(hist4[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 44: si-shc-crand
###################################################
par(mar=c(4, 3, 0.2, 0.2), mgp=c(1.5,0.7,0))
df1 <- data.frame('SHC(seq.)'=res1$cRand,'SHC(nm=2)'=res2$cRand,
'SHC(nm=3)'=res3$cRand,'SHC(nm=4)'=res4$cRand,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns,decreasing=T)]
boxplot(df1, las=2, ylab="Corrected Rand", ylim=c(0,1), cex.names=.6, cex.axis=.6, cex.lab=.8)


###################################################
### code chunk number 45: si-shc-qt
###################################################
par(mar=c(4, 3, 0.2, 0.2))
df1 <- data.frame('SHC(seq.)'=res1$queryTime,'SHC(nm=2)'=res2$queryTime,
'SHC(nm=3)'=res3$queryTime,'SHC(nm=4)'=res4$queryTime,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, ylab="Query time (ms)", 
        ylim=c(0,(max(mns)*1.1)), mgp=c(1.9,0.7,0))


###################################################
### code chunk number 46: si-shc-nc
###################################################
par(mar=c(4, 3.5, 0.2, 0.2))
df1 <- data.frame('SHC(seq.)'=res1$nodeCount,'SHC(nm=2)'=res2$nodeCount,
'SHC(nm=3)'=res3$nodeCount,'SHC(nm=4)'=res4$nodeCount,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, mgp=c(2.6,0.7,0), ylab="Nodes visited/calculated", ylim=c(0,(max(mns)*1.1)))


###################################################
### code chunk number 47: si-shc-ccr
###################################################
par(mar=c(4, 3, 0.2, 0.2))
df1 <- data.frame('SHC(seq.)'=res1$computationCostReduction,'SHC(nm=2)'=res2$computationCostReduction,
'SHC(nm=3)'=res3$computationCostReduction,'SHC(nm=4)'=res4$computationCostReduction,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, mgp=c(1.5,0.7,0), ylab="Comp. cost reduction (%)", ylim=c(0,(max(mns)*1.1)))


###################################################
### code chunk number 48: si-shc-ut
###################################################
par(mar=c(4, 2.8, 0.2, 0.2))
df1 <- data.frame('SHC(seq.)'=res1$updateTime,'SHC(nm=2)'=res2$updateTime,
'SHC(nm=3)'=res3$updateTime,'SHC(nm=4)'=res4$updateTime,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns,decreasing=T)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, mgp=c(1.5,0.7,0), ylab="Update time (ms)", ylim=c(0,(max(mns)*1.1)))


###################################################
### code chunk number 49: si-shc-pt
###################################################
par(mar=c(4, 3, 0.2, 0.2))
df1 <- data.frame('SHC(seq.)'=res1$processTime,'SHC(nm=2)'=res2$processTime,
'SHC(nm=3)'=res3$processTime,'SHC(nm=4)'=res4$processTime,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns,decreasing=T)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, mgp=c(1.8,0.7,0), ylab="Processing time (ms)", ylim=c(0,(max(mns)*1.1)))


###################################################
### code chunk number 50: SigmaIndex.Rnw:1132-1133
###################################################
set.seed(1000)


###################################################
### code chunk number 51: SigmaIndex.Rnw:1135-1139
###################################################
dsd <- DSD_Gaussians(k=50,outliers=50,separation_type="Mahalanobis",
    separation=2,space_limit=c(0,90),variance_limit=8,
    outlier_options=list(outlier_horizon=20000,
                         outlier_virtual_variance=0.3))


###################################################
### code chunk number 52: si-dsdplot-notseparated
###################################################
par(mar=c(3, 3, 0.8, 0.8), mgp=c(1.8,0.7,0), cex.axis=1.5, cex.lab=1.5)
plot(dsd, 20000)
reset_stream(dsd)


###################################################
### code chunk number 53: SigmaIndex.Rnw:1147-1150
###################################################
shc_c <- DSC_SHC.man(2, 3.2, 0.3, cbVarianceLimit = 0, cbNLimit = 0, 
decaySpeed = 0, sharedAgglomerationThreshold = 700, sigmaIndex = T, 
sigmaIndexNeighborhood = 3)


###################################################
### code chunk number 54: SigmaIndex.Rnw:1152-1153
###################################################
shc_c$RObj$setPseudoOfflineCounter(500)


###################################################
### code chunk number 55: SigmaIndex.Rnw:1156-1158
###################################################
evaluate(shc_c, dsd, n=20000, type="macro", measure=c("crand",
"outlierjaccard"), single_pass_update=T)


###################################################
### code chunk number 56: SigmaIndex.Rnw:1161-1181
###################################################
n <- 20000
reset_stream(dsd)
c1 <- DSC_SHC.man(2, 3.2, 0.3, cbVarianceLimit = 0, cbNLimit = 0, decaySpeed = 0, sharedAgglomerationThreshold = 700, sigmaIndex = F)
c1$RObj$setPseudoOfflineCounter(500)
res1 <- evaluate_with_callbacks(c1, dsd, n=n, measure = c("cRand", "queryTime", "updateTime", "processTime","nodeCount",  "computationCostReduction","OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()),single_pass_update=T, use_outliers=T)

reset_stream(dsd)
c2 <- DSC_SHC.man(2, 3.2, 0.3, cbVarianceLimit = 0, cbNLimit = 0, decaySpeed = 0, sharedAgglomerationThreshold = 700, sigmaIndex = T, sigmaIndexNeighborhood = 2)
c2$RObj$setPseudoOfflineCounter(500)
res2 <- evaluate_with_callbacks(c2, dsd, n=n, measure = c("cRand", "queryTime", "updateTime", "processTime","nodeCount", "computationCostReduction", "OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()), single_pass_update=T, use_outliers=T)

reset_stream(dsd)
c3 <- DSC_SHC.man(2, 3.2, 0.3, cbVarianceLimit = 0, cbNLimit = 0, decaySpeed = 0, sharedAgglomerationThreshold = 700, sigmaIndex = T, sigmaIndexNeighborhood = 3)
c3$RObj$setPseudoOfflineCounter(500)
res3 <- evaluate_with_callbacks(c3, dsd, n=n, measure = c("cRand", "queryTime","updateTime", "processTime","nodeCount", "computationCostReduction", "OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()), single_pass_update=T, use_outliers=T)

reset_stream(dsd)
c4 <- DSC_SHC.man(2, 3.2, 0.3, cbVarianceLimit = 0, cbNLimit = 0, decaySpeed = 0, sharedAgglomerationThreshold = 700, sigmaIndex = T, sigmaIndexNeighborhood = 4)
c4$RObj$setPseudoOfflineCounter(500)
res4 <- evaluate_with_callbacks(c4, dsd, n=n, measure = c("cRand", "queryTime","updateTime", "processTime","nodeCount", "computationCostReduction", "OutlierJaccard"), type="macro", callbacks=list(shc=SHCEvalCallback()), single_pass_update=T, use_outliers=T)


###################################################
### code chunk number 57: si10
###################################################
hist2 <- getHistogram(c2)
par(mar=c(4, 3, 0.5, 0.5), mgp=c(1.5,0.5,0))
barplot(unname(hist2[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 58: si11
###################################################
hist3 <- getHistogram(c3)
par(mar=c(4, 3, 0.5, 0.5), mgp=c(1.5,0.5,0))
barplot(unname(hist3[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 59: si12
###################################################
hist4 <- getHistogram(c4)
par(mar=c(4, 3, 0.5, 0.5), mgp=c(1.5,0.5,0))
barplot(unname(hist4[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 60: si-shc-crand2
###################################################
par(mar=c(4, 3, 0.2, 0.2), mgp=c(1.5,0.7,0))
df1 <- data.frame('SHC(seq.)'=res1$cRand,'SHC(nm=2)'=res2$cRand,'SHC(nm=3)'=res3$cRand,'SHC(nm=4)'=res4$cRand,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns,decreasing=T)]
boxplot(df1, las=2, ylab="Corrected Rand", ylim=c(0,1), cex.names=.6, cex.axis=.6, cex.lab=.8)


###################################################
### code chunk number 61: si-shc-oji2
###################################################
par(mar=c(4, 3, 0.2, 0.2), mgp=c(1.5,0.7,0))
df1 <- data.frame('SHC(seq.)'=res1$OutlierJaccard,'SHC(nm=2)'=res2$OutlierJaccard,'SHC(nm=3)'=res3$OutlierJaccard,'SHC(nm=4)'=res4$OutlierJaccard,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns,decreasing=T)]
boxplot(df1, las=2, ylab="Outlier Jaccard", ylim=c(0,1), cex.names=.6, cex.axis=.6, cex.lab=.8)


###################################################
### code chunk number 62: si-shc-qt2
###################################################
par(mar=c(4, 3, 0.2, 0.2),mgp=c(1.9,0.7,0))
df1 <- data.frame('SHC(seq.)'=res1$queryTime,'SHC(nm=2)'=res2$queryTime,'SHC(nm=3)'=res3$queryTime,'SHC(nm=4)'=res4$queryTime,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, ylab="Query time (ms)", ylim=c(0,(max(mns)*1.1)))


###################################################
### code chunk number 63: si-shc-ccr2
###################################################
par(mar=c(4, 3, 0.2, 0.2),mgp=c(1.5,0.7,0))
df1 <- data.frame('SHC(seq.)'=res1$computationCostReduction,'SHC(nm=2)'=res2$computationCostReduction,'SHC(nm=3)'=res3$computationCostReduction,'SHC(nm=4)'=res4$computationCostReduction,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, ylab="Comp. cost reduction (%)", ylim=c(0,(max(mns)*1.1)))


###################################################
### code chunk number 64: si-shc-pt2
###################################################
par(mar=c(4, 3, 0.2, 0.2),mgp=c(1.8,0.7,0))
df1 <- data.frame('SHC(seq.)'=res1$processTime,'SHC(nm=2)'=res2$processTime,
'SHC(nm=3)'=res3$processTime,'SHC(nm=4)'=res4$processTime,check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns,decreasing=T)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, ylab="Processing time (ms)", ylim=c(0,(max(mns)*1.1)))


###################################################
### code chunk number 65: SigmaIndex.Rnw:1285-1290
###################################################
shc_res1 <- readRDS("sensors_shc1.RDS")
shc_res2 <- readRDS("sensors_shc2.RDS")
shc_res3 <- readRDS("sensors_shc3.RDS")

pts <- shc_res1[,"points"]/1000


###################################################
### code chunk number 66: si-shc-sensors-qt
###################################################
shc_res1_qt <- aggregate(shc_res1[,"queryTime"], by = list(pts %/% 10 *10), FUN=mean)
shc_res2_qt <- aggregate(shc_res2[,"queryTime"], by = list(pts %/% 10 *10), FUN=mean)
shc_res3_qt <- aggregate(shc_res3[,"queryTime"], by = list(pts %/% 10 *10), FUN=mean)

par(mar=c(4,3,0.5,0.5), mgp=c(1.8,0.7,0), cex.lab=1.5)
plot(shc_res1_qt, type="l", lty=3, xlab="Position in Stream (1000s)", ylab="Query time (ms)")
lines(shc_res2_qt, type="l", lty=2)
lines(shc_res3_qt, type="l", lty=1, lwd=2)
legend(x="topleft", legend=c("SHC(sequential)","SHC(nm=3)","SHC(nm=6)"), lty=c(3,2,1), lwd=c(1,1,2), bty="n")


###################################################
### code chunk number 67: si-shc-sensors-ut
###################################################
shc_res1_ut <- aggregate(shc_res1[,"updateTime"], by = list(pts %/% 10 *10), FUN=mean)
shc_res2_ut <- aggregate(shc_res2[,"updateTime"], by = list(pts %/% 10 *10), FUN=mean)
shc_res3_ut <- aggregate(shc_res3[,"updateTime"], by = list(pts %/% 10 *10), FUN=mean)

par(mar=c(4,3,0.5,0.5), mgp=c(1.8,0.7,0), cex.lab=1.5)
plot(shc_res2_ut, type="l", lty=2, xlab="Position in Stream (1000s)", ylab="Update time (ms)", ylim=c(0,max(max(shc_res2_ut[2]),max(shc_res3_ut[2]))))
lines(shc_res3_ut, type="l", lty=3)
legend(x="topleft", legend=c("SHC(nm=3)","SHC(nm=6)"),
lty=c(2,3), lwd=c(1,1), bty="n")


###################################################
### code chunk number 68: si-shc-sensors-pt
###################################################
shc_res1_pt <- aggregate(shc_res1[,"processTime"], by = list(pts %/% 10 *10), FUN=mean)
shc_res2_pt <- aggregate(shc_res2[,"processTime"], by = list(pts %/% 10 *10), FUN=mean)
shc_res3_pt <- aggregate(shc_res3[,"processTime"], by = list(pts %/% 10 *10), FUN=mean)

par(mar=c(4,3,0.5,0.8), mgp=c(1.8,0.7,0), cex.lab=1.5)
plot(shc_res1_pt, type="l", lty=3, xlab="Position in Stream (1000s)", ylab="Processing time (ms)")
lines(shc_res2_pt, type="l", lty=2)
lines(shc_res3_pt, type="l", lty=1, lwd=2)
legend(x="topleft", legend=c("SHC(sequential)","SHC(nm=3)","SHC(nm=6)"),
lty=c(3,2,1), lwd=c(1,1,2), bty="n")


###################################################
### code chunk number 69: si-shc-sensors-nc
###################################################
shc_res1_nc <- aggregate(shc_res1[,"nodeCount"], by = list(pts %/% 10 *10), FUN=mean)
shc_res2_nc <- aggregate(shc_res2[,"nodeCount"], by = list(pts %/% 10 *10), FUN=mean)
shc_res3_nc <- aggregate(shc_res3[,"nodeCount"], by = list(pts %/% 10 *10), FUN=mean)

par(mar=c(4,3,0.5,0.5), mgp=c(1.8,0.7,0), cex.lab=1.5)
plot(shc_res1_nc, type="l", lty=3, xlab="Position in Stream (1000s)", ylab="Nodes visited/calculated")
lines(shc_res2_nc, type="l", lty=2)
lines(shc_res3_nc, type="l", lty=1, lwd=2)
legend(x="topleft", legend=c("SHC(sequential)","SHC(nm=3)","SHC(nm=6)"),
lty=c(3,2,1), lwd=c(1,1,2), bty="n")


###################################################
### code chunk number 70: si-shc-sensors-ccr
###################################################
shc_res2_ccr <- aggregate(shc_res2[,"computationCostReduction"], by = list(pts %/% 10 *10), FUN=mean)
shc_res3_ccr <- aggregate(shc_res3[,"computationCostReduction"], by = list(pts %/% 10 *10), FUN=mean)

par(mar=c(4,3,0.5,0.5), mgp=c(1.8,0.7,0), cex.lab=1.5)
plot(shc_res2_ccr, type="l", lty=2, xlab="Position in Stream (1000s)", ylab="Comp. cost reduction (%)", ylim=c(0,100))
lines(shc_res3_ccr, type="l", lty=3)
legend(x="bottomright", legend=c("SHC(nm=3)","SHC(nm=6)"),
lty=c(2,3), lwd=c(1,1), bty="n")


###################################################
### code chunk number 71: si-shc-sensors-crand
###################################################
shc_res1_a <- aggregate(shc_res1[,"cRand"], by = list(pts %/% 10 *10), FUN=mean)
shc_res2_a <- aggregate(shc_res2[,"cRand"], by = list(pts %/% 10 *10), FUN=mean)
shc_res3_a <- aggregate(shc_res3[,"cRand"], by = list(pts %/% 10 *10), FUN=mean)

par(mar=c(4, 3, 0.2, 0.2), mgp=c(1.5,0.7,0))
df1 <- data.frame('SHC(seq.)'=shc_res1[,"cRand"],'SHC(nm=3)'=shc_res2[,"cRand"],'SHC(nm=6)'=shc_res3[,"cRand"],check.names=F)
mns <- colMeans(df1)
df1 <- df1[,order(mns,decreasing=T)]
barplot(colMeans(df1), las=2, cex.names=.6, cex.axis=.6, cex.lab=.8, ylab="Avg. Corr. Rand", ylim=c(0,1))


###################################################
### code chunk number 72: SigmaIndex.Rnw:1475-1483
###################################################
sigma <- matrix(data=c(1,0,0,1),nrow=2,ncol=2)
mu <- list(P1=c(20,20), P2=c(27,20), P3=c(20,30), P4=c(5,20), P5=c(20,5))
si1 <- SigmaIndex(theta = 3.2, neighborhood = 6.4)
for(mu_n in names(mu)) addPopulation(si1,mu_n,mu[[mu_n]],sigma,200)
si2 <- SigmaIndex(theta = 3.2, neighborhood = 9.6)
for(mu_n in names(mu)) addPopulation(si2,mu_n,mu[[mu_n]],sigma,200)
si3 <- SigmaIndex(theta = 3.2, neighborhood = 12.9)
for(mu_n in names(mu)) addPopulation(si3,mu_n,mu[[mu_n]],sigma,200)


###################################################
### code chunk number 73: si-theorems2
###################################################
hist1 <- getHistogram(si1)
par(mar=c(4, 3, 0.5, 0.7), mgp=c(1.5,0.5,0))
barplot(unname(hist1[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 74: si-theorems3
###################################################
hist2 <- getHistogram(si2)
par(mar=c(4, 3, 0.5, 0.7), mgp=c(1.5,0.5,0))
barplot(unname(hist2[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


###################################################
### code chunk number 75: si-theorems4
###################################################
hist3 <- getHistogram(si3)
par(mar=c(4, 3, 0.5, 0.7), mgp=c(1.5,0.5,0))
barplot(unname(hist3[1,]), cex.names=.7, cex.axis=.7, cex.lab=.8, ylab="Density", xlab="Comp. cost reduction (%)",space=0, beside=T)
axis(side=1,at=c(0,20,40,60,80,100),cex.axis=.7)


