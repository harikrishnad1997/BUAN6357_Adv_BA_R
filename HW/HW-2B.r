##Libraries
if(demo) {rm(list = ls())}
if(require(tidyverse)) {install.packages("tidyverse")}
if(require(broom)) {install.packages("broom")}
if(require(data.table)) {install.packages("data.table")}
library(tidyverse)  # or require(), must be installed
library(broom)
library(data.table)

# parameterization
demo               <- F   # expand outputs when using for demonstration/teaching

maxClusters        <- 10   # maximum number of clusters to be considered

byRows             <- 1    # for use with apply()
byCols             <- 2    # for use with apply()

seed               <- 893571057    # for demo purposes ONLY!

alpha              <- 0.05 # classic default

# initialize RNG
set.seed(seed)

# training set
kd       <- fread("HW_2_train.csv")
train <- kd
## Adding grps to help with comparing to Walter's code
kd$grps <- 0

if(demo) { str(kd) 
summary(kd)}
# testing set
kd2      <- fread("HW_2_test.csv")
test <- kd2
kd2$grps <- 0
if(demo) { str(kd2) 
summary(kd2)}

#--------------------------------------------------------------------------------------------------------------------------------#

## Kmeans 

# pre-allocate space
twssKM   <- rep(-1, maxClusters)
if(demo) { print(twssKM) }

# cluster the training set under varying scenarios, keep only TWGSS
tmp      <- kd[,c("v1","v2","v3","v4","v5")]        # working copy, w/o GRP info
for (i in 1:maxClusters) {
  set.seed(seed)
  twssKM[i] <- kmeans(tmp,                   # variables
                      i,                     # number of clusters
                      iter.max=15,           # max iterations
                      nstart=5)$tot.withinss # warm-up and TWGSS
}
if(demo) { print(twssKM) }

d1Calc        <- function(v) {               # v: vector
  n       <- length(v)
  d1      <- v[1:(n-1)] - v[2:n]         # first order differences
  d1scale <- d1/max(d1)                  # relative scale
  return (list(d1      = d1,
               d1scale = d1scale) )
}

(d1KM       <- d1Calc(twssKM) )              # d1scale has big changes at 3 and 5 
twssKM[4]

nClust     <- 4

# re-fit and validate
set.seed(seed)
tmp        <- kd[,c("v1","v2","v3","v4","v5")]    # working copy of "kd" w/o GRP info
kmClust    <- kmeans(tmp,
                     nClust,
                     iter.max=15,
                     nstart=5)$cluster

pcAgree   <- function(tbl) {
  return ( sum(apply(tbl, byRows, max))/sum(tbl) )
}                                     # compare values numerically

trn       <- kd            # working copy
trn$clust <- kmClust

str(trn)

tmp       <- kd[,c("v1","v2","v3","v4","v5")]  
summary(tmp)  # working copy w/o GRP info
d         <- dist(tmp)
hc        <- hclust(d^2) 
hcObj <- hc                   # sq Euclidean == WSS; default method="complete" (merge tree)
# plot(hc)                                    # dendrogram (rendering of merge tree in hc)

# objective 3.2: evaluate hclust scenarios

# use data.table() form for penalty fn calc
twssHC    <- rep(-1,maxClusters)                    # pre-allocate space
for (i in 1:maxClusters) {
  t0        <- cutree(hc,i)                       # define memberships for -i- clusters
  t1        <- data.table(c=t0,i=1:length(t0))    # transform to DT for processing
  t2        <- t1[,
                  .(twss=sum(scale(tmp[i], center=T, scale=F)^2) ),
                  by=c]                           # WGSS vector for clusters
  twssHC[i] <- sum(t2$twss)                       # TWGSS across clusters
}

if(demo) { print(twssHC) }

hcPen     <- twssHC
twssHC[4]

(d1HC     <- d1Calc(twssHC) )

# plot(1:length(d1HC$d1), 
#      d1HC$d1,
#      type="b",
#      main="Total Within-Group SS by Number of Clusters",
#      xlab="Number of Clusters",
#      ylab="Total With-Group SS (Hclust, D1)")

# plot(1:length(d1HC$d1scale), 
#      d1HC$d1scale,
#      type="b",
#      main="Total Within-Group SS by Number of Clusters",
#      xlab="Number of Clusters",
#      ylab="Total With-Group SS (Hclust, D1scale)")

# plot(4:length(d1HC$d1scale), 
#      d1HC$d1scale[4:length(d1HC$d1scale)],
#      type="b",
#      main="Total Within-Group SS by Number of Clusters",
#      xlab="Number of Clusters",
#      ylab="Total With-Group SS (Hclust, D1scale)")

# compare to grp and newClust w/ table
nClust <- 4
(hcClust   <- cutree(hc,nClust))

(tblHcKm  <- table(kmClust,hcClust, dnn=c("kmClust","hcClust" ) ) )

pcAgree(tblHcKm)

# objective 4.1: MHD calculation prep for testing data

# 
# prep for mahalanobis distance
# uses scale() to get means and std.dev; no loops
variability <- function(df) {
  df$clust   <- NULL
  df$grps    <- NULL
  n          <- nrow(df)
  df2        <- scale(df, center=T, scale=T)
  vcvinv     <- solve(cov(df2))                      # inverse of variance-covariance matrix
  return( list(n      = n,                           # create list to store heterogeneous data
               avg    = attr(df2, "scaled:center"),  # from meta-data: averages
               sdev   = attr(df2, "scaled:scale"),   # from meta-data: standard deviations
               vcvinv = vcvinv )
  )
}

trn$clust <- kmClust
mhWork   <- trn                       %>%
  group_by(clust)           %>%
  do( desc=variability(select(., v1, v2, v3, v4, v5) ) )

if(demo) { 
  print(mhWork)
  str(mhWork) }

(clusters <- mhWork$clust)
(desc     <- mhWork$desc)

(mhDf     <- length(desc[[1]]$avg))             # degrees of freedom for chi-sq
# observe list access use [[ ]]

# objective 4.3: reference set from training data

###
# calculate mahalanobis distances for training set to check behavior
###
# storage space for distances
d.train  <- matrix( -1,
                    nrow=nrow(trn),
                    ncol=length(clusters) )
str(d.train)
# clear out "grps", but hold onto it for later use

kdGrps   <- kd$grps
kd$grps  <- NULL
str(kd)
# objective 4.4: MHD to each cluster by row (observation)

# collect Mahalanobis distance for training data across clusters
for ( i in seq_along(clusters) ) {
  t            <- desc[[i]]
  tdf          <- scale(kd, 
                        center=t$avg, 
                        scale=t$sdev)   # center and scale w/ orig. cluster params
  d.train[, i] <- mahalanobis(tdf, 
                              center=F,       # data already centered
                              cov=t$vcvinv,   # covariance matrix pre-computed
                              inverted=T)     #   and already inverted
}

# objective 4.5: find closest cluster w/ probability 

# which cluster for each training observation?
mhdClust    <- apply(d.train, 
                     byRows, 
                     which.min)    # where is the min located? (col num is cluster)
str(mhdClust)
(kd$grps     <- trn$clust)
kd$cluster  <- mhdClust
(vtblkm        <- table(kmClust, 
                     kd$cluster, 
                     dnn=c("grps", "cluster") ))
pcAgree(vtblkm)

# chi-squared cdf value for each training observation
minStat  <- apply(d.train, 
                  byRows, 
                  min)             # what is the min value? (closest cluster)
chiSq    <- pchisq(minStat, 
                   df=mhDf, 
                   lower.tail=F)
summary(chiSq)

sum(chiSq<=alpha)/length(chiSq)    # proportion <= alpha (diagnostic) should be ~alpha

kd$minStat   <- minStat
(kd$chiSq     <- chiSq)
###
# storage space for distances
d2       <- matrix( -1,
                    nrow=nrow(kd2),
                    ncol=length(clusters) )

# hold on to "grps" info for accuracy check but remove from kd2 structure

kd2Grps  <- kd2$grps
kd2$grps <- NULL

# collect Mahalanobis distance for test data across clusters
for ( i in seq_along(clusters) ) {
  t       <- desc[[i]]
  tdf     <- scale(kd2, 
                   center=t$avg, 
                   scale=t$sdev)   # scale to Z w/ orig. cluster dist
  d2[, i] <- mahalanobis(tdf, 
                         center=F, 
                         cov=t$vcvinv, 
                         inverted=T)
}

# which cluster for each testing observation?
newClust    <- apply(d2, 
                     byRows, 
                     which.min) # location of min value by row
kd2$cluster <- newClust

set.seed(seed)
tmptst        <- kd2[,c("v1","v2","v3","v4","v5")]    # working copy of "kd" w/o GRP info
# d         <- dist(kd2)
# hc        <- hclust(d^2)   
# (kmtst   <- cutree(hc,nClust))
kmtst    <- kmeans(tmptst,
                     nClust,
                     iter.max=15,
                     nstart=5)$cluster
kd2$grps <- kmtst
table(kd2$grps, 
      kd2$cluster, 
      dnn=c("tstGrps", "tstCluster") )

# chi-squared cdf value for each testing observation
minStat  <- apply(d2, 
                  byRows, 
                  min)     # min value by row
chiSq    <- pchisq(minStat,
                   df=mhDf, 
                   lower.tail=F)
summary(chiSq)
sum(chiSq<=alpha)/length(chiSq) # diagnostic; should be ~alpha

kd2$minStat   <- minStat
kd2$chiSq     <- chiSq


#------------------------------------------------------------------------------------------------------------------------------
## hclust 

trnh <- trn
trnh$hclust <- hcClust
mhWorkhc   <- trnh                       %>%
  group_by(hclust)           %>%
  do( desc=variability(select(., v1, v2, v3, v4,v5) ) )

if(demo) { 
  print(mhWorkhc)
  str(mhWorkhc) }

(clusters <- mhWorkhc$hclust)
(desc     <- mhWorkhc$desc)

(mhDf     <- length(desc[[1]]$avg))             # degrees of freedom for chi-sq
# observe list access use [[ ]]

d.train  <- matrix( -1,
                    nrow=nrow(trnh),
                    ncol=length(clusters) )

kdh       <- fread("HW_2_train.csv")
kdh$grps <- 0

if(demo) { str(kdh) 
summary(kdh)}
# testing set
kdh2      <- fread("HW_2_test.csv")
kdh2$grps <- 0
if(demo) { str(kdh2) 
summary(kdh2)}

kdGrps   <- kdh$grps
kdh$grps  <- NULL
# objective 4.4: MHD to each cluster by row (observation)

# collect Mahalanobis distance for training data across clusters
for ( i in seq_along(clusters) ) {
  t            <- desc[[i]]
  tdf          <- scale(kdh, 
                        center=t$avg, 
                        scale=t$sdev)   # center and scale w/ orig. cluster params
  d.train[, i] <- mahalanobis(tdf, 
                              center=F,       # data already centered
                              cov=t$vcvinv,   # covariance matrix pre-computed
                              inverted=T)     #   and already inverted
}

# objective 4.5: find closest cluster w/ probability 

# which cluster for each training observation?
(mhdClust    <- apply(d.train, 
                     byRows, 
                     which.min) )   # where is the min located? (col num is cluster)
kdh$cluster  <- mhdClust
kdh$grps     <- trnh$hclust

(vtblHC        <- table(hcClust, 
                     kdh$cluster, 
                     dnn=c("hClust", "cluster") ))
pcAgree(vtblHC)

# chi-squared cdf value for each training observation
minStat  <- apply(d.train, 
                  byRows, 
                  min)             # what is the min value? (closest cluster)
chiSq    <- pchisq(minStat, 
                   df=mhDf, 
                   lower.tail=F)
summary(chiSq)

sum(chiSq<=alpha)/length(chiSq)    # proportion <= alpha (diagnostic) should be ~alpha

kd$minStat   <- minStat
kd$chiSq     <- chiSq

# calculate mahalanobis distances for testing set to compare behavior
###
# storage space for distances
d2       <- matrix( -1,
                    nrow=nrow(kdh2),
                    ncol=length(clusters) )

# hold on to "grps" info for accuracy check but remove from kd2 structure

kdh2Grps  <- kdh2$grps
kdh2$grps <- NULL

# collect Mahalanobis distance for test data across clusters
for ( i in seq_along(clusters) ) {
  t       <- desc[[i]]
  tdf     <- scale(kdh2, 
                   center=t$avg, 
                   scale=t$sdev)   # scale to Z w/ orig. cluster dist
  d2[, i] <- mahalanobis(tdf, 
                         center=F, 
                         cov=t$vcvinv, 
                         inverted=T)
}

# which cluster for each testing observation?
newClust    <- apply(d2, 
                     byRows, 
                     which.min) # location of min value by row
kdh2$cluster <- newClust

tmptst        <- kdh2[,c("v1","v2","v3","v4","v5")]    # working copy of "kd" w/o GRP info

d         <- dist(tmptst)
hc        <- hclust(d^2) 
hctst <- cutree(hc,nClust)
kdh2$grps <- hctst
table(hctst, 
      kdh2$cluster )

# chi-squared cdf value for each testing observation
minStat  <- apply(d2, 
                  byRows, 
                  min)     # min value by row
chiSq    <- pchisq(minStat,
                   df=mhDf, 
                   lower.tail=F)
summary(chiSq)

sum(chiSq<=alpha)/length(chiSq) # diagnostic; should be ~alpha



