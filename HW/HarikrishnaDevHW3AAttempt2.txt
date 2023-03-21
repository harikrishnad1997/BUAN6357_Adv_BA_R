setwd("c:/data/BUAN6357/HW_3"); source("prep.txt", echo=T)

require(data.table)
require(partykit)       
require(tidyverse)
# for ctree()
# rm(list = ls())
demo <- F

# constants
cols      <- 7          # number of readout segments
byRows    <- 1          # used by apply()
byCols    <- 2          # used by apply()
seed      <- 975773311          # not a great choice

## parameters
n         <- c(25,50,100,250,500,1000,2500,5000);n    # number of replications per digit (small for fast test)
p         <- 0.9   # probability a segment works correctly

# flags for various uses
debug     <- T
plots     <- F
verbose   <- F
demo3     <- T

classes   <- c(0,1,2,3,4,5,6,7,8,9)
minDigit  <- min(classes)
maxDigit  <- max(classes)
numDigits <- length(classes)

# setClass(classifications.analysis,slots = list(sample_size = "integer",actual_digits = "data.table",logit_classification = ) )

for (ite in n) {
set.seed(seed)

# digits (copies of each digit)
t1 <- rep(classes, ite)
# segment patterns for each digit
t2 <- c(1,1,1,0,1,1,1,
        0,0,1,0,0,1,0,
        1,0,1,1,1,0,1,
        1,0,1,1,0,1,1,
        0,1,1,1,0,1,0,
        1,1,0,1,0,1,1,
        0,1,0,1,1,1,1,
        1,0,1,0,0,1,0,
        1,1,1,1,1,1,1,
        1,1,1,1,0,1,0)

t3 <- rep(t2, ite)                 # noiseless data, reapeated "n" times
t4 <- rbinom(length(t3),         # decide where "noise" happens
             1,                  # number of Bernoulli trials
             1-p)                # probability of failure event (noise)

# flip the bits (randomly) to add noise [ t4 == 1 designates failure event ]
t5 <- ifelse(t4 == 1, 
             1-t3, 
             t3)

# reshape vector into matrix
t5                  <- matrix(data=t5,
                              nrow=length(classes)*ite, 
                              ncol=cols, 
                              byrow=T)
dim(t1)             <- c(length(t1), 1)
t6                  <- cbind(t1, t5);
simDigits           <- as.data.frame(t6)
# describe(t4)
colnames(simDigits) <- c("digit", "s1", "s2", "s3", "s4", "s5", "s6", "s7")

# step 1.3: user-defined function to classify and return P(correct classification)

# MBR classification w/ optional normalization
#
# df:      indexable structure with relative frequencies
# classes: category labels, in order with df
# scale:   normalize for MBR calculation (boolean) [ default is TRUE ]
#
# BR is 1-pcat when scale==T

mbr <- function(df, classes, scale=T) {
#  byRows <- 1 [ global ]
   idx     <- apply(df, byRows, which.max) # find category number
   cat     <- classes[idx]                 # translate to category label
   pcat    <- apply(df, byRows, max)
   if (scale) {                            # can be turned off (e.g. ctree w/ multinomial)
      sc   <- apply(df, byRows, sum)       # row normalization value (margin sum)
      pcat <- pcat/sc                      # normalized row max value
      }
    return (data.table(cat=cat,p.value=pcat) )
   }

# step 2.0: classification models 

# 
#
# fit classification models (3 approaches for now)
#
# 10x logit:   each v. other (indicator)
# 10x tree:    each v. other (indicator)
#  1x tree:    multinomial (factor)
#
#

# step 2.1: 10 logits, 1 per outcome 

###
#
# setup: 10x logit (each v. other)
#
td            <- simDigits  # create temp copy
fitted.logit  <- matrix(rep(NA,nrow(td)*numDigits), 
                        nrow=nrow(td) )
digits        <- td$digit   # save for later use (and re-use)
td$digit      <- NULL       # strip this out for convenience in model spec

# fit: get predictions from individual models
for ( i in 1:length(classes) ) {
    d                 <- classes[i]
    td$y              <- 0         # initialize
    td$y[digits == d] <- 1         # indicator for -each- digit
    m                 <- glm(y ~ ., 
                             data=td, 
                             family=binomial())
    fitted.logit[,i]  <- m$fitted.values
    }
    
if (debug) summary(m)  # model for last category (digit: 9)

# classify
t           <- mbr(fitted.logit, classes)
summary(t)
class.logit <- t$cat     # classifications
p.logit     <- t$p.value
risk.logit  <- 1-p.logit # Bayes Risk 


(hits.logit <- table(class.logit, 
                     digits,
                     dnn=c("classif","actual") ) )
(pc.logit   <- sum(diag(hits.logit))/sum(hits.logit)) # percent correct             
# 0.754

# cleanup
td$y        <- NULL

# quick query
if (verbose) (td9 <- td[digits == 9, ])

# step 2.2: 10 trees, 1 per outcome 

###
#
# 10x tree: each v. other
#
###
#
# setup
#
fitted.tree10 <- matrix(rep(NA,nrow(td)*numDigits), 
                        nrow=nrow(td) )

# fit: get predictions from individual models
for ( i in 1:length(classes) ) {
    d                  <- classes[i]
    td$y               <- 0         # initialize
    td$y[digits == d]  <- 1         # indicator for -each- digit
    m                  <- ctree(y ~ ., 
                                data=td)
    fitted.tree10[,i]  <- predict(m)
    }

m # tree structure for last category (d==9)

# quick look
if (plots) plot(m)

# cleanup
td$y         <- NULL


# classify
t            <- mbr(fitted.tree10, classes)
class.tree10 <- t$cat      # classifications
p.tree10     <- t$p.value
risk.tree10  <- 1-p.tree10 # Bayes Risk 


(hits.tree10 <- table(class.tree10, 
                      digits,
                      dnn=c("classif","actual") ) )
(pc.tree10   <- sum(diag(hits.tree10))/sum(hits.tree10)) # percent correct             
# 0.726

# step 2.3: 1 tree, factor outcome, multinomial model

###
#
# 1 tree: multinomial (factor variable)
#
###
#
# setup
#
td$fDigits    <- as.factor(digits)  # triggers classification

m             <- ctree(fDigits~.,
                       data=td)

fitted.tree1  <- predict(m)              # classification (caution: factor)
pprob.tree1   <- predict(m,type="prob")  # individual class probabilities 
# pprob.tree1 row stochastic? try: summary(apply(pprob.tree1,byRows,sum))
# find min Bayes Risk classification based on probabilities [ modified 20190926 ]
t            <- mbr(pprob.tree1, classes, scale=F)
class.tree1  <- t$cat     # classifications
p.tree1      <- t$p.value
risk.tree1   <- 1-p.tree1 # Bayes Risk 


(hits.tree1   <- table(fitted.tree1, 
                       digits,
                       dnn=c("classif","actual") ) )
(pc.tree1     <- sum(diag(hits.tree1))/sum(hits.tree1)) # percent correct             
# 0.716

# step 3.1: conditional Bayes Risk calculation

# comment: 1-pc.tree1 is overall BR (pc == proportion correct, i.e. empirical probability estimate
#          conditioning sets a pre-condition, do we know the "actual" (underlying) true classification or
#            do we know the classification assigned?  Both are valid questions to be answered (and easily 
#            answered at this point).  This will only work with data (training or testing) where
#            actual AND classification information is known, hence these are empirical probability
#            estimates which require representative data samples.

# working data: hist.tree1 (confusion matrix for multinomial tree)
if(demo)  {hits.tree1}

# correct classifications
(d           <- diag(hits.tree1) )

rSums        <- apply(hits.tree1, byRows, sum) # classification outcomes
cSums        <- apply(hits.tree1, byCols, sum) # actual underlying values

(rPC          <- d/rSums )                     # classification "pc" by digit
(cPC          <- d/cSums )                     # classification "pc" by actual underlying

(rBR          <- 1-rPC )                       # Bayes Risk by digit
(cBR          <- 1-cPC )                       # Bayes Risk by actual underlying

# classifications.analysis <- list(sample_size = ite,actual_digits = simDigits,)

if(ite == 25) {
    s25 <- data.frame(matrix(nrow = nrow(t1), ncol = 0))
    s25$t1Pr <- p.tree1
    s25$t1Cl <- class.tree1 
    s25$t10Pr <- p.tree10
    s25$t10Cl <- class.tree10
    s25$lPr <- p.logit
    s25$lCl <- class.logit
    s25$digit <- c(t1)
}

if(ite == 50) {
    s50 <- data.frame(matrix(nrow = nrow(t1), ncol = 0))
    s50$t1Pr <- p.tree1
    s50$t1Cl <- class.tree1 
    s50$t10Pr <- p.tree10
    s50$t10Cl <- class.tree10
    s50$lPr <- p.logit
    s50$lCl <- class.logit
    s50$digit <- c(t1)
}
if(demo) {
print("Loop kinda worked for ite:")
print(ite)
}
}
if(demo) {
summary(s25)
summary(s50)
}
source("validate.txt", echo=T)