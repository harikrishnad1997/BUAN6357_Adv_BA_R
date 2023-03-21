setwd("c:/data/BUAN6357/HW_3"); source("prep.txt", echo=T)

require(data.table)
require(partykit)       
require(tidyverse)
# for ctree()
rm(list = ls())
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

class_function <- function(seed,ite,classes) {
set.seed(seed)
minDigit  <- min(classes)
maxDigit  <- max(classes)
numDigits <- length(classes)
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
colnames(simDigits) <- c("digit", "s1", "s2", "s3", "s4", "s5", "s6", "s7")

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
#
td            <- simDigits  # create temp copy
fitted.logit  <- matrix(rep(NA,nrow(td)*numDigits), 
                        nrow=nrow(td) )
digits        <- td$digit   # save for later use (and re-use)
td$digit      <- NULL     

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

td$y        <- NULL

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


if (plots) plot(m)

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

td$fDigits    <- as.factor(digits)  # triggers classification

m             <- ctree(fDigits~.,
                       data=td)

fitted.tree1  <- predict(m)              # classification (caution: factor)
pprob.tree1   <- predict(m,type="prob")  # individual class probabilities 


t            <- mbr(pprob.tree1, classes, scale=F)
class.tree1  <- t$cat     # classifications
p.tree1      <- t$p.value
risk.tree1   <- 1-p.tree1 # Bayes Risk 


(hits.tree1   <- table(fitted.tree1, 
                       digits,
                       dnn=c("classif","actual") ) )
(pc.tree1     <- sum(diag(hits.tree1))/sum(hits.tree1)) # percent correct             


probabilty_function <- function(tbl) {
(d           <- diag(tbl) )

rSums        <- apply(tbl, byRows, sum) # classification outcomes
cSums        <- apply(tbl, byCols, sum) # actual underlying values

(rPC          <- d/rSums )                     # classification "pc" by digit
(cPC          <- d/cSums )                     # classification "pc" by actual underlying

(rBR          <- 1-rPC )                       # Bayes Risk by digit
(cBR          <- 1-cPC )   
results_BR <- list(r_prob = rPC,cProb = cPC,rBR=rBR,cBR = cBR)
return(results_BR)
}                    # Bayes Risk by actual underlying

prob.hist.logit <- probabilty_function(hits.logit)
prob.hist.tree.10 <- probabilty_function(hits.tree10)
prob.hist.tree.1 <- probabilty_function(hits.tree1)
dataframe <- data.frame(matrix(nrow = nrow(t1), ncol = 0))
dataframe$t1Pr <- p.tree1
dataframe$t1Cl <- class.tree1 
dataframe$t10Pr <- p.tree10
dataframe$t10Cl <- class.tree10
dataframe$lPr <- p.logit
dataframe$lCl <- class.logit
dataframe$digit <- c(t1)

final_results <- list(ite=ite,s=dataframe,logit.10.per_correct = pc.logit,tree.10.per_correct = pc.tree10,tree.1.per_correct = pc.tree1, 
                      logit.10.hits = hits.logit,tree.10.hits = hits.tree10,tree.1.hits = hits.tree1,
                      prob.hist.logit=prob.hist.logit,prob.hist.tree.10=prob.hist.tree.10,prob.hist.tree.1=prob.hist.tree.1)
return(final_results)
}

for(ite in n) {
  assign(paste0("f",ite),class_function(seed,ite,classes))
  assign(paste0("s",ite),class_function(seed,ite,classes)$s)
  print(paste("The loop has run for n = ",ite))
}
if(demo) {
for(ite in n) {
  print(get(paste0("f",ite)))
}
}
source("validate.txt", echo=T)

