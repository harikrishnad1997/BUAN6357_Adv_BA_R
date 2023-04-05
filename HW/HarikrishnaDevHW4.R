setwd("c:/data/BUAN6357/HW_4"); source("prep.txt", echo=T)

rm(list=ls())
demo <- T
require(MASS)
data("airquality")

require(data.table)
require(tidyverse)

# convert to data.table

raw     <- setDT(na.omit(airquality))[,-c("Day")]
raw$Month <- as.factor(raw$Month)

seed    <- 504737137

(n      <- nrow(raw))

fmla <- Ozone ~ .

base.resid <- lm(fmla,data=raw)$residuals

if(demo){summary(lm(fmla,data=raw))}

set.seed(seed)
frac <- 0.1
tst <- sample(n,ceiling(frac*n))
m1  <- lm(fmla,data=raw[-tst,])
if(demo){summary(m1)}
cv.trn <- m1$fitted.values
if(demo){summary(cv.trn)}
cv.tst <- predict(m1,raw[tst,])
if(demo){summary(cv.tst)}
cv.resid <- data.table(loc=tst,diff=raw$Ozone[tst]-cv.tst)
if(demo){summary(cv.resid)}
fitPred <- function(trn,tst,i) {
  m <- lm(fmla,data=trn)  
  y <- tst$Ozone    
  p <- predict(m,tst)     
  r <- y - p              
  return(data.table(idx=i,resid=r,y=y,yhat=p))
}

c2      <- data.table(grp=1:n, tst=1:n)
dt2     <- c2[,fitPred(raw[-tst,],raw[tst,],tst),by=.(grp)]
jk.resid      <- data.table(loc=dt2$idx,diff=dt2$resid)
if(demo){summary(jk.resid)}


set.seed(seed)
kf  <- 10
t       <- rep(1:kf,ceiling(n/kf))[1:n]
c3      <- data.table(k=t ,idx=sample(n))  
dt3     <- c3[,fitPred(raw[-idx,],raw[idx,],idx),by=.(k)]
kf.resid <- data.table(k=dt3$k,loc=dt3$idx,diff=dt3$resid)
if(demo){summary(kf.resid)}

ci <- function(res,a) {
alpha <- a 
# Parametric             
(z     <- qnorm(1-(alpha/2)) ) 
(t     <- mean(res) )             
(s     <- sd(res)   )            
(lbP   <- t-z*s   )           
(ubP   <- t+z*s   ) 
# Non parametric
(lo    <- alpha/2)  # lower bound definition
(hi    <- 1-lo)     # upper bound definition
(ci    <- quantile(res, c(lo,hi)) )
results <- list(par_ci=data.table("lb"=lbP,"ub"=ubP),non_par_ci=ci)
return(results)
}

ci(base.resid,0.05)
ci(cv.resid$diff,0.05)
ci(jk.resid$diff,0.05)
ci(kf.resid$diff,0.05)
source("validate.txt", echo=T)