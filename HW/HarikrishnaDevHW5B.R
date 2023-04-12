setwd("c:/data/BUAN6357/HW_5"); source("prep.txt", echo=T)

rm(list=ls())
demo <- T
require(MASS)
data("airquality")

require(data.table)
require(tidyverse)
require(broom)

# convert to data.table

raw     <- setDT(na.omit(airquality))[,-c("Day")]
raw$Month <- as.factor(raw$Month)

seed    <- 337337168

n      <- nrow(raw) 
b <- 500 

fmla <- Ozone ~ .

set.seed(seed)
t2 <- sample(n, size=b*n, replace=T) 
c2 <- data.table(g=rep(1:b,each=n),idx=t2)
cl.dt <- c2[,tidy(lm(fmla,data=raw[idx,]),conf.int=T), by=g]
if(demo) {summary(cl.dt)
View(cl.dt)}


set.seed(seed)
t3 <- sample(rep(1:n,b))
c3 <- data.table(g=rep(1:b,each=n),idx=t3)
bal.dt <- c3[,tidy(lm(fmla,data=raw[idx,]),conf.int=T),by=g] 
if(demo) {summary(bal.dt)
View(bal.dt)}

ciPNP    <- function(v,ci) {  # data vector, 0 < CI < 1
         alpha <- 1-ci                 # area in the tails
         z     <- qnorm(1-alpha/2)     # Z for start of upper tail
         t     <- mean(v)              # for Paramteric
         s     <- sd(v)                # for Parametric
         lP    <- t-z*s                # lower bound for Parametric CI
         uP    <- t+z*s                # upper bound for Parametric CI
         bNP   <- c(alpha/2,1-alpha/2) # bounds for Non-Parametric CI
         ci    <- quantile(v, bNP)     # actual Non-Parametric CI bounds values
         return(data.table(avg=t,stDev=s,loP=lP,     hiP=uP,
                                         loNP=ci[1], hiNP=ci[2]) )
         }


coverage <- 0.95
ci2      <- cl.dt[,ciPNP(estimate,coverage),by=term]

ci3      <- bal.dt[,ciPNP(estimate,coverage),by=term]

m0 <- lm(fmla,data=raw)  
p0 <- m0$fitted.values  
r0 <- m0$residuals      
bs.dt <- as.data.table(tidy(m0,conf.int=T)) 
summary(m0)

source("validate.txt", echo=T)