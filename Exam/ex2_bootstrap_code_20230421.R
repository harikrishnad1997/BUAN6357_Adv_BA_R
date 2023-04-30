###
# 
# BUAN 6357 Spring 2023 (Johnston) - Exam 2: Bootstrap
#
###
options(scipen=10,width=70)

wd  <- "c:/data/BUAN6357/exams/exam2" # change as needed
setwd(wd)
raw <- read.csv("olsData.csv")
raw <- raw[complete.cases(raw),]

require(tidyverse)
require(broom)
require(data.table)

s   <- 410350638

mdl <- V13~.

uF1 <- function (lbl,v,alpha=0.05) {
  z    <- qnorm(1-(alpha/2));  m  <- mean(v);   s  <- sd(v)
  lbP  <- m-z*s;     ubP <- m+z*s
  lo   <- alpha/2;   hi  <- 1-lo
  ci   <- quantile(v, c(lo,hi)); mu  <- mean((v-m)^2)
  tnm  <- quantile(v, c(.25, .5, .75))
  return (data.table(lbl=lbl, m=m, s=s, z=z, lbP=lbP, ubP=ubP, 
                     lbNP=ci[1], ubNP=ci[2],
                     q1= tnm[1], q2=tnm[2], q3=tnm[3],  
                     mse=mu, rmse=sqrt(mu)))
}

r0  <- lm(mdl,data=raw)
s0  <- tidy(r0,conf.int=T)
as.data.frame(s0)

(n  <- nrow(raw) )
(b  <- 750       )

set.seed(s)
t   <- data.table(grp=rep(1:b,each=n), idx=sample(n,n*b,replace=T))
b1  <- t[,tidy(lm(mdl,data=raw[idx,]),conf.int=T), by=grp]
(s1 <- b1[term!="(Intercept)",uF1("b1",estimate),by=term] )

set.seed(s)
t   <- data.table(grp=rep(1:b, each=n), idx=sample(rep(1:n,b)))
b2  <- t[,tidy(lm(mdl,data=raw[idx,]),conf.int=T), by=grp]
(s2 <- b2[term!="(Intercept)",uF1("b2",estimate),by=term] )

