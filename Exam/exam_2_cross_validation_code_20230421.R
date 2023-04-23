###
# 
# BUAN 6357 Spring 2023 (Johnston) - Exam 2: Cross Validation
#
###
options(scipen=10, width=70)

require(data.table)

wd  <- "c:/data/BUAN6357/exams/exam2" # change as needed
setwd(wd)

raw <- read.csv("olsData.csv")
raw <- raw[complete.cases(raw),]

s   <- 171057756
n   <- nrow(raw)
mdl <- V13~.

r0  <- lm(mdl,data=raw)

uF1 <- function(df,i) {
  m  <- lm(mdl,data=df[-i,])
  t  <- df$V13[i]-predict(m,df[i,])
  return(list(loc=i, diff=t) )
}

t   <- data.table(grp=1:n, i=1:n)
r1  <- t[,uF1(raw,i), by=.(grp)]

set.seed(s)
k   <- 10 

uF2 <- function(df, grps, tst) {
  m  <- lm(mdl,data=df[-tst,])
  t  <- df$V13[tst]-predict(m,df[tst,])
  return(list(k=grps, loc=tst, diff=t) )
}

g   <- rep(1:k,ceiling(n/k))[1:n]
t   <- data.table(idx=g, k=g, i=sample(1:n))
r2  <- t[,uF2(raw, k, i), by=.(idx)]

uF3 <- function (lbl,v,alpha=0.05) {
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

t0  <- uF3("0",r0$residuals)
t1  <- uF3("1",r1$diff)
t2  <- uF3("2",r2$diff)

(resids.ci <- rbindlist(list(t0, t1, t2)))
