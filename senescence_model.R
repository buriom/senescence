#Model for immunoscenesense via hayflick limit
#Author: Gershom Buri
#Date: 14th Jan 2019
#*******************************************************************************
library(parallel)
library("minpack.lm")


#******************************  functions **********************************

#Implementing the Prod symbol in the equation f_t
mltpy <- function(j, x, beta, h, n){
  total <- prod(sapply(1:j,function(k) 2 * beta * h ^n/ (h^n + (x-k)^n)))
  return(total)
}

#function for the fraction of T cells with x divisions at time t; f_t
#experimental parameters: delta = turn over rate (per day), C = initial number
# at t = 0 
f <- function(x, t, beta, delta=0.0015, h=50, n=6, C=0.0017, a=1/(33.33*12), mu=9){
   d <- 1 - beta - delta
  #sum of f_0's
  total1 <- dpois(x, mu)*(1 + C / d * sum(sapply(0:x, 
         function(i) exp(-a * i) *  d^(- i)))) 
  
  total2 <-  total1 +  sum(sapply(1:x, function(j) d^(-j) * dpois(x - j, mu) * 
          mltpy(j, x, beta, h, n) *  (choose(t,j) +  C / d * 
        sum(sapply(0:(x - j), function(i) choose(x - i, j) * 
                     exp(-a * i) * d^(-i))))))
  
    return(d^(t) * total2)
}

#Calculating aging 
getPred <- function(t, beta){
  #Predictions at time point t
  coll <- sapply(0:100, f, t = t*365, beta)
  #normalising
  f_t <- coll/sum(coll)
  #summing f_t >= h
  return(sum(f_t[50:100]))
}

#residual function for fitting; based on log-transformed values
residFun <- function(par, observed, yr){
  cl <- makeCluster(3)
  clusterExport(cl, c("f","mltpy", "getPred", "par", "t"))
  prdctns <- parSapply(cl, (yr - 20), FUN = getPred, par$beta)
  stopCluster(cl)
 #prdctns <- sapply((t - 20), getPred, par$beta)
  resids <- log10(par$M * prdctns) - log10(observed)#par$M * prdctns - observed#
  return(ifelse(is.nan(resids), 1e6, resids))
}

inits <- list(beta =  3.10975e-05, M = 1.88351e+20)
L <-  c(0, 0) 
U = c(0.001,1e23)