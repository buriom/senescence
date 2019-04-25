#*******************************************************************************
#In this implementation, I assumed a constant beta and fitted 3 parameters, 
#i.e beta,  mu and tau
#*******************************************************************************
library("minpack.lm")

#list of fixed parameters from literature
prs <- list(
  delta = 0.001, h=50, n=1, a=1/(33.33 *365)
)

#fraction of T cells with x no. of divisions at time t
f <- function(x, t, beta, mu, tau, parms=prs) with(parms, {
  d <- 1 - beta - delta
  #sum of f_0's
  sumj0 <- dpois(x, mu) + tau / d *dpois(x, mu)* sum((d*exp(a))^(0:-x))
  #  coeffs <- d^(-(1:x))
  sumj_all <-  sumj0 +  sum(sapply(1:x,
     function(j) d^(-j)*prod(2*beta*(h^n)/(h^n + (x - (1:j))^n))*(choose(t,j)*
                dpois(x - j, mu) +  tau/d*dpois(x - j, mu)*sum(
       (d*exp(a))^-(0:(x - j))*choose(x - 0:(x - j), j)
     )))
  )
  return(d^(t) * sumj_all)
})

#calculating the immunoscence predictor
getPred <- function(mdl,t, parms, fxdParms){
  prdctns <- sapply(1:length(t), function(i) {
    sapply(1:length(t[[i]]), function(j) {
      coll <- sapply(0:100, mdl, t = (t[[i]][j]-20)*365, parms$beta,parms$mu,
                     parms$tau, fxdParms)
      #normalising
      f_t <- coll/sum(coll)
      sum(f_t[50:100])})}, simplify = FALSE)
}

#calculating the disparity between predictions and observed data
residFun <- function(par, mdl,  observed, t, wghts, fxdParms){
  unWghtdPreds <- getPred(mdl, t, par, fxdParms)
  prdctns <-  mapply(function(x,y) sum(x * y)/sum(y), unWghtdPreds,wghts)
  resids <- log10(prdctns) - log10(observed)
  return(ifelse(is.nan(resids), 1e6, resids))
}

#********************************* data fitting ******************************
.args <- c("preProcessed/StomachCancer_data.rds")

.args <- commandArgs(trailingOnly = TRUE)

#load the data and calculate the population weights in each class
incidenceData <- readRDS(.args[1])
#incidenceData <- readRDS("preProcessed/StomachCancer_data.rds")
#incidenceData[dim(incidenceData)[1],2] <- 100
tpts <-  mapply(function(i,j) i:j, incidenceData$Age.lb, incidenceData$Age.ub, 
                SIMPLIFY = FALSE)
#extract weight: 4.18 and 0.04 obtained from 
#https://www.census.gov/prod/cen2010/briefs/c2010br-03.pdf
wght <-  mapply(function(i,j) {things <- i:j; ifelse( 
  things<50, 1, 4.18-0.04*(things))},incidenceData$Age.lb, incidenceData$Age.ub, 
  SIMPLIFY = FALSE)

#initial values for fitting
inits <- list(beta=0.0005560091,mu=39.218,tau=246);a<-c(0,0,0);b<-c(1,100,10000)

obs <- incidenceData$`Rate per 100,000`/max(incidenceData$`Rate per 100,000`)
#fitting after normalising
fit <- nls.lm(par = inits, lower = a, upper = b,
              fn = residFun, mdl = f, fxdParms = prs,observed = obs,t = tpts,
              wghts = wght, control = nls.lm.control(nprint=1,  ftol = 1e-4))
#predictions using the fitted values
unWghtdPreds <- getPred(f, tpts, fit$par, prs)
prdctns <-  mapply(function(x,y) sum(x * y)/sum(y), unWghtdPreds,wght)
#plotting
age <- incidenceData$Age.ub
logObs <- log10(obs)
preds <- log10(prdctns)
cancer <- gsub(".*preProcessed/\\s*|_data.rds*", "", .args)
plot(age, logObs,  col = "blue",
     xlab = "Age", ylab = bquote('log'[10]*' Incidence'), 
     main = bquote(paste(.(cancer) ," data with ",tau," fitted")))
lines(age, preds)
#Calculating R^2
Rsquare <- cor(obs,preds); Rsquare^2

#Calculating the AIC
logL <- function(object, REML = FALSE, ...) { 
  
  res <- object$fvec 
  
  N <- length(res) 
  
  val <-  -N * (log(2 * pi) + 1 - log(N) + log(sum(res^2)))/2 
  
  ## the formula here corresponds to estimating sigma^2. 
  
  attr(val, "df") <- 1L + length(coef(object)) 
  
  attr(val, "nobs") <- attr(val, "nall") <- N 
  
  class(val) <- "logLik" 
  
  val 
  
}

perfom <- 2 * (length(fit$par) + 1) - 2 * logL(fit)

least <- min(obs)

text(70, (least + 0.05), cex = 1, bquote(paste('R'^2*' = ', .(Rsquare^2))))
text(70, (least + 0.13), cex = 1, bquote(paste('AIC = ', .(perfom))))
