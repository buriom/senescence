#*******************************************************************************
#In this implementation, I assumed a constant beta and  I fitted 2 parameters, 
#i.e beta and  mu, with tau fixed
#*******************************************************************************
library("minpack.lm")

#list of fixed parameters from literature
prs <- list(
  delta = 0.0015, h=50, n=1, tau=0.0017, a=1/(33.33 *365)#fixed paramaters
)

#fraction of T cells with x no. of divisions at time t
f <- function(x, t, beta, mu, parms=prs) with(parms, {
  d <- 1 - beta - delta
  #sum of f_0's
  sumj0 <- dpois(x, mu) + tau / d *dpois(x, mu)* sum((d*exp(a))^(0:-x))
  #  coeffs <- d^(-(1:x))
  sumj_all <-  sumj0 +  sum(sapply(1:x,
   function(j) d^(-j)*prod(2*beta*(h^n)/(h^n + (x - (1:j))^n))*(choose(t,j) *
        dpois(x - j, mu) +  tau/d*dpois(x - j, mu)*sum(
     (d*exp(a))^-(0:(x - j))*choose(x - 0:(x - j), j)
   )))
  )
  return(d^(t) * sumj_all)
})

#calculating the immunoscence predictor
getPred <- function(t, mdl, beta,mu, parms){
  #Predictions at time point t
  coll <- sapply(0:100, mdl, t = t*365, beta,mu, parms)
  #normalising
  f_t <- coll/sum(coll)
  #summing f_t >= h
  return(sum(f_t[50:100]))
}

#calculating the disparity between predictions and observed data
residFun <- function(par, mdl,  observed, t, wghts, parms){
  prdctns <- sapply(1:length(t), function(i) sum(wghts[[i]] * 
  sapply(1:length(t[[i]]), function(j) getPred((t[[i]][j] - 20), mdl, 
                                 par$beta, par$mu, parms))
  )/sum(wghts[[i]]))   
  resids <- log10(prdctns) - log10(observed)
  return(ifelse(is.nan(resids), 1e6, resids))
}

#********************************* data fitting ******************************
.args <- c("preProcessed/LipCancer_data.rds")

.args <- commandArgs(trailingOnly = TRUE)

#load the data and calculate the population weights in each class
incidenceData <- readRDS(.args[1])
#incidenceData <- readRDS("preProcessed/StomachCancer_data.rds")
incidenceData[dim(incidenceData)[1],2] <- 100
tpts <-  mapply(function(i,j) i:j, incidenceData$Age.lb, incidenceData$Age.ub)
# 4.18 and 0.04 obtained from https://www.census.gov/prod/cen2010/briefs/c2010br-03.pdf
wght <-  mapply(function(i,j) {things <- i:j; ifelse( 
  things<50, 1, 4.18-0.04*(things))},incidenceData$Age.lb, incidenceData$Age.ub)

#initial values for fitting
inits <- list(beta=0.00154582,mu=32.3559); a <- c(0,0); b <- c(1,100)
obs <- incidenceData$`Rate per 100,000`/ max(incidenceData$`Rate per 100,000`)
#fitting after normalising
fit <- nls.lm(par = inits, lower = a, upper = b,
  fn = residFun, mdl = f, parms = prs, observed = obs, t = tpts, wghts = wght,
  control = nls.lm.control(nprint=1,  ftol = 1e-4))
#predictions using the fitted values
prdctns <- sapply(1:length(tpts), 
            function(i) sum(wght[[i]] * sapply(1:length(tpts[[i]]),
    function(j) getPred((tpts[[i]][j] - 20), f, fit$par$beta,fit$par$mu, prs)))/
                    sum(wght[[i]])) 
#plotting
age <- incidenceData$Age.ub
logobs <- log10(obs)
preds <- log10(prdctns)
plot(age, logobs,  col = "blue",
     xlab = "Age", ylab = expression('log'[10]*' Incidence'), main = 
       expression(paste("Bones and joints Cancer Data with ",tau," fixed")))
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
text(70, (least + 0.43), cex = 1, bquote(paste('R'^2*' = ', .(Rsquare^2))))
text(70, (least + 0.13), cex = 1, bquote(paste('AIC = ', .(perfom))))
