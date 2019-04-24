#*******************************************************************************
#In this implementation, I consider a time varying beta, B_t and fit 2 parameters
#, i.e   mu and tau
#*******************************************************************************

library("minpack.lm")

#******************************  functions **********************************
#function for the fraction of T cells with x divisions at time t; f_t
#experimental parameters: delta = turn over rate (per day), C = initial number
# at t = 0 

#list of fixed parameters from literature
prs <- list(
  delta = 0.0017, h=50, n=1, a=1/(33.33*365)
)

#propotion of the T cell population with x no. of divisions at time t
f <- function(x, t, tau, mu, beta, fxdParms=prs) with(fxdParms, {
  d <- 1 - beta - delta
  #sum of f_0's
  sumj0 <- dpois(x, mu) + tau / d *dpois(x, mu)* sum((d*exp(a))^(0:x))
  #  coeffs <- d^(-(1:x))
  sumj_all <- sumj0+sum(sapply(1:x,function(j) d^(-j)*prod(2*beta*(h^n)/(h^n + 
                        (x - (1:j))^n))*(choose(t,j) * dpois(x - j, mu) +  
                        tau/d*dpois(x - j, mu)*sum((d*exp(a))^-(0:(x - j))*
                                                     choose(x - 0:(x - j), j)
                                   )))
  )
  return(d^(t) * sumj_all)
})

#*****Calculating aging whilst updating B_t at each time step (1 yr steps) *****
prdct <- function(mdl, t, parms, fxdParms){
  prb <- dpois(0:100,parms$mu)
  mylist1 <- list()
  for (i in 1:nrow(t)){
    mylist2 <- c()
    for(j in 1:ncol(t)){
      #updating beta at each time step, i.e after a year
      beta <- (fxdParms$delta - parms$tau * exp(-fxdParms$a*t[i,j]))/(1 - 2 * 
        sum(((0:100)^fxdParms$n/((0:100)^fxdParms$n + fxdParms$h^fxdParms$n)) * prb))
      #calculating the propotions for all 0...100 cell divisions 
      prb <- sapply(0:100, mdl, t = t[i,j], parms$tau, parms$mu, beta, fxdParms)
      #calculation immunoscence predictor by normalising the propotions from above
      mylist2[j] <- sum((prb/sum(prb))[50:100])
    }
    mylist1[[i]] <- mylist2 
  }
  return(mylist1)
}

#test
#system.time(b <- prdct(f, t(tpts),  inits, prs))
#b

bestPrdctns <- function(mdl, t, parms, fxdParms){
  prb <- dpois(0:100,parms$mu)
  propnsList <- list()
  betaList <- list()
  for (i in 1:nrow(t)){
    mylist2 <- c()
    mylist3 <- c()
    for(j in 1:ncol(t)){
      beta <- (fxdParms$delta - parms$tau * exp(-fxdParms$a*t[i,j]))/
        (1 - 2 * sum(((0:100)^fxdParms$n/((0:100)^fxdParms$n + 
        fxdParms$h^fxdParms$n)) * prb))
      mylist3[j] <- beta
      #calculating the propotions for all 0...100 cell divisions 
      prb <- sapply(0:100, mdl, t = t[i,j], parms$tau, parms$mu, beta, fxdParms)
      mylist2[j] <- sum((prb/sum(prb))[50:100])
    }
    propnsList[[i]] <- mylist2 
    betaList[[i]] <- mylist3
  }
  return(list(propnsList,betaList))
}


#calculating the disparity between predictions and observed data
residFun <- function(par, mdl,  observed, t, wghts, fxdParms){
  unWghtdPreds <-  prdct(mdl, t, par, fxdParms)
  unWghtdMat <- sapply(unWghtdPreds, function(x) as.numeric(unlist(x)))
  WgtdMat <- unWghtdMat*wghts
  prdctns <- sapply(1:ncol(WgtdMat),function(i) sum(WgtdMat[,i])/sum(wghts[,i]))
  resids <- log10(prdctns) - log10(observed)
  return(ifelse(is.nan(resids), 1e6, resids))
}

#Read in data and extract timepoints and weights
datafile <- "preProcessed/LipCancer_data.rds"
incidenceData <- readRDS(datafile)
incidenceData[nrow(incidenceData),2] <- 89#to make ub for last class = 89
hk <- head(seq(0,5,1),-1)
#extract timepoints
tpts <-  t(sapply(incidenceData$Age.lb, function(i) hk + (i-20))*365)
#extract weight
wght <-  mapply(function(i,j) {things <- i:j; ifelse( things<50, 1, 4.18-0.04*
                      (things))}, incidenceData$Age.lb, incidenceData$Age.ub)

#initial values for fitting
inits <- list(mu=44.48019, tau=0.003697389); a <- c(0,0); b <- c(100,1)

#fitting after normalising
fit <- nls.lm(par = inits, lower = a, upper = b,
              fn = residFun, mdl = f, fxdParms = prs, observed = 
              incidenceData$`Rate per 100,000`/max(incidenceData$`Rate per 100,000`),
              t = tpts, wghts = wght, control = nls.lm.control(nprint=1,
                                                               ftol = 1e-4))

#predictions using the fitted values
unWghtdPreds <-  bestPrdctns(f, tpts, fit$par, prs)[[1]]
unWghtdMat <- sapply(unWghtdPreds, function(x) as.numeric(unlist(x)))
WgtdMat <- unWghtdMat*wght
prdctns <- sapply(1:ncol(WgtdMat),function(i) sum(WgtdMat[,i])/sum(wght[,i]))
#plotting
age <- incidenceData$Age.ub
obs <- log10(incidenceData$`Rate per 100,000`/max(incidenceData$`Rate per 100,000`))
preds <- log10(prdctns)
cancer <- gsub(".*preProcessed/\\s*|_data.rds*", "", datafile)
plot(age, obs,  col = "blue",
     xlab = "Age", ylab = bquote('log'[10]*' Incidence'), main =
       bquote(paste(.(cancer) ," data with ",tau," fitted")))
lines(age, preds)
#Calculating R^2
Rsquare <- cor(obs,preds); Rsquare^2

#function for calculating the AIC
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


#calculate AIC
perfom <- 2 * (length(fit$par) + 1) - 2 * logL(fit)

#anotating the graph
least <- min(obs)
text(70, (least + 0.05), cex = 1, bquote(paste('R'^2*' = ', .(Rsquare^2))))
text(70, (least + 0.13), cex = 1, bquote(paste('AIC = ', .(perfom))))
