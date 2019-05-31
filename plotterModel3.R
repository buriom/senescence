#list of fixed parameters from literature
prs <- list(
  h=50, n=1, a=1/(33.33*365)#delta = 0.001, 
)

#proportion of the T cell population with x no. of divisions at time t
f <- function(x, t, tau, mu, beta, delta, fxdParms=prs) with(fxdParms, {
  d <- 1 - beta - delta
  #sum of f_0's
  sumj0 <- dpois(x, mu) + tau / d *dpois(x, mu)* sum((d*exp(a))^(0:x))
  #  coeffs <- d^(-(1:x))
  sumj_all <- sumj0+sum(sapply(1:x,function(j) 
    d^(-j)*prod(2*beta*(h^n)/(h^n+(x-(1:j))^n))*
      (choose(t,j)*dpois(x - j, mu) + tau/d*dpois(x - j, mu)*
         sum((d*exp(a))^-(0:(x - j))* choose(x - 0:(x - j), j))))
  )
  return(d^(t) * sumj_all)
})

#*****Calculating aging whilst updating B_t at each time step (1 yr steps) *****
getPred <- function(mdl, t, parms, fxdParms){
  t <- (t-min(t))*365
  beta <- parms$delta
  #prb <- dpois(0:100,parms$mu)
  mylist2 <- c()
  for(j in 1:length(t)){
    #calculating the propotions for all 0...100 cell divisions 
    prb <- sapply(0:100, mdl, t = t[j], parms$tau, parms$mu,
                  beta, parms$delta, fxdParms)
    #updating beta at each time step, i.e after a year
    beta <- (parms$delta-parms$tau*exp(-fxdParms$a*(t[j])))/
      (1 - 2 * sum(((0:100)^fxdParms$n/((0:100)^fxdParms$n + 
                                          fxdParms$h^fxdParms$n)) * prb))
    #calculation immunoscence predictor by normalising the propotions from above
    mylist2[j] <- sum((prb/sum(prb))[50:100])
  }
  return(mylist2)
}


Fit2ParModel <- function(FittedData, a=NULL, b=NULL, k_mid=0.044, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1), k_mid = 0.044)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out2_m <-
    nls(
      log(y) ~ log(a / (exp(
        b * exp(-k_mid * x)
      ) - 1)),
      data = FittedData,
      start = StartList,
      algorithm = "port",
      lower = LowVect,
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out2_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out2_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  if(Rounds > 0){
    for (i in 1:Rounds) {
      
      out2_m <-
        nls(
          log(y) ~ log(a / (exp(
            b * exp(-k_mid * x)
          ) - 1)),
          data = FittedData,
          start = StartList,
          lower = LowVect,
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      
      PEst <- summary(out2_m)$parameters[, "Estimate"]
      names(PEst) <- rownames(summary(out2_m)$parameters)
      StartList[names(PEst)] <- PEst
      
    }
  }
  
  R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA, NA)
  names(PVals) <- c("a", "b", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, b, k_mid)
  
  ReturnVect <- c(PVals, AIC(out2_m), R2Mb, deviance(out2_m), sum(summary(out2_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}

Fit1ParModel <- function(FittedData, a=NULL, k_mid=0.044, Rounds=10){
  
  StartList <- list(a = runif(1), k_mid = 0.044)
  
  StartList <- StartList[c(is.null(a), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out1_m <-
    nls(
      log(y) ~ log(a * exp(0.044 * x)),
      data = FittedData,
      start = StartList,
      lower = LowVect,
      algorithm = "port",
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out1_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out1_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  for (i in 1:Rounds) {
    
    
    out1_m <-
      nls(
        log(y) ~ log(a * exp(k_mid * x)),
        data = FittedData,
        start = StartList,
        lower = LowVect,
        algorithm = "port",
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    
    PEst <- summary(out1_m)$parameters[, "Estimate"]
    names(PEst) <- rownames(summary(out1_m)$parameters)
    StartList[names(PEst)] <- PEst
    
  }
  
  
  R1Mb <- 1 - sum((summary(out1_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA)
  names(PVals) <- c("a", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, k_mid)
  
  ReturnVect <- c(PVals, AIC(out1_m), R1Mb, deviance(out1_m), sum(summary(out1_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}



Predict.Model <- function(Model, Pars, NewX){
  
  if(Model == "A"){
    Model.Pred <- Pars["a"] * exp(Pars["k"] * NewX)
  }
  
  if(Model == "B"){
    Model.Pred <- Pars["a"] / (exp(Pars["b"] * exp(-Pars["k"] * NewX)) - 1)
  }
  
  if(Model == "C"){
    Model.Pred <- exp(Pars["b"])*NewX^Pars["a"]
  }
  
  if(Model == "TM"){
    Model.Pred <- Pars["a"]/(exp(exp(Pars["k"]*(NewX-Pars["t_0"])))-1)
  }
  
  return(Model.Pred)
  
}

#********************************* data fitting ******************************

#.args <- commandArgs(trailingOnly = TRUE)

d <- allCancers
thu <- split(d, ceiling(seq_along(d)/4))

pdf("B_t.pdf",title = "Changing Beta_t")
for (i in thu) {
  par(mfrow = c(2,2))
  for(cancer in i){
    #load the data and calculate the population weights in each class
    FittedData <- readRDS(paste0("preProcessed/", gsub("\\s+","\\",cancer), "_data.rds")[1])
    #incidenceData <- readRDS("preProcessed/StomachCancer_data.rds")
    mod2 <- Fit2ParModel(FittedData)
    mod1 <- Fit1ParModel(FittedData)
    
    obs <- FittedData$y/ max(FittedData$y)
    #fitting after normalising
    
    fit <- readRDS(paste0("fits/model3fits/",gsub("\\s+","\\",cancer),".rds")[1])
    
    prdctns <- getPred(f, FittedData$x, list(delta = fit[1], mu = fit[2], tau = fit[3]), prs)
    
    #jpeg(paste0('figures/model1Figs/', cancer,'.jpg'))
    
    age <- FittedData$x
    logObs <- log10(obs)
    preds <- log10(prdctns/max(prdctns))
    plot(age, logObs,  col = "blue",
         xlab = "Age", ylab = bquote('log'[10]*' Incidence'), cex = 0.8,type = "p",pch =17,
         main = bquote(bold(paste(.(cancer) ," with ",beta[t]))),
         cex.lab=.8, cex.axis=.8, cex.main=.8, cex.sub=.8)
    lines(age, preds, col = "red", lwd = 2)
    thing <- Predict.Model("A",mod1[1:2],age)
    lines(age,log10(thing/max(thing)), col = "green", lwd = 2)
    thing2 <- Predict.Model("B",c("a" = 375.216829, "b" = 16.2304424, "k" = 0.044),age)
    lines(age,log10(thing2/max(thing2)), col = "orange", lwd = 2)
    vertSpan <- range(logObs)
    least <- vertSpan[1]
    space <- (vertSpan[2]-vertSpan[1])/12
    
    legend(x=55,y=(least + 4*space),c("IM-I","IM-II","Our 2nd Model"),cex=.8,col=c("green", "orange", "red"),lwd=2)
  }
}
dev.off()
