library("SEER2R")

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
  #diffrce <- (20+(14-length(t))*5)
  beta <- parms$delta
  #prb <- dpois(0:100,parms$mu)
  mylist1 <- list()
  for (i in 1:length(t)){
    mylist2 <- c()
    for(j in 1:length(t[[i]])){
      #calculating the propotions for all 0...100 cell divisions 
      prb <- sapply(0:100, mdl, t = t[[i]][j], parms$tau, parms$mu,
                    beta, parms$delta, fxdParms)
      #updating beta at each time step, i.e after a year
      n <- fxdParms$n
      beta <- (parms$delta-parms$tau*exp(-fxdParms$a*t[[i]][j]))/
        (1 - 2 * sum(1/(1 + (fxdParms$h/(0:100))^n) * prb))
      #calculation immunoscence predictor by normalising the propotions from above
      mylist2[j] <- sum((prb/sum(prb))[50:100])
    }
    mylist1[[i]] <- mylist2 
  }
  return(mylist1)
}

#test
#inits <- list(mu=33.8672, tau=0.0011, delta = 0.00336) 
#system.time(b <- getPred(f, tpts,  inits, prs))
#unlist(lapply(b, mean))
#b

bestPrdctns <- function(mdl, t, parms, fxdParms){
  #diffrce <- (20+(14-length(t))*5)
  #prb <- dpois(0:100,parms$mu)
  beta <- parms$delta
  propnsList <- list()
  betaList <- list()
  for (i in 1:length(t)){
    mylist2 <- c()
    mylist3 <- c()
    for(j in 1:length(t[[i]])){
      mylist3[j] <- beta
      #calculating the propotions for all 0...100 cell divisions 
      prb <- sapply(0:100, f, t = t[[i]][j], parms$tau, parms$mu,
                    beta, parms$delta, fxdParms)
      n <- fxdParms$n
      beta <- (parms$delta-parms$tau*exp(-fxdParms$a*t[[i]][j]))/
        (1 - 2 * sum(1/(1 + (fxdParms$h/(0:100))^n) * prb))
      mylist2[j] <- sum((prb/sum(prb))[50:100])
    }
    propnsList[[i]] <- mylist2 
    betaList[[i]] <- mylist3
  }
  return(list(propnsList,betaList))
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
  
  if(Model == "PLIM"){
    Model.Pred <- NewX^Pars["c"]*(Pars["a"] / (exp(Pars["b"] * exp(-Pars["k"] * NewX)) - 1)) #Pars["a"]/(exp(exp(Pars["k"]*(NewX-Pars["t_0"])))-1)
  }
  
  return(Model.Pred)
  
}

AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic",
                         TXTfileName = "Incidence-Full-1Y-18.txt",
                         UseVarLabelsInData = TRUE)

allCancers <- unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))

monthlyB_t <-read.csv("resultsTable.csv")
allFits <- monthlyB_t$Cancer.Site


allCancers2 <- allCancers[allFits]

d <- (allCancers2[c(88,97)])
d[65] <- "Myeloid and Monocytic Leukemia"
thu <- split(d, ceiling(seq_along(d)/4))
mdl3Prdctns2 <- list()
pdf("newB_t3.pdf",title = "Immune models")
for (i in thu[22:25]) {
  par(mfrow = c(2,2))
    for(cancer in i){
    #load the data and calculate the population weights in each class
    FittedData <- readRDS(paste0("preProcessed/", gsub("\\s+","\\",cancer), "_data.rds")[1])
    jug <- FittedData$x
    lu <- jug -1
    FittedData$lb <- (lu-lu[1])*365+365/12
    FittedData$ub <- (jug-lu[1])*365
    tpts <-  mapply(function(i,j) seq(i,j,365/12), FittedData$lb, FittedData$ub,
    SIMPLIFY = FALSE)
    #fitting after normalising
    monthlyB_t <-read.csv("resultsTable.csv")
    fit <- monthlyB_t[which(monthlyB_t$Cancer.Site==cancer),]
    fit1 <- fit[3:5]
    names(fit1) <- c('delta',"mu","tau")
    #predictions using the fitted values
    prdctns1 <-  bestPrdctns(f, tpts, fit1, prs)
    prdctns <- unlist(lapply(prdctns1[[1]], mean))
    #mdl3Prdctns2[[cancer]] <- prdctns
    obs <- FittedData$y/ max(FittedData$y)
    age <- FittedData$x
    logObs <- log10(obs)
    preds <- log10(prdctns)
    plot(age, logObs,  col = "blue",
    xlab = "Age", ylab = bquote('log'[10]*' Normalised Incidence'), cex = 0.6,type = "p",pch =17,
    main = bquote(bold(paste(.(cancer)))),
    cex.lab=.8, cex.axis=.8, cex.main=.8, cex.sub=.8)
    lines(age, preds, col = "red", lwd = 2)
    
    fit2 <- unlist(fit[8:9])
    names(fit2) <- c("a","k")
    thing <- Predict.Model("A",fit2,age)
    lines(age,log10(thing), col = "green", lwd = 2)
    
    fit3 <- unlist(fit[14:16])
    names(fit3) <- c("a","b","k")
    thing <- Predict.Model("B",fit3,age)
    lines(age,log10(thing), col = "orange", lwd = 2)
    
    fit4 <- as.numeric(SumData[SumData[,1] == cancer,2:5])
    names(fit4) <- c("a","b","c","k")
    thing <- Predict.Model("PLIM",fit4,age)
    lines(age,log10(thing), col = "purple", lwd = 2)
    
    vertSpan <- range(logObs)
    least <- vertSpan[1]
    space <- (vertSpan[2]-vertSpan[1])/12
    legend(x=43,y=(least + 4*space),
    c(paste('IM-I: AIC=', round(unlist(fit[[10]]),1),', R^2=',round(unlist(fit[[11]]),2)),
    paste("IM-II: AIC=", round(unlist(fit[[17]]),1),', R^2=',round(unlist(fit[[18]]),2)),
    paste("IM-III: AIC=", round(unlist(fit[[6]]),1),', R^2=',round(unlist(fit[[7]]),2))),
    cex=.6,col=c("green", "orange", "red"),lwd=2)  }
}
dev.off()

