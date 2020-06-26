library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library("minpack.lm")

# Support Functions ------------------------------------------------

# extract log scale ticks
GetTick <- function(Min, Max) {
  LowTick <- 10 ^ floor(log10(Min))
  LowTick <- floor(Min / LowTick) * LowTick
  
  HighTick <- 10 ^ floor(log10(Max))
  HighTick <- ceiling(Max / HighTick) * HighTick
  
  Pwr <- floor(log10(LowTick))
  Idx <- LowTick / (10 ^ Pwr)
  
  Ticks <- NULL
  TickLab <- NULL
  
  while (Idx * 10 ^ Pwr <= HighTick) {
    Ticks <- c(Ticks, Idx * 10 ^ Pwr)
    
    if (Idx == 1) {
      TickLab <- c(TickLab, bquote("" ^  ~ .(Pwr)))
    } else {
      TickLab <- c(TickLab, "")
    }
    
    if (Idx < 9) {
      Idx <- Idx + 1
    } else {
      Idx <- 1
      Pwr <- Pwr + 1
    }
    
  }
  
  AllInfo <- list()
  AllInfo$Ticks <- Ticks
  AllInfo$TickLab <- as.vector(TickLab)
  
  return(AllInfo)
  
}

#fCalculate AIC
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

# Models--------------------------------------------------------------

# CIM1********************************************
# 
#list of fixed parameters from literature
prs_CIM1 <- list(
  delta = 0.0015, h=50, n=1, tau=0.0017, a=1/(33.33 *365)#fixed paramaters
)

#proportion of the T cell population with x no. of divisions at time t
f_CIM1 <- function(x, t, beta, mu, parms=prs_CIM) with(parms, {
  d <- 1 - beta - delta
  #sum of f_0's
  sumj0 <- dpois(x, mu) + tau / d *dpois(x, mu)* sum((d*exp(a))^(0:-x))
  #  coeffs <- d^(-(1:x))
  sumj_all <-  sumj0 +  sum(sapply(1:x,function(j) 
    d^(-j)*prod(2*beta*(h^n)/(h^n+(x - (1:j))^n))*
      (choose(t,j)*dpois(x - j, mu)+tau/d*dpois(x - j, mu)*
         sum((d*exp(a))^-(0:(x - j))*choose(x - 0:(x - j), j))))
  )
  return(d^(t) * sumj_all)
})

#calculating the immunoscence predictor
getPred_CIM1 <- function(mdl,t, parms, fxdParms){
  diffrce <- t[[1]][1]
  prdctns <- sapply(1:length(t), function(i) {
    sapply(1:length(t[[i]]), function(j) {
      coll <- sapply(0:100, mdl, t = (t[[i]][j]-diffrce)*365, parms$beta,parms$mu, 
                     fxdParms)
      #normalising
      f_t <- coll/sum(coll)
      sum(f_t[50:100])})}, simplify = FALSE)
}

#calculating the disparity between predictions and observed data
residFun_CIM1 <- function(par, mdl,  observed, t, fxdParms){
  prdctns1 <- getPred_CIM1(mdl, t, par, fxdParms)
  prdctns <- sapply(prdctns1, mean)# mapply(function(x,y) sum(x * y)/sum(y), unWghtdPreds,wghts)
  resids <- log10(prdctns) - log10(observed)
  return(ifelse(is.nan(resids), 1e6, resids))
}


# CIM2 *************************************
 
#list of fixed parameters from literature
prs_CIM2 <- list(
  h=50, n=1, a=1/(33.33*365)#delta = 0.001, 
)

#proportion of the T cell population with x no. of divisions at time t
f_CIM2 <- function(x, t, tau, mu, beta, delta, fxdParms=prs_CIM2) with(fxdParms, {
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

#*****Calculate aging metric whilst updating B_t bi-annually *****
getPred_CIM2 <- function(mdl, t, parms, fxdParms){
  beta <- parms$delta
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

# Using best fit parameters to predict model output
bestPrdctns <- function(mdl, t, parms, fxdParms){
  beta <- parms$delta
  propnsList <- list()
  betaList <- list()
  for (i in 1:length(t)){
    mylist2 <- c()
    mylist3 <- c()
    for(j in 1:length(t[[i]])){
      mylist3[j] <- beta
      #calculating the propotions for all 0...100 cell divisions 
      prb <- sapply(0:100, mdl, t = t[[i]][j], parms$tau, parms$mu,
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

#calculating the disparity between predictions and observed data
residFun_CIM2 <- function(par, mdl,  observed, t, fxdParms){
  unWghtdPreds <-  getPred_CIM2(mdl, t, par, fxdParms)
  #prdctns <-  mapply(function(x,y) sum(x * y)/sum(y), unWghtdPreds,wghts)
  #resids <- log10(prdctns) - log10(observed)
  resids <- log(unlist(lapply(unWghtdPreds, mean))) - log(observed)
  return(ifelse(is.nan(resids), 1e6, resids))
}

# Infectios Diseases --------------------------------------------------

allDiseases=c(  
  "MRSA", "West Nile Virus Disease", "COVID-19", "group A Streptococcus",
  "Streptococcus pneumoniae",  "group B Streptococcus", "Legionellosis"
  ,"Haemophilus influenzae", "Influenza A", "Neisseria meningitidis"
)

viralDiseases <- c("COVID-19","West Nile Virus Disease","Influenza A")

bacterialDiseases <- c(  
  "MRSA",  "group A Streptococcus",
  "Streptococcus pneumoniae",  "group B Streptococcus", "Legionellosis"
  ,"Haemophilus influenzae", "Neisseria meningitidis"
)

# to read-in previously computed CIM predictions
# prevPreds_CIM1 <- readRDS("prevInfectiousCIM1Preds.rds")
# prevFits_CIM1 <- readRDS("prevInfectiousCIM1Fits.rds")
# 
# prevPreds_CIM2 <- readRDS("prevInfectiousCIM2Preds.rds")

myData1 <- NULL; myData2 <- NULL

prevPreds_CIM1 <- NULL
for (disease in viralDiseases){
  BgkCol <- "green" # BgkCol <- "yellow" for Bacterial diseases
  
  FittedData <- read.csv(paste0("infectiousDisease/", disease,".csv"))  
  
  # For disease data without confidence intervals
  if (disease %in% c("MRSA","Influenza A")){
    if (disease == "MRSA"){
      above19 <- FittedData$Age.lb>=18
    }else{
      above19 <- FittedData$Age.lb>=5
    }
    # filtering age to 18 or 20 yrs and above, when we iniate the model
    FittedData2 <- FittedData[above19,]
    # time points to run the model
    tpts <-  mapply(function(i,j) i:j, FittedData2$Age.lb, FittedData2$Age.ub, 
                    SIMPLIFY = FALSE)
    # observations to fit against
    obs <- FittedData2$`Rate.per.100.000`/ max(FittedData2$`Rate.per.100.000`)
    
    # Since model takes time, some values were previously run and saved
    if(disease %in% names(prevPreds_CIM1)){
      fit_CIM1 <- prevFits_CIM1[[disease]]
      prdctns_CIM1 <- prevPreds_CIM1[[disease]]
      
      fit_CIM2 <- try(readRDS(paste0("infectiousDisease/fits/model3fits/",disease,".rds")[1]))

      prdctns_CIM2 <- prevPreds_CIM2[[disease]]
    }
    else{ # For non saved model outputs, we reproduce by initially fitting
      
      # CIM1 fitting and predictions
      inits_CIM1 <- list(beta=0.00119252,mu=35.2693); a_CIM1 <- c(0,0); b_CIM1 <- c(1,100)
      
      fit_CIM1 <- nls.lm(par = inits_CIM1, lower = a_CIM1, upper = b_CIM1,
                    fn = residFun_CIM1, mdl = f_CIM1, fxdParms = prs_CIM1, observed = obs,
                    t = tpts,  control = nls.lm.control(nprint=1, ftol = 1e-3))
      
      prdctns1 <- getPred_CIM1(f_CIM1, tpts, fit_CIM1$par, prs_CIM1)
      prdctns_CIM1 <- unlist(lapply(prdctns1, mean))
      
      # CIM2 fitting and predictions
      lb <- FittedData2$Age.lb
      ub <- FittedData2$Age.ub
      
      lb2 <- (lb-lb[1])*365+365/2
      ub2 <- (ub-lb[1])*365
      
      tpts_CIM2 <-  mapply(function(i,j) seq(i,j,365/2), lb2, ub2,
                      SIMPLIFY = FALSE)
      
      inits_CIM2 <- list(mu=35.1138, tau=0.00255166, delta = 0.00124666)
      a_CIM2 <- c(0,0.00005,0.00005); b_CIM2 <- c(100,1, 1)
      
      fit_CIM2 <- nls.lm(par = inits_CIM2, lower = a_CIM2, upper = b_CIM2,
                    fn = residFun_CIM2, mdl = f_CIM2, fxdParms = prs_CIM2, observed = obs,
                    t = tpts_CIM2,  control = nls.lm.control(nprint=1, ftol = 1e-3))
      
      prdctns2 <-  bestPrdctns(f_CIM2, tpts_CIM2, fit_CIM2$par, prs_CIM2)
      prdctns_CIM2 <- unlist(lapply(prdctns2[[1]], mean))
      
    }
    # extract and store fit metrics; AIC and Rsquared
    AIC_CIM1 <- 2 * (length(fit_CIM1$par) + 1) - 2 * logL(fit_CIM1)
    Rsquare_CIM1 <- (cor(log(obs),log(prdctns_CIM1)))^2
    myData1 <- rbind(myData1,c(disease, unlist(fit_CIM1$par), fit_CIM1$deviance, AIC_CIM1, Rsquare_CIM1))
    
    AIC_CIM2 <- 2 * (length(fit_CIM2$par) + 1) - 2 * logL(fit_CIM2)
    Rsquare_CIM2 <- (cor(log(obs),log(prdctns_CIM2)))^2
    myData2 <- rbind(myData2,c(disease, unlist(fit_CIM2$par), fit_CIM2$deviance, AIC_CIM2, Rsquare_CIM2))
    
    FittedData$obs <- FittedData$Rate.per.100.000/max(FittedData$Rate.per.100.000)
    
    FittedData$CIM1 <- c(rep(NA, sum(!above19)),prdctns_CIM1)
    FittedData$CIM2 <- c(rep(NA, sum(!above19)),prdctns_CIM2)
    
    # get ticks for the log scale
    Tks <- GetTick(Min = max(min(c(0.01, FittedData$obs, na.rm = TRUE))),
                   Max = max(FittedData$obs, na.rm = TRUE))
    
    # plot
    p <- ggplot(
      data = FittedData,
      mapping = aes(
        x = Age.ub,
        y = obs
      )
    )
    
    p <- p +  theme(
      plot.title=element_text(face="bold", size=25, hjust = 0.5),
      axis.text=element_text(size=25),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = adjustcolor(BgkCol, alpha.f = 0.1)),
      plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
    ) + labs(x = "Age", y = bquote('log'[10]*' Normalised Incidence')) +
      ggtitle(trimws(disease)) +
      geom_point(color="blue", size = 3) + geom_line(color="blue", size = 2)
    
    p <- p + scale_y_log10(
      breaks = Tks$Ticks,
      labels = Tks$TickLab,
      limits = range(Tks$Ticks)
    ) + scale_x_continuous(
      breaks = c(0, 45, 90),
      ##limits = c(0, 90)
    )
    
    p <- p +  geom_line(
      aes(x = Age.ub, y = CIM2),
      color = adjustcolor("#FF3333", alpha.f = 1),
      size = 2)+  geom_line(
        aes(x = Age.ub, y = CIM1),
        color = adjustcolor("#FFCC00", alpha.f = 0.7),
        size = 2)
    print(p)
  }
  else{ # For disease data with confidence intervals
    n <- FittedData$Sample.Size
    z <- 1.96 #95% confidence level
    
    # calculate and store CI
    if (disease %in% c("COVID-19")){
      prp <- FittedData$No..of.Cases/n
      FittedData$Lower_Confidence_Interval <- (prp+z^2/(2*n)-z*sqrt((prp*(1-prp)/n)+z^2/(4*n^2)))/(1+z^2/n)*1e2/max(FittedData$Rate.per.100.000)
      FittedData$Upper_Confidence_Interval <- (prp+z^2/(2*n)+z*sqrt((prp*(1-prp)/n)+z^2/(4*n^2)))/(1+z^2/n)*1e2/max(FittedData$Rate.per.100.000)
      
      above19 <- FittedData$Age.lb>=20
    }
    else if (disease %in% c("West Nile Virus Disease","Legionellosis")){
      prp <- FittedData$Rate.per.100.000*1e-5
      FittedData$Lower_Confidence_Interval <- (prp+z^2/(2*n)-z*sqrt((prp*(1-prp)/n)+z^2/(4*n^2)))/(1+z^2/n)*1e5/max(FittedData$Rate.per.100.000)
      FittedData$Upper_Confidence_Interval <- (prp+z^2/(2*n)+z*sqrt((prp*(1-prp)/n)+z^2/(4*n^2)))/(1+z^2/n)*1e5/max(FittedData$Rate.per.100.000)
      
      above19 <- FittedData$Age.lb>=20
    }
    else{
      prp <- FittedData$Rate.per.100.000*1e-5
      FittedData$Lower_Confidence_Interval <- (prp+z^2/(2*n)-z*sqrt((prp*(1-prp)/n)+z^2/(4*n^2)))/(1+z^2/n)*1e5/max(FittedData$Rate.per.100.000)
      FittedData$Upper_Confidence_Interval <- (prp+z^2/(2*n)+z*sqrt((prp*(1-prp)/n)+z^2/(4*n^2)))/(1+z^2/n)*1e5/max(FittedData$Rate.per.100.000)
      
      above19 <- FittedData$Age.lb>=18
    }
    
    FittedData2 <- FittedData[above19,]
    
    tpts <-  mapply(function(i,j) i:j, FittedData2$Age.lb, FittedData2$Age.ub, 
                    SIMPLIFY = FALSE)
    
    obs <- FittedData2$`Rate.per.100.000`/ max(FittedData2$`Rate.per.100.000`)
    
    if(disease %in% names(prevPreds_CIM1)){
      fit_CIM1 <- prevFits_CIM1[[disease]]
      prdctns_CIM1 <- prevPreds_CIM1[[disease]]
      
      fit_CIM2 <- try(readRDS(paste0("infectiousDisease/fits/model3fits/",disease,".rds")[1]))
      
      prdctns_CIM2 <- prevPreds_CIM2[[disease]]
      
      
    }
    else{
      # CIM1 fitting and predictions
      inits_CIM1 <- list(beta=0.00119252,mu=34.2693); a_CIM1 <- c(0,0); b_CIM1 <- c(1,100)
      fit_CIM1 <- nls.lm(par = inits_CIM1, lower = a_CIM1, upper = b_CIM1,
                         fn = residFun_CIM1, mdl = f_CIM1, fxdParms = prs_CIM1, observed = obs,
                         t = tpts,  control = nls.lm.control(nprint=1, ftol = 1e-3))
      prdctns1 <- getPred_CIM1(f_CIM1, tpts, fit_CIM1$par, prs_CIM1)
      prdctns_CIM1 <- unlist(lapply(prdctns1, mean))
      
      # CIM2 fitting and predictions
      lb <- FittedData2$Age.lb
      ub <- FittedData2$Age.ub
      
      lb2 <- (lb-lb[1])*365+365/2
      ub2 <- (ub-lb[1])*365
      
      tpts_CIM2 <-  mapply(function(i,j) seq(i,j,365/2), lb2, ub2,
                           SIMPLIFY = FALSE)
      
      inits_CIM2 <- list(mu=35.1138, tau=0.00255166, delta = 0.00124666)
      a_CIM2 <- c(0,0.00005,0.00005); b_CIM2 <- c(100,1, 1)
      
      fit_CIM2 <- nls.lm(par = inits_CIM2, lower = a_CIM2, upper = b_CIM2,
                         fn = residFun_CIM2, mdl = f_CIM2, fxdParms = prs_CIM2, observed = obs,
                         t = tpts_CIM2,  control = nls.lm.control(nprint=1, ftol = 1e-3))
      
      prdctns2 <-  bestPrdctns(f_CIM2, tpts_CIM2, fit_CIM2$par, prs_CIM2)
      prdctns_CIM2 <- unlist(lapply(prdctns2[[1]], mean))
      
    }
    # Calculate fit metrics
    AIC_CIM1 <- 2 * (length(fit_CIM1$par) + 1) - 2 * logL(fit_CIM1)
    Rsquare_CIM1 <- (cor(log(obs),log(prdctns_CIM1)))^2
    myData1 <- rbind(myData1,c(disease,unlist(fit_CIM1$par), fit_CIM1$deviance, AIC_CIM1, Rsquare_CIM1))
    
    AIC_CIM2 <- 2 * (length(fit_CIM2$par) + 1) - 2 * logL(fit_CIM2)
    Rsquare_CIM2 <- (cor(log(obs),log(prdctns_CIM2)))^2
    myData2 <- rbind(myData2,c(disease, unlist(fit_CIM2$par), fit_CIM2$deviance, AIC_CIM2, Rsquare_CIM2))    
    
    FittedData$obs <- FittedData$Rate.per.100.000/max(FittedData$Rate.per.100.000)
    
    FittedData$CIM1 <- c(rep(NA, sum(!above19)),prdctns_CIM1)
    FittedData$CIM2 <- c(rep(NA, sum(!above19)),prdctns_CIM2)
    
    # Get ticks for log plot
    Tks <- GetTick(Min = max(c(0.001, min(FittedData$Lower_Confidence_Interval, na.rm = TRUE))),
                   Max = max(FittedData$Upper_Confidence_Interval, na.rm = TRUE))
    
    # Plot
    p <- ggplot(
      data = FittedData,
      mapping = aes(
        x = Age.ub,
        y = obs
      )
    )
    
    p <- p +  theme(
      plot.title=element_text(face="bold", size=25, hjust = 0.5),
      axis.text=element_text(size=25),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = adjustcolor(BgkCol, alpha.f = 0.1)),
      plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
    ) + labs(x = "Age", y = bquote('log'[10]*' Normalised Incidence')) +
      ggtitle(trimws(disease)) +
      geom_point(color="blue", size = 3) + geom_line(color="blue", size = 2)
    
    p <- p + scale_y_log10(
      breaks = Tks$Ticks,
      labels = Tks$TickLab,
      #limits = range(Tks$Ticks)
    ) + scale_x_continuous(
      breaks = c(0, 45, 90),
      #limits = c(0, 90)
    )
    
    p <- p  + geom_pointrange(color="blue",
                              aes(
                                ymin = Lower_Confidence_Interval,
                                ymax = Upper_Confidence_Interval
                              )
    )
    p <- p +  geom_line(
      aes(x = Age.ub, y = CIM2),
      color = adjustcolor("#FF3333", alpha.f = 1),
      size = 2)+  geom_line(
        aes(x = Age.ub, y = CIM1),
        color = adjustcolor("#FFCC00", alpha.f = 0.7),
        size = 2)
    print(p)
  }
}

colnames(myData1) <- c("Disease","beta","mu","RSS","AIC","Rsquare")
colnames(myData2) <- c("Disease","mu","tau", "delta","RSS","AIC","Rsquare")
myData1
myData2

# Pub plots
NewPlotList <- list(
  PlotList$`Chronic Myeloid Leukemia`,
  PlotList$`Soft Tissue including Heart`,
  PlotList$Brain,
  PlotList$`Colon and Rectum`,
  PlotList$Lip,
  PlotList$Stomach)

ml=plot_grid(plotlist=NewPlotList, ncol = 3)

ggsave(filename = "meso.pdf", plot = p, width = 8,  scale = 2)