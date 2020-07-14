library(dplyr)
library(ggplot2)
#library(xlsx)
library(SEER2R)
library("minpack.lm")


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

#calculating the disparity between predictions and observed data
residFun <- function(par, mdl,  observed, t, wghts, fxdParms){
  prdctns <-  getPred(mdl, t, par, fxdParms)
  #prdctns <-  mapply(function(x,y) sum(x * y)/sum(y), unWghtdPreds,wghts)
  resids <- log(unlist(lapply(prdctns, mean))) - log(observed)
  return(ifelse(is.nan(resids), 1e6, resids))
}

#function to Calculate the AIC
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

#********************************* data fitting ******************************



# Define support functions

Convert1Year <- function(Vect) {
  as.integer(trimws(strsplit(
    x = Vect, split = "years", fixed = TRUE
  )))
}

AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic",
                         TXTfileName = "Incidence-Full-1Y-18.txt",
                         UseVarLabelsInData = TRUE)
AllData$AgeAdjusted_Rate <- as.numeric(levels(AllData$AgeAdjusted_Rate))[AllData$AgeAdjusted_Rate]
AllData$Lower_Confidence_Interval <- as.numeric(levels(AllData$Lower_Confidence_Interval))[AllData$Lower_Confidence_Interval]
AllData$Upper_Confidence_Interval <- as.numeric(levels(AllData$Upper_Confidence_Interval))[AllData$Upper_Confidence_Interval]


#SumData <- NULL

# for(Cancer in unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))){
#   
#   if(Cancer == "All Sites"){
#     next()
#   }
#   
#   print(Cancer)
genderList <- c("Male", "Female")
Cancer <- "Salivary Gland"  
CancerList <- unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))[-(52:66)]
monthlyB_t <-read.csv("resultsTable.csv")
allFits <- monthlyB_t$Cancer.Site


CancerList <- list("Rectum and Rectosigmoid Junction","Rectosigmoid Junction") #as.character(resToFit) # CancerList[CancerList %in% allFits]
for (gender in genderList){
  for (Cancer in CancerList){
    if (Cancer %in% c("Skin excluding Basal and Squamous","Lip", "Esophagus")){
      next()
    }
    DataSubset <- AllData[trimws(AllData$`Site_recode_ICDO3/WHO_2008`) == Cancer &
                            AllData$Race_recode_W_B_AI_API == "White" &
                            trimws(AllData$Sex) == gender,]
    
    if(sum(DataSubset$AgeAdjusted_Rate, na.rm = TRUE)==0){
      next
    }
    
    Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
    Converted.Age[Converted.Age==84] <- NA
    
    DataSubset <- cbind(DataSubset, Converted.Age)
    
    XFil <- DataSubset$Converted.Age
    YFil <- DataSubset$AgeAdjusted_Rate
    
    ToFil <- YFil == 0 | XFil < 20 | is.na(XFil) | is.na(YFil)
    
    XFil <- XFil[!ToFil]
    YFil <- YFil[!ToFil]
    
    FittedData <- data.frame(cbind(XFil, YFil))
    colnames(FittedData) <- c("x", "y")
    # rescaled this way because min-max would result in log(o) = Inf
    #FittedData$y <-  FittedData$y/ max(FittedData$y)
    ting <- try(readRDS(paste0("fits/model3fits/",gsub("\\s+","\\",Cancer),".rds")[1]))
  
    # ting <- try(readRDS(paste0("Betas/Male/fits/",gsub("\\s+","\\",Cancer),".rds")[1])) 
    # ting <- unname(readRDS(paste0("Betas/",gender,"/fits/",Cancer,".rds")))
    inits <- list(mu=ting[2], tau=ting[3], delta = ting[1]); a <- c(0,1e-7,1e-7); b <- c(100,1, 1)
    
    jug <- FittedData$x
    lu <- jug -1
    FittedData$lb <- (lu-lu[1])*365+365/12
    FittedData$ub <- (jug-lu[1])*365
    
    tpts <-  mapply(function(i,j) seq(i,j,365/12), FittedData$lb, FittedData$ub, 
                    SIMPLIFY = FALSE)
    
    a <- c(0,1e-7,0.0005); b <- c(100,1, 1)
    obs <- FittedData$y/ max(FittedData$y)
    
    system.time(fit <- nls.lm(par = inits, lower = a, upper = b,
                              fn = residFun, mdl = f, fxdParms = prs, observed = obs, t = tpts,
                              control = nls.lm.control(nprint=1,  ftol = 1e-2)))
    
    #predictions using the fitted values
    prdctns1 <-  bestPrdctns(f, tpts, fit$par, prs)
    # betating  <- unlist(prdctns1[[2]])
    # plot((unlist(tpts)/365)[-(1:240)], betating[-(1:240)], type = "l")
    
    # add tpts to the stored data
    prdctns1[[3]] <- tpts
    prdctns1[[4]] <- FittedData$x
    prdctns1[[5]] <- FittedData$y
    saveRDS(prdctns1,paste0("Betas/",gender,"/predictions/",Cancer,".rds"))
    
    prdctns <- unlist(lapply(prdctns1[[1]], mean))
    
    Rsquare <- (cor(log(obs),log(prdctns)))^2
    
    #Calculating AIC
    perfom <- 2 * (length(fit$par) + 1) - 2 * logL(fit)
    
    fitResults <- c(delta = fit$par$delta, mu = fit$par$mu, tau = fit$par$tau, AIC = perfom, Rsqr = Rsquare)
    
    saveRDS(inits, paste0("Betas/",gender,"/fits/",Cancer,".rds"))  
  }
}

  # # incidencies: y - female, z - male
# plot(FittedData$x,FittedData$z, col = 'red', type = "l")
# lines(FittedData$x,FittedData$y, col = 'blue', type = "l")
# 
# # normalised incidences
# plot(FittedData$x,FittedData$z/max(FittedData$z), col = 'red', type = "l")
# lines(FittedData$x,FittedData$y/max(FittedData$y), col = 'blue', type = "l")
# 
# # normalised and then log-transformed
# plot(FittedData$x,log10(FittedData$z/max(FittedData$z)), col = 'red', type = "l")
# lines(FittedData$x,log10(FittedData$y/max(FittedData$y)), col = 'blue', type = "l")

# with fits
malePredictions1 <- readRDS(paste0("Betas/Male/predictions/",Cancer,".rds"))
malePredictions <- unlist(lapply(malePredictions1[[1]], mean))
plot(malePredictions1[[4]],log10(malePredictions1[[5]]/max(malePredictions1[[5]])), col = 'red')
lines(malePredictions1[[4]],log10(malePredictions/max(malePredictions)), col = 'red')

femalePredictions1 <- readRDS(paste0("Betas/Female/predictions/",Cancer,".rds"))
femalePredictions <- unlist(lapply(femalePredictions1[[1]], mean))
points(femalePredictions1[[4]],log10(femalePredictions1[[5]]/max(femalePredictions1[[5]])), col = 'blue')
lines(femalePredictions1[[4]],log10(femalePredictions/max(femalePredictions)), col = 'blue')

# Beta change
plot((unlist(malePredictions1[[3]])/365+20)[-(1:360)], unlist(malePredictions1[[2]])[-(1:360)], type = "l", col='red',
     lwd=2, xlab = "Age", ylab = bquote(beta), main = bquote(.(Cancer)))
lines((unlist(femalePredictions1[[3]])/365+20)[-(1:360)], unlist(femalePredictions1[[2]])[-(1:360)], col='blue', lwd = 2)
legend("bottomright", inset = 0.05, legend=c("Male", "Female"), col=c("red", "blue"), lty=c(1,1), 
       lwd = c(2,2), cex=0.8)

# Test for significance of the interaction
malePredictions1 <- readRDS(paste0("Betas/Male/predictions/",Cancer,".rds"))
betating  <- unlist(malePredictions1[[2]])
femalePredictions1 <- readRDS(paste0("Betas/Female/predictions/",Cancer,".rds"))
betating2  <- unlist(femalePredictions1[[2]])

# Years to prune
prune <- 10*12#months times years
betas <- data_frame(timeStep = c((unlist(malePredictions1[[3]])/365+20)[-(1:prune)],
                                 (unlist(femalePredictions1[[3]])/365+20)[-(1:prune)]), 
                    Beta = c(betating[-(1:prune)],betating2[-(1:prune)]),
                    Gender = c(rep("M", length(betating)-prune),
                               rep("F", length(betating2)-prune)))

ggplot(data = betas, mapping = aes(x = timeStep, y = Beta, color = Gender)) +
  geom_line()+ labs(title = Cancer, x = "Age",y = "Proliferation rate") + 
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))

linearModel <- lm(Beta ~ timeStep*Gender, data=betas) 
summary(linearModel)

koefs <- linearModel$coefficients
(koefs[2]+koefs[4])/koefs[2]

d <- CancerList[1:28]
thu <- split(d, ceiling(seq_along(d)/4))
pdf("all_betas.pdf",title = "Proliferation rate")
for (i in thu) {
  par(mfrow = c(2,2))
  for(Cancer in i){
    # Test for significance of the interaction
    malePredictions1 <- readRDS(paste0("Betas/Male/predictions/",Cancer,".rds"))
    betating  <- unlist(malePredictions1[[2]])
    femalePredictions1 <- readRDS(paste0("Betas/Female/predictions/",Cancer,".rds"))
    betating2  <- unlist(femalePredictions1[[2]])
    
    # Years to prune
    prune <- 30*12#months times years
    betas <- data_frame(timeStep = c((unlist(malePredictions1[[3]])/365+20)[-(1:prune)],
                                     (unlist(femalePredictions1[[3]])/365+20)[-(1:prune)]), 
                        Beta = c(betating[-(1:prune)],betating2[-(1:prune)]),
                        Gender = c(rep("M", length(betating)-prune),
                                   rep("F", length(betating2)-prune)))
    
    p <- ggplot(data = betas, mapping = aes(x = timeStep, y = Beta, color = Gender)) +
      geom_line()+ labs(title = Cancer, x = "Age",y = "Proliferation rate") + 
      theme_bw()+theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()  
