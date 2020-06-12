library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

# Support Functions ------------------------------------------------

# format year column
Convert1Year <- function(Vect) {
  as.integer(trimws(strsplit(
    x = Vect, split = "years", fixed = TRUE
  )))
}

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

# Predictions for the CIM model

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

# Process Infection Data --------------------------------------------------




setwd("~/Downloads")

library(xlsx)


Dis1 <- read.xlsx("Diseases2.xlsx", sheetIndex = 1,stringsAsFactors=FALSE)

Dis1$Age=as.numeric(as.character(Dis1$Age))
Dis1$Inc=as.numeric(as.character(Dis1$Inc))
Dis1$LowerCI=as.numeric(as.character(Dis1$LowerCI))
Dis1$UpperCI=as.numeric(as.character(Dis1$UpperCI))
Dis1$Pop=as.numeric(as.character(Dis1$Pop))


InfToPlot=c(  
  # "MRSA"           ,             "West Nile Virus Disease"   ,           
  # "Streptococcus pneumoniae"   ,  "group B Streptococcus "    
  # ,    "Legionellosis"                 
  # ,"Haemophilus influenzae "  ,
)






# Process Cancer Data -----------------------------------------------------

library("SEER2R")

AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic", TXTfileName = "Incidence-Full-1Y-18.txt", UseVarLabelsInData = TRUE)

# unique(AllData$`Site_recode_ICDO3/WHO_2008`)

AllData$AgeAdjusted_Rate <- as.numeric(levels(AllData$AgeAdjusted_Rate))[AllData$AgeAdjusted_Rate]
AllData$Lower_Confidence_Interval <- as.numeric(levels(AllData$Lower_Confidence_Interval))[AllData$Lower_Confidence_Interval]
AllData$Upper_Confidence_Interval <- as.numeric(levels(AllData$Upper_Confidence_Interval))[AllData$Upper_Confidence_Interval]

# Plot all Cancers --------------------------------------------------------------

k_mid=0.044

# to read-in previously computed CIM predictions
prevPreds <- readRDS("predictionsComparison.rds")

# to hold parameter fits
results <- NULL
# to hold plots
PlotList <- list()

#for(Cancer in c("Appendix","Uterus, NOS","Female Genital System","Other Lymphocytic Leukemia"))
for(Cancer in unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))){
  
  if(Cancer %in% c("All Sites","Other Myeloid/Monocytic Leukemia")){
    next()
  }
  
  
  DataSubset <- AllData[trimws(AllData$`Site_recode_ICDO3/WHO_2008`) == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "Male and female",]
  dim(DataSubset)
  
  print(Cancer)
  
  Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
  Converted.Age[Converted.Age==84] <- NA
  
  DataSubset <- cbind(DataSubset, Converted.Age)
  
  XFil <- DataSubset$Converted.Age
  YFil <- DataSubset$AgeAdjusted_Rate
  
  DataSubset <- DataSubset[-which(is.na(XFil)),]
  
  
  ToFil <- YFil == 0 | XFil < 18 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  
  
  DataSubset$Lower_Confidence_Interval <- DataSubset$Lower_Confidence_Interval/ max(YFil)
  DataSubset$Upper_Confidence_Interval <- DataSubset$Upper_Confidence_Interval/ max(YFil)
  DataSubset$AgeAdjusted_Rate <- DataSubset$AgeAdjusted_Rate/ max(YFil)
  
  Tks <- GetTick(Min = max(c(0.001, min(DataSubset$Lower_Confidence_Interval, na.rm = TRUE))),
                 Max = max(DataSubset$Upper_Confidence_Interval, na.rm = TRUE))
  
  YFil <- YFil/max(YFil)
  
  #FittedData <- data.frame(cbind(XFil, YFil))
  #colnames(FittedData) <- c("x", "y")
  skip_to_next <- FALSE
  # read-in model3 fits
  fit <- try(readRDS(paste0("fits/model3fits/",gsub("\\s+","\\",gsub("\\s+","",Cancer)),".rds")[1]))
  
  if(Cancer %in% names(prevPreds)){
    CIMpreds1 <- prevPreds[[Cancer]]
    CIMpreds <- c(rep(NA,dim(DataSubset)[1]-length(CIMpreds1)),CIMpreds1)
    DataSubset <- cbind(DataSubset, CIMpreds)
  }else{
    FittedData <- readRDS(paste0("preProcessed/", gsub("\\s+","\\",Cancer), "_data.rds")[1])
    ageColumn <- FittedData$x
    lu <- ageColumn -1
    FittedData$lb <- (lu-lu[1])*365+365/12
    FittedData$ub <- (ageColumn-lu[1])*365
    tpts <-  mapply(function(i,j) seq(i,j,365/12), FittedData$lb, FittedData$ub,
                    SIMPLIFY = FALSE)
    
    fit1 <- as.list(fit[1:3])
    #extract CIM  predictions
    prdctns1 <-  bestPrdctns(f, tpts, fit1, prs)
    CIMpreds1 <- unlist(lapply(prdctns1[[1]], mean))
    CIMpreds <- c(rep(NA,dim(DataSubset)[1]-length(CIMpreds1)),CIMpreds1)
    DataSubset <- cbind(DataSubset, CIMpreds)
  }
  
  # IM-II fitting 
  out2_m <-
    nls(
      log(YFil) ~ log(a / (exp(
        b * exp(-k_mid * XFil)
      ) - 1)),
      start = list(a = 1, b = 1),
      algorithm = "port",
      lower = c(1e-10, 1e-10),
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  a2_m <- summary(out2_m)$parameters["a", "Estimate"]
  b2_m <- summary(out2_m)$parameters["b", "Estimate"]

  for (i in 1:10) {
    out2_m <-
      nls(
        log(YFil) ~ log(a / (exp(
          b * exp(-k_mid * XFil)
        ) - 1)),
        start = list(a = a2_m, b = b2_m),
        lower = c(1e-10, 1e-10),
        algorithm = "port",
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8/i,
          warnOnly = TRUE
        )
      )
    a2_m <- summary(out2_m)$parameters["a", "Estimate"]
    b2_m <- summary(out2_m)$parameters["b", "Estimate"]
  }
  
  # IM-II predictions
  ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
  ModelB.Pred[DataSubset$Converted.Age<head(XFil,1) | DataSubset$Converted.Age>83] <- NA
  DataSubset <- cbind(DataSubset, ModelB.Pred)
  
  R2MdlB <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
  fitsData <- c(a=a2_m,b=b2_m,RSS=deviance(out_plim),AIC = AIC(out2_m), Rsqr = R2MdlB)
  
  # PLIM fitting 
  myMat <- NULL
  for (i in 1:500){
    tryCatch({
      out_plim <-
        nls(
          log(YFil) ~ w*log(XFil)+log(a / (exp(
            b * exp(-k_mid * XFil)
          ) - 1)),
          #data = FittedData,
          start = list(a = runif(1), b = runif(1,0,13),w = runif(1)),
          algorithm = "port",
          lower = rep(1e-11, times=3),
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8,
            warnOnly = TRUE))
      R2MdlD <- 1 - sum((summary(out_plim)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
      myMat <- rbind(myMat,c(summary(out_plim)$coefficients[,1], RSS = deviance(out_plim), AIC = AIC(out_plim), Rsqr = R2MdlD))
    },
    error = function(e) { skip_to_next <- TRUE})
    if(skip_to_next) { next }
    
  }
  
  indx <-  which(myMat[,4] == min(myMat[,4]))[1]
  myMat[indx,]
  a_plim <- myMat[indx,1]
  b_plim <- myMat[indx,2]
  w_plim <- myMat[indx,3]

  # PLIM predictions 
  
  ModelD.Pred <- DataSubset$Converted.Age^w_plim*(a_plim / (exp(b_plim * exp(-k_mid * DataSubset$Converted.Age)) - 1))
  ModelD.Pred[DataSubset$Converted.Age<head(XFil,1) | DataSubset$Converted.Age>83] <- NA
  DataSubset <- cbind(DataSubset, ModelD.Pred)
  
  # matrix of model fits and measures of goodness of fits
  results <- rbind(results,c(Cancer, fit, fitsData, myMat[indx1,]))

  BgkCol <- "white"
  
  p <- ggplot(
    data = DataSubset,
    mapping = aes(
      x = Converted.Age,
      y = AgeAdjusted_Rate
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
    ggtitle(trimws(Cancer)) +
    geom_point(color="blue")
  
  p <- p + scale_y_log10(
    breaks = Tks$Ticks,
    labels = Tks$TickLab,
    #limits = range(Tks$Ticks)
  ) + scale_x_continuous(
    breaks = c(0, 45, 90),
    limits = c(0, 90)
  )
  
  p <- p  + geom_pointrange(color="blue",
                            aes(
                              ymin = Lower_Confidence_Interval,
                              ymax = Upper_Confidence_Interval
                            )
  )
  
  p <- p +  geom_line(
    mapping = aes(x = Converted.Age, y = ModelB.Pred),
    color = adjustcolor("green", alpha.f = 1),
    size = 3)+
    geom_line(
      mapping = aes(x = Converted.Age, y = CIMpreds),
      color = adjustcolor("red", alpha.f = 1),
      size = 3)+
    geom_line(
      mapping = aes(x = Converted.Age, y = ModelD.Pred),
      color = adjustcolor("purple", alpha.f = 0.8),
      size = 2
    )
  
  p
  # list of plots
  PlotList[[trimws(Cancer)]] <- p
  
}


NewPlotList <- list(
  PlotList$`Chronic Myeloid Leukemia`,
  PlotList$`Soft Tissue including Heart`,
  PlotList$Brain,
  PlotList$`Colon and Rectum`,
  PlotList$Lip,
  PlotList$Stomach)

ml=plot_grid(plotlist=NewPlotList, ncol = 3)

ggsave(filename = "Fig3.pdf", plot = ml, width = 8,  scale = 2)
