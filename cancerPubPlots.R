# Finalised on 10th of July, this script aims at producing plots for CIM fits to all
# cancers included in the model paper. A number of functions (Convert1Year and GetTick) were borrowed and modified from
# https://github.com/Albluca/ImmuneModelSEER/blob/master/
# 
# 
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

#***************************** Implementation of  CIM**************************

# list of fixed parameters from literature
prs <- list(
  h=50, n=1, a=1/(33.33*365) #Free parameters: tau, mu, beta, delta
)


# Fraction of the T cell population with x no. of divisions at time t; Equantion 3 in pub
f <- function(x, t, tau, mu, beta, delta, fxdParms=prs) with(fxdParms, {
  d <- 1 - beta - delta
  # sum of f_0's
  sumj0 <- dpois(x, mu) + tau / d *dpois(x, mu)* sum((d*exp(a))^(0:x))
  # rest of the sum
  sumj_all <- sumj0+sum(sapply(1:x,function(j) 
    d^(-j)*prod(2*beta*(h^n)/(h^n+(x-(1:j))^n))*
      (choose(t,j)*dpois(x - j, mu) + tau/d*dpois(x - j, mu)*
         sum((d*exp(a))^-(0:(x - j))* choose(x - 0:(x - j), j))))
  )
  return(d^(t) * sumj_all)
})

# Calculating aging whilst updating B_t at each time step (monthly updates)
getPred <- function(mdl, t, parms, fxdParms){
  beta <- parms$delta # assumption that beta == delta at 20 yrs
  mylist1 <- list() 
  for (i in 1:length(t)){
    mylist2 <- c()
    for(j in 1:length(t[[i]])){
      # calculating f for x = 0...100 cell divisions 
      prb <- sapply(0:100, mdl, t = t[[i]][j], parms$tau, parms$mu,
                    beta, parms$delta, fxdParms)
      # updating beta at each time step, i.e after 1 month (Equation 7 of Pub)
      n <- fxdParms$n
      beta <- (parms$delta-parms$tau*exp(-fxdParms$a*t[[i]][j]))/
        (1 - 2 * sum(1/(1 + (fxdParms$h/(0:100))^n) * prb))
      # Calculating immunoscence predictor after normalising the f(x) from above (Equation 6)
      mylist2[j] <- sum((prb/sum(prb))[50:100])
    }
    mylist1[[i]] <- mylist2 
  }
  return(mylist1)
}

# Generate model predictions using best fits and collecting B_t
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

# Process Cancer Data ----------------------------------------------------------

library("SEER2R")

# Load Cancer incidence data from SEER upto 2016
AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic", TXTfileName = "Incidence-Full-1Y-18.txt", UseVarLabelsInData = TRUE)

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

#c("Appendix","Uterus, NOS","Female Genital System","Other Lymphocytic Leukemia")){
for(Cancer in unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))){
  
  if(Cancer %in% c("All Sites","Other Myeloid/Monocytic Leukemia")){
    next()
  }
  
  # Cancer-data processing
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
  
  # Used 20 yrs and above
  ToFil <- YFil == 0 | XFil < 20 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  # get confidence intervals
  DataSubset$Lower_Confidence_Interval <- DataSubset$Lower_Confidence_Interval/ max(YFil)
  DataSubset$Upper_Confidence_Interval <- DataSubset$Upper_Confidence_Interval/ max(YFil)
  DataSubset$AgeAdjusted_Rate <- DataSubset$AgeAdjusted_Rate/ max(YFil)
  # get ticks for the log scale
  Tks <- GetTick(Min = max(c(0.001, min(DataSubset$Lower_Confidence_Interval, na.rm = TRUE))),
                 Max = max(DataSubset$Upper_Confidence_Interval, na.rm = TRUE))
  
  # Scale incidence for fitting
  YFil <- YFil/max(YFil)
  
  skip_to_next <- FALSE
  # read-in model3 fits
  fit <- try(readRDS(paste0("fits/model3fits/",gsub("\\s+","\\",gsub("\\s+","",Cancer)),".rds")[1]))
  
  # use saved CIM predictions or obtain new predictions for each cancer
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
  fitsData <- c(a=a2_m,b=b2_m,RSS=deviance(out2_m),AIC = AIC(out2_m), Rsqr = R2MdlB)
  
  # PLIM fitting 
  myMat <- NULL
  for (i in 1:1000){
    tryCatch({
      out_plim <-
        nls(
          log(YFil) ~ w*log(XFil)+log(a / (exp(
            b * exp(-k_mid * XFil)
          ) - 1)),
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
  #results <- rbind(results,c(Cancer, fit, fitsData, myMat[indx,]))
  
  #--------------------------------- Plotting ------------------------------
  
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
    axis.title=element_text(size=20),
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
  
  p <- p +
    geom_line(
      mapping = aes(x = Converted.Age, y = ModelD.Pred),
      color = adjustcolor("green", alpha.f = 1),
      size = 3
    ) +  geom_line(
    mapping = aes(x = Converted.Age, y = ModelB.Pred),
    color = adjustcolor("purple", alpha.f = 0.6),
    size = 2)+
    geom_line(
      mapping = aes(x = Converted.Age, y = CIMpreds),
      color = adjustcolor("#FF3333", alpha.f = 1),
      size = 3)
  
  p
  # list of plots
  PlotList[[trimws(Cancer)]] <- p
  
}


# Fig 3 A
NewPlotList <- list(
  PlotList$`Cervix Uteri`,
  PlotList$`Thyroid`,
  PlotList$Tonsil,
  PlotList$`Endocrine System`
)

ml=plot_grid(plotlist=NewPlotList, ncol = 2)

ggsave(filename = "weBetterPlots.pdf", plot = ml, width = 8, height = 8,
       scale = 2)
# Fig 3 B
NewPlotList2 <- list(
PlotList$'Salivary Gland',
PlotList$'Colon and Rectum',
PlotList$'Melanoma of the Skin',
PlotList$'Non-Hodgkin Lymphoma')

ml=plot_grid(plotlist=NewPlotList2, ncol = 2)

ggsave(filename = "agreeingPlots.pdf", plot = ml, width = 8, height = 8,
       scale = 2)


# Fig S2
thu <- split(names(PlotList), ceiling(seq_along(names(PlotList))/12))
for(i in 1:length(thu)){
  ml=plot_grid(plotlist=PlotList[thu[[i]]], ncol = 3, nrow = 4)
  # Individual pdf pages that were manually merged there after
  pdfName <-  paste("allCancerPlots",i,".pdf")
  ggsave(filename = pdfName, plot = ml, width = 8, height = 8,
         scale = 2)
}
