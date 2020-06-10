library(dplyr)
library(ggplot2)
library(xlsx)
library(SEER2R)


k_mid=0.044


# Define support functions

Convert1Year <- function(Vect) {
  as.integer(trimws(strsplit(
    x = Vect, split = "years", fixed = TRUE
  )))
}



Fit2ParModel <- function(FittedData, a=NULL, b=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(k_mid))]
  LowVect <- rep(-Inf, times=length(StartList))
  
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


Fit1ParModel <- function(FittedData, a=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
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


Fit3ParModel <- function(FittedData, a=NULL, b=NULL, c=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1),c = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(c))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out2_m <-
    nls(
      log(y) ~ c*log(x)+log(a / (exp(
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
          log(y) ~ c*log(x)+log(a / (exp(
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
  
  PVals <- c(NA, NA, NA, NA)
  names(PVals) <- c("a", "b","c", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, b, c, k_mid)
  
  ReturnVect <- c(PVals, AIC(out2_m), R2Mb, deviance(out2_m), sum(summary(out2_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "c", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}


FitPLModel <- function(FittedData, Rounds=10){
  
  out_pl <-
    nls(
      log(y) ~ a*log(x) + b,
      start = list(a = 1, b = 1),
      data = FittedData,
      algorithm = "port",
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  a_pl <- summary(out_pl)$parameters["a", "Estimate"]
  b_pl <- summary(out_pl)$parameters["b", "Estimate"]
  
  
  for (i in 1:Rounds) {
    out_pl <-
      nls(
        log(y) ~ a*log(x) + b,
        data = FittedData,
        start = list(a = a_pl, b = b_pl),
        algorithm = "port",
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8/i,
          warnOnly = TRUE
        )
      )
    a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    
  }
  
  RPLMb <- 1 - sum((summary(out_pl)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  ReturnVect <- c(a_pl, b_pl, AIC(out_pl), RPLMb, deviance(out_pl), sum(summary(out_pl)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}

TryFitting <- function(Model, FittedData, a = NULL, b = NULL, k_mid = NULL,
                       t_0 = NULL, t_0_low = -Inf, t_0_high = Inf,
                       Rounds = 20, Retires = 10){
  
  AllFits <- NULL
  
  if(Model == "A"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit1ParModel(FittedData = FittedData, a = a, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "B"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit2ParModel(FittedData = FittedData, a = a, b = b, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "C"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- FitPLModel(FittedData = FittedData, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  if(Model == "PLIM"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit3ParModel(FittedData = FittedData, a = NULL, b = NULL,c = NULL, k_mid = NULL, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "TM"){
    for(i in 1:Retires){
      if(is.null(k_mid)){
        TR <- try(ModelPars <- FitTimeModel(FittedData = FittedData, a = a, t_0 = t_0, k = NULL,
                                            t_0_low = t_0_low, t_0_high = t_0_high, Rounds = Rounds))
      } else{
        TR <- try(ModelPars <- FitTimeModel(FittedData = FittedData, a = a, t_0 = t_0, k = -k_mid,
                                            t_0_low = t_0_low, t_0_high = t_0_high, Rounds = Rounds))
      }
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  return(AllFits[which.min(AllFits[,"Resid"]),])
  
}

AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic",
                         TXTfileName = "Incidence-Full-1Y-18.txt",
                         UseVarLabelsInData = TRUE)
AllData$AgeAdjusted_Rate <- as.numeric(levels(AllData$AgeAdjusted_Rate))[AllData$AgeAdjusted_Rate]
AllData$Lower_Confidence_Interval <- as.numeric(levels(AllData$Lower_Confidence_Interval))[AllData$Lower_Confidence_Interval]
AllData$Upper_Confidence_Interval <- as.numeric(levels(AllData$Upper_Confidence_Interval))[AllData$Upper_Confidence_Interval]


SumData <- NULL

for(Cancer in unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))){
  
  if(Cancer == "All Sites"){
    next()
  }
  
  print(Cancer)
  
  DataSubset <- AllData[trimws(AllData$`Site_recode_ICDO3/WHO_2008`) == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "Male and female",]
  
  if(sum(DataSubset$AgeAdjusted_Rate, na.rm = TRUE)==0){
    next
  }
  
  Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
  Converted.Age[Converted.Age==84] <- NA
  
  DataSubset <- cbind(DataSubset, Converted.Age)
  
  XFil <- DataSubset$Converted.Age
  YFil <- DataSubset$AgeAdjusted_Rate
  
  #    Tks <-
  #      GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
  #              Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))
  
  ToFil <- YFil == 0 | XFil < 20 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  
  FittedData <- data.frame(cbind(XFil, YFil))
  colnames(FittedData) <- c("x", "y")
  # rescaled this way because min-max would result in log(o) = Inf
  #FittedData$y <-  FittedData$y/ max(FittedData$y)
  
  # ModData.A <- TryFitting(Model = "A", FittedData = FittedData,k_mid = k_mid,
  #                         Rounds = 4, Retires = 50)
  # 
  # ModData.B <- TryFitting(Model = "B", FittedData = FittedData,k_mid = k_mid,
  #                         Rounds = 9, Retires = 100)
  # ModData.C <- TryFitting(Model = "C", FittedData = FittedData,
  #                         Rounds = 9, Retires = 100)
  ModData.D <- TryFitting(Model = "PLIM", FittedData = FittedData,
                          Rounds = 9, Retires = 10)
  
  
  
  # names(ModData.A) <- paste("IM-I",names(ModData.A), sep = ":")
  # names(ModData.B) <- paste("IM-II",names(ModData.B), sep = ":")
  # names(ModData.C) <- paste("PLM",names(ModData.C), sep = ":")
  names(ModData.D) <- paste("PLIM",names(ModData.D), sep = ":")
  SumData <- rbind(SumData, c(Cancer, ModData.D))#,  ModData.A,ModData.B, ModData.C))
  
  
}

summary(bestFits[bestFits$Cancer.Site %in% fit3names,7])
sum(bestFits[bestFits$Cancer.Site %in% fit3names,7]>0.95)
summary(as.numeric(SumData3[,7]))
sum(as.numeric(SumData3[,7]>0.95))

SumData[,1] <- trimws(SumData[,1])
colnames(SumData)[1] <- "Cancer/Site"


AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic",
                         TXTfileName = "Incidence-Full-1Y-18.txt",
                         UseVarLabelsInData = TRUE)

allCancers <- unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))
  #d <- sort(allCancers)[c(4,37,48,51,61,66,69)]
#allCancers[100] <- "Male Genital System"

sumData2 <- NULL
for(Cancer in allCancers[-c(96,100)]){
  #load the data and calculate the population weights in each class
  #fit <- try(readRDS(paste0("fits/model3fits/",gsub("\\s+","\\",Cancer),".rds")[1]))# model3 fits
  fit <- try(readRDS(paste0("fits/model1fits/",gsub("\\s+","\\",Cancer),".rds")[1])) # model1 fits
  
  names(fit) <- paste("Sene-Mdl",names(fit), sep = ":")
  sumData2 <- rbind(sumData2, c(Cancer, fit))
}
colnames(sumData2)[1] <- "Cancer.Site"


cmplteData <- left_join(as.data.frame(sumData2),as.data.frame(SumData), "Cancer.Site")
finalData <- cmplteData[-4,]
write.csv(finalData,"resultsTable.csv")
