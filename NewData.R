library("SEER2R")
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

allCancers <- unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))


for(Cancer in allCancers){
  #Cancer <- "Anus, Anal Canal and Anorectum"
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
  
  ToFil <- YFil == 0 | XFil < 18 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  
  FittedData <- data.frame(cbind(XFil, YFil))
  colnames(FittedData) <- c("x", "y")
  saveRDS(FittedData, paste0("Preprocessed/", gsub("\\s+","\\",Cancer),"_data.rds"))
}
