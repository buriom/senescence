library(tidyverse)


.args <- c("/home/buri/cancerData/lip.csv", "lip.rds")
.args <- commandArgs(trailingOnly = TRUE)


#********************************** data load ********************************
df <- read.csv(.args[1], na.strings = "-")
df <- df[,c("Age","Rate.per.100.000")]

df1 <- as.tibble(df)

#remove NA rows
df <- df[!is.na(df$Rate.per.100.000),]
#replace age ranges by mid-points
df$Age <- seq(22,90,5)[(15 -(dim(df)[1])):14]#mid points no the best index for
#the ranges
View(df)
saveRDS(df, tail(.args,1))