library(data.table)

.args <- c("cancerData/BonesandJointsCancer_proc.csv", "BonesandJointsCancer_data.rds")
.args <- commandArgs(trailingOnly = TRUE)


#********************************** data load ********************************
df <- fread(.args[1], na.strings = "-")

df <- df[Sex == "Both Sexes" & `Race/Ethnicity` == "All Races (includes Hispanic)" & Age != "All Ages"]

# df <- df[,.SD,.SDcols=c("Age", "Year", grep("per 100,000", colnames(df), value = T)) ]
# setnames(df, c("Age","Year","Rate per 100,000"))
# I deleted the year column
df <- df[,.SD,.SDcols=c("Age", grep("per 100,000", colnames(df), value = T)) ]
setnames(df, c("Age","Rate per 100,000"))
# this gives us a sample, but might want to ensure it's detrended for our purposes

df[Age=="All Ages", Age := "0-89" ]
df[, Age := gsub("^[[:alpha:][:blank:]]+", "", Age) ]

df <- df[,  Age := gsub(".*<[[:blank:]]*",   "0-",   Age)][,  Age := gsub("\\+.*$", "-89", Age)][ Age != "" ]

df[,c("Age.lb","Age.ub") := {
  .(as.integer(sub("-.+", "", Age)), as.integer(sub(".*-", "", Age)))
}]

#remove NA rows
df <- df[
  Age.lb >= 20 & !is.na(`Rate per 100,000`),
  .(`Rate per 100,000`=mean(`Rate per 100,000`)),
  by=.(Age.lb, Age.ub)
]

#replace age ranges by mid-points
#df$Age <- seq(22,90,5)[(15 -(dim(df)[1])):14]#mid points no the best index for
#the ranges
# View(df)
saveRDS(df, tail(.args,1))
