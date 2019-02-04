library(data.table)


.args <- c("/home/buri/cancerData/lip.csv", "lip.rds")
.args <- commandArgs(trailingOnly = TRUE)


#********************************** data load ********************************
df <- fread(.args[1], na.strings = "-")
df <- df[,c("Age","Rate per 100,000")]
df <- within(df, {
  Age <- gsub("<",   "0-",   Age)
  Age <- gsub("\\+", "-100", Age)
})

#str_split_fixed(df$Age,"-",2)
df <- within(df, {
  Age.lb <- as.integer(sub("-.+", "", Age))
  Age.ub <- as.integer(sub(".*-", "", Age))
})

#remove NA rows
df <- df[
  Age.ub >= 20 &
  !is.na(`Rate per 100,000`)
]

#replace age ranges by mid-points
#df$Age <- seq(22,90,5)[(15 -(dim(df)[1])):14]#mid points no the best index for
#the ranges
# View(df)
saveRDS(df, tail(.args,1))
