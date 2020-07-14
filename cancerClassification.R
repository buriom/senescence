library(ggplot2)
library("SEER2R")

AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic",
                         TXTfileName = "Incidence-Full-1Y-18.txt",
                         UseVarLabelsInData = TRUE)

allCancers <- unique(trimws(AllData$`Site_recode_ICDO3/WHO_2008`))

d <- sort(allCancers)
d <- d[-c(4,37,48,51,61,66,69)]
d[47] <- "Male Genital System"

sumry <- function(i){
  fit <- readRDS(paste0("fits/model3fits/",gsub("\\s+","\\",d[i]),".rds")[1])
  return(fit)
  }

resultSumry <- sapply(1:95, sumry)
y <- as.data.frame(t(resultSumry))
rownames(y) <- d
#View(y)
pca <- prcomp(y, scale=FALSE)
km <- kmeans(scale(y[,1:3]),5)
km2 <- kmeans(scale(y[,4:5]),1)
y$`parm clstrs` <- as.factor(km$cluster)
y$`fit clstrs` <- as.factor(km2$cluster)

stdData1 <- cbind(y,pca$x[,1:2])

ggplot(as.data.frame(stdData1), aes(x=PC1, y=PC2, shape=`fit clstrs`, color = `parm clstrs`)) + 
  geom_jitter(position = "jitter",size=3) + theme_bw() + 
  #ggtitle('First 2 prinipal components clustering w/ kmeans') + 
  xlab("Principal Component 1") + ylab("Principal Component 2")
rownames(stdData2)[which(stdData1[,"fit clstrs"] %in% c(1))]

rownames(stdData2)[which(stdData1[,"fit clstrs"] %in% c(2)  & stdData1[,"parm clstrs"] == 2)]
