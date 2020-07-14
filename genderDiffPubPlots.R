library(dplyr)
# To gather gender-specific proliferation rates of 31 best fits
noGenderDifference <- c("Female Genital System",
                        "Cervix Uteri",
                        "Corpus and Uterus, NOS",
                        "Corpus Uteri",
                        "Uterus, NOS",
                        "Ovary",
                        "Vagina",
                        "Vulva",
                        "Other Female Genital Organs",
                        "Male Genital System",
                        "Prostate",
                        "Testis",
                        "Penis",
                        "Other Male Genital Organs")

bestFits <- read.csv("bestFits.csv", header = F)
d <-  system("ls Betas/Female/predictions/", TRUE)
alreadyCalcd <- sapply(d,function(x) substr(x, 1, nchar(x)-4))
cleanBestFits <- bestFits %>% filter(!(V1 %in% noGenderDifference))%>% 
  filter((V1 %in% alreadyCalcd))
resToFit <- (as.character(cleanBestFits$V1))

thu <- split(resToFit[-c(14,22)], ceiling(seq_along(resToFit[-c(15,24)])/4))


pdf("annotated_betas.pdf",title = "Proliferation rate")

my_ratios <- list()
my_signf <- list()
for (i in thu[[1]]) {
  df <- data.frame()
  for(Cancer in i){
    # if (Cancer %in% c("Skin excluding Basal and Squamous","Lip", "Esophagus", "Tonsil")){
    #   next()
    # }
    # Test for significance of the interaction
    malePredictions1 <- readRDS(paste0("Betas/Male/predictions/",Cancer,".rds"))
    betating  <- unlist(malePredictions1[[2]])
    femalePredictions1 <- readRDS(paste0("Betas/Female/predictions/",Cancer,".rds"))
    betating2  <- unlist(femalePredictions1[[2]])
    
    # Years to prune
    prune <- 20*12#months times years
    betas <- data_frame(timeStep = c((unlist(malePredictions1[[3]])/365+20)[-(1:prune)],
                                     (unlist(femalePredictions1[[3]])/365+20)[-(1:prune)]), 
                        Beta = c(betating[-(1:prune)],betating2[-(1:prune)]),
                        Gender = c(rep("M", length(betating)-prune),
                                   rep("F", length(betating2)-prune)),
                        cancer=Cancer)
    linearModel <- lm(Beta ~ timeStep*Gender, data=betas)
    smry <- summary(linearModel)
    koefs <- smry$coefficients
    my_ratios[Cancer] <- (koefs[2,1]+koefs[4,1])/koefs[2,1]
    my_signf[Cancer] <- koefs[4,4]
    df <- rbind.data.frame(df,betas)
  }

  p <- ggplot(data = df, mapping = aes(x = timeStep, y = Beta, color = Gender)) +
    geom_line()+ labs(title = "Gender seperated proliferation", x = "Age",y = "Proliferation rate") +
    facet_wrap( ~  cancer, ncol=2,nrow = 2) + 
    #theme_pander() +
    #theme(legend.position="none",) + 
    theme_bw()+theme(plot.title = element_text(hjust = 0.5),strip.background = element_rect( fill="white"))
  print(p)
}
        
dev.off()  

summary((unlist(my_ratios)))

plotGraph <- function(resToFit){
  for(Cancer in resToFit){
    # if (Cancer %in% c("Skin excluding Basal and Squamous","Lip", "Esophagus", "Tonsil")){
    #   next()
    # }
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
    # linearModel <- lm(Beta ~ timeStep*Gender, data=betas) 
    # koefs <- linearModel$coefficients
    # my_ratios[Cancer] <- (koefs[2]+koefs[4])/koefs[2]
    p <- ggplot(data = betas, mapping = aes(x = timeStep, y = Beta, color = Gender)) +
      geom_line()+ labs(title = Cancer, x = "Age",y = "Proliferation rate") + 
      theme_bw()+theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
