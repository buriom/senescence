# Finalised on 10th of July, this script contains all necessary functions for fitting and plotting
# CIM to T-lymphoblastic leukemia data
# A number of functions were borrowed and modified from  https://github.com/Albluca/ImmuneModelSEER/blob/master/

library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

# Support Functions ------------------------------------------------

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

#list of fixed parameters from literature
prs <- list(
  h=50, n=1, a=1/(33.33*365)#Free parameters: tau, mu, beta, delta
)

# Fraction of the T cell population with x no. of divisions at time t; Equantion 3 in pub
f <- function(x, t, tau, mu, beta, delta, fxdParms=prs) with(fxdParms, {
  d <- 1 - beta - delta
  #sum of f_0's
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
  beta <- parms$delta# assumption that beta == delta at 20 yrs
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

residFun <- function(par, mdl,  observed, t, wghts, fxdParms){
prdctns <-  getPred(mdl, t, par, fxdParms)
   resids <- log(unlist(lapply(prdctns, mean))) - log(observed)
   return(ifelse(is.nan(resids), 1e6, resids))
}

#--------------------------------- data fitting ------------------------------
# Load previous fits
#fit <- readRDS("fits/model3fits/T cell Lymphoblastic Leukemia.rds")

# load data for fitting
FittedData <- read.csv("T_Lymphoblastic.csv", header = TRUE)

# adapt the time points for monthly updates of B_t
jug <- FittedData$x[-(1:20)]
lu <- jug -1
lb <- (lu-lu[1])*365+365/12
ub <- (jug-lu[1])*365

tpts <-  mapply(function(i,j) seq(i,j,365/12), lb, ub, 
                SIMPLIFY = FALSE)

# scale incidence rates
obs <- FittedData$y/ max(FittedData$y)

# #initial values and parameter bounds for fitting
inits <- list(mu=42.7372, tau=5.00578e-05 , delta = 1e-07); a <- c(0,1e-5,1e-5); b <- c(100,1, 1)

# fitting after normalising
system.time(fit <- nls.lm(par = inits, lower = a, upper = b,
                          fn = residFun, mdl = f, fxdParms = prs, observed = obs[-(1:20)], t = tpts,
                          control = nls.lm.control(nprint=1,  ftol = 1e-1)))

#extract CIM  predictions
prdctns1 <-  bestPrdctns(f, tpts, fit$par, prs)
CIMpreds1 <- unlist(lapply(prdctns1[[1]], mean))
CIMpreds <- c(rep(NA,sum(FittedData$x<20)),CIMpreds1)

FittedData$CIMpreds <- CIMpreds

# scale the confidence intervals and incidence rate
FittedData$Lower.Confidence.Interval <- FittedData$Lower.Confidence.Interval/max(FittedData$y)

FittedData$Upper.Confidence.Interval <- FittedData$Upper.Confidence.Interval/max(FittedData$y)

FittedData$y <- FittedData$y/max(FittedData$y)

# Get ticks for the logarithm scale
Tks <- GetTick(Min = max(c(0.001, min(FittedData$Lower.Confidence.Interval, na.rm = TRUE))),
               Max = max(FittedData$Upper.Confidence.Interval, na.rm = TRUE))

#--------------------------------- Plotting ------------------------------
BgkCol <- "white"

p <- ggplot(
  data = FittedData,
  mapping = aes(
    x = x,
    y = y
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
  ggtitle(trimws("T-Lymphoblastic Leukemia")) +
  geom_point(color="blue")

p <- p + scale_y_log10(
  breaks = Tks$Ticks,
  labels = Tks$TickLab,
  limits = range(Tks$Ticks)
) + scale_x_continuous(
  breaks = c(0, 45, 90),
  limits = c(0, 90)
)

p <- p  + geom_pointrange(color="blue",
                          aes(
                            ymin = Lower.Confidence.Interval,
                            ymax = Upper.Confidence.Interval
                          ))

p <- p +
  geom_line(
    aes(x = x, y = CIMpreds),
    color = adjustcolor("#FF3333", alpha.f = 1),
    size = 3)

p

# Save plot
ggsave(filename = "T-Lymphoblastic Leukemia.pdf", plot = p, width = 4, height = 4,
       scale = 2)

