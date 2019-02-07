getPred <- function(t, mdl, beta, parms){
  coll <- sapply(0:100, mdl, t = t*365, beta, parms)
  f_t <- coll/sum(coll)
  return(sum(f_t[50:100]))
}

#********************** fit summary ******************************
.args <- c("/fits/fit_simple_brainONS.rds",
           "/model_simple.rds",
           "/preProcessed/brainONS.rds")

.args <- commandArgs(trailingOnly = TRUE)

#argument assignment
fit <- readRDS(.args[1])
mdl <- readRDS(.args[2])
incidenceData <- readRDS(.args[3])
t <- fit$tpts
wghts <- fit$wght

#fitted model predicitions
prdctns <- sapply(1:length(t), 
  function(i) sum(wghts[[i]] * sapply(1:length(t[[i]]),
  function(j) getPred((t[[i]][j] - 20), mdl$model, fit$fit$par$beta, mdl$prs)))/
    (length(wghts[[i]])*sum(wghts[[i]]))) 

#********************************** Plotting ********************************
title <- sub(".*preProcessed/ *(.*?) *.rds.*", "\\1", .args[3])

jpeg(paste0('figures/', title,'.jpg'))
age <- incidenceData$Age.lb
preds <- log10(prdctns*fit$fit$par$M)
obs <- log10(incidenceData$`Rate per 100,000`)

plot(age, preds, type = "l", col = "blue",
xlab = "Age", ylab = expression('log'[10]*' Incidence'),main =title)

points(age, obs)
least <- min(preds)
Rsquare <- cor(obs,preds)

text(70, (least + 0.2), cex = 1, bquote(paste('R'^2*' = ', .(Rsquare^2))))
text(70, (least + 0.4), cex = 1, bquote(paste('Beta = ', .(fit$fit$par$beta))))

dev.off()
