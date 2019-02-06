#********************** fit summary ******************************
getPred <- function(t, mdl, beta, parms){
  coll <- sapply(0:100, mdl, t = t*365, beta, parms)
  f_t <- coll/sum(coll)
  return(sum(f_t[50:100]))
}

.args <- c("/fits/fit_simple_brainONS.rds",
           "/model_simple.rds",
           "/preProcessed/brainONS.rds")

.args <- commandArgs(trailingOnly = TRUE)


fit <- readRDS(.args[1])
mdl <- readRDS(.args[2])
incidenceData <- readRDS(.args[3])
t <- fit$tpts
wghts <- fit$wght

prdctns <- sapply(1:length(t), function(i) sum(wghts[[i]] * 
  sapply(1:length(t[[i]]), function(j) getPred((t[[i]][j] - 20), mdl$model, 
  fit$fit$par$beta, mdl$prs)))/(length(wghts[[i]])*sum(wghts[[i]]))) 

title <- sub(".*preProcessed/ *(.*?) *.rds.*", "\\1", .args[3])

jpeg(paste0('figures/', title,'.jpg'))
plot(incidenceData$Age.lb,log10(prdctns*fit$fit$par$M), type = "l", col = "blue",
xlab = "Age", ylab = expression('log'[10]*' Incidence'),main =title)

#lines(df$Age, log10(prdctns*fit$par$M), col = "blue")

#plot(incidenceData$Age.ub, log10(prdctns*fit$par$M), t)
points(incidenceData$Age.lb, log10(incidenceData$`Rate per 100,000`))
least <- min(log10(prdctns*fit$fit$par$M))
Rsquare <- cor(log10(incidenceData$`Rate per 100,000`),log10(prdctns*fit$fit$par$M))

text(70, (least + 0.2), cex = 1, bquote(paste('R'^2*' = ', .(Rsquare^2))))
dev.off()
