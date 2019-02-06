library("minpack.lm")

#******************************** Calculating aging ************************
getPred <- function(t, mdl, beta, parms){
  #Predictions at time point t
  coll <- sapply(0:100, mdl, t = t*365, beta, parms)
  #normalising
  f_t <- coll/sum(coll)
  #summing f_t >= h
  return(sum(f_t[50:100]))
}

#********************************* Residual function ************************
residFun <- function(par, mdl,  observed, t, wghts, parms){
  prdctns <- sapply(1:length(t), function(i) sum(wghts[[i]] * 
  sapply(1:length(t[[i]]), function(j) getPred((t[[i]][j] - 20), mdl, par$beta, parms))
                    )/(length(wghts[[i]])*sum(wghts[[i]])))   
  resids <- log10(par$M * prdctns) - log10(observed)
  return(ifelse(is.nan(resids), 1e6, resids))
}

#********************************* data fitting ******************************
.args <- c("/home/buri/senescenceModel/model_simple.rds", 
           "/home/buri/senescenceModel/preProcessed/lip.rds","fit1.rds")

.args <- commandArgs(trailingOnly = TRUE)

#model <- readRDS(.args[1])
incidenceData <- readRDS(.args[2])
#prs <- model$defPars
tpts <-  mapply(function(i,j) i:j, incidenceData$Age.lb, incidenceData$Age.ub)
wght <-  mapply(function(i,j) {things <- i:j; ifelse( things<50, 1, 4-0.04*(things))},
                 incidenceData$Age.lb, incidenceData$Age.ub)

fit <- with(readRDS(.args[1]), {
  nls.lm(par = inits, lower = L, upper = U,
    fn = residFun, mdl = model, parms = prs, observed = incidenceData$`Rate per 100,000`,
    t = tpts, wghts = wght)
  })#, control = nls.lm.control(nprint=1)

fit_results <- list(
  fit = fit,
  tpts = tpts,
  wght = wght
)
saveRDS(fit_results, tail(.args,1))