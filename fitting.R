
#********************************* data fitting ******************************

fit <- nls.lm(par = inits, lower = L, upper = U,
              fn = residFun, observed = incidenceData$Rate.per.100.000,
              t = incidenceData$Age,
              control = nls.lm.control(nprint=1))#,  ftol = 1e-4
