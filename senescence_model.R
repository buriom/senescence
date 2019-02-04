#Model for immunoscenesense via hayflick limit
#Author: Gershom Buri
#Date: 14th Jan 2019
#*******************************************************************************

#******************************  functions **********************************
#function for the fraction of T cells with x divisions at time t; f_t
#experimental parameters: delta = turn over rate (per day), C = initial number
# at t = 0 

# Default parameters
prs <- list(
  delta=0.0015, h=50, n=6, C=0.0017, a=1/(33.33*12), mu=9
)

f <- function(x, t, beta, parms=prs) with(prs, {
   d <- 1 - beta - delta
  #sum of f_0's
  sumj0 <- dpois(x, mu)*(
    1 + C / d * sum((d*exp(a))^(0:-x))
  ) 
  
  coeffs <- d^(-(1:x))*dpois(x - (1:x), mu)
  
  sumj_all <-  sumj0 +  sum(coeffs*sapply(1:x,
    function(j) prod(2*beta*(h^n)/(h^n + (x - (1:j))^n))*(choose(t,j) +  C/d*sum(
      (d*exp(a))^-(0:(x - j))*choose(x - 0:(x - j), j)
    )))
  )
  
    return(d^(t) * total2)
})

