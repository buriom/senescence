library(ggplot2)
library(reshape2)
# Hill function
hillFun <- function(x,n,threshold=50) x^n/(threshold^n+x^n)
plot(hillFun(1:100,3,30),type = "l")


#list of fixed parameters from literature
prs <- list(
  delta = 0.0015, h=50, n=3, tau=0.0017, a=1/(33.33*365)#fixed paramaters
)

#proportion of the T cell population with x no. of divisions at time t
f <- function(x, t, beta, mu, parms=prs) with(parms, {
  d <- 1 - beta - delta
  #sum of f_0's
  sumj0 <- dpois(x, mu) + tau / d *dpois(x, mu)* sum((d*exp(a))^(0:-x))
  #  coeffs <- d^(-(1:x))
  sumj_all <-  sumj0 +  sum(sapply(1:x,function(j) 
    d^(-j)*prod(2*beta*(h^n)/(h^n+(x - (1:j))^n))*
      (choose(t,j)*dpois(x - j, mu)+tau/d*dpois(x - j, mu)*
         sum((d*exp(a))^-(0:(x - j))*choose(x - 0:(x - j), j))))
  )
  return(d^(t) * sumj_all)
})

thing <- sapply(0:100, f,7*365,0.008,34,prs)
plot(thing, type = 'l')
# normalised distro
plot(thing/sum(thing), type = 'l')

# Generate theoretical distro
generations <- seq(0,60,2)
thing <- sapply(generations, function(x) sapply(0:100, f,x*365,0.0037,10,prs))
thing3 <- t(thing)/apply(thing, 2, sum)
thing2 <- as.data.frame(t(thing3))
colnames(thing2) <- 20+ generations
thing2$Divisions <- 0:100

#tito <- "Evolution of  through time"
df <- melt(thing2 ,  id.vars = 'Divisions', variable.name = 'years')
#colnames(df) <- c('Divisions', 'years', 'value')
ggplot(df, aes(Divisions,value)) + geom_line(aes(colour = years))
