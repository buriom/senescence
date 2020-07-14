# In this script, we generate data for Fig 1 in the pub. The final plot however is 
# generated in python
library(ggplot2)
library(reshape2)

# Hill function
hillFun <- function(x,n,HL=50) x^n/(HL^n+x^n)
n=3
HL=55
plot(hillFun(1:100,n,HL),type = "l", col = "red", lwd = 3,
     main = paste("Hill function with n=",n, "and Hayflick limit = ", HL),
     xlab = "Hayflick Limit", ylab = "Probability of replicative immunosenescence")


#list of fixed parameters from literature
prs <- list(
  delta = 0.0015, h=50, n=3, tau=0.0017, a=1/(33.33*365)
)

# Fraction of the T cell population with x no. of divisions at time t, Pub Eqn 3
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

# Plot of No. of divisions distribution
divisionsDistro <- sapply(0:100, f,7*365,0.008,34,prs)

# normalising to add up to 1
plot(divisionsDistro/sum(divisionsDistro)*100, 
     type = "l", col = "red", lwd = 3,
     main = "Distribution of cell Divisions at one time point",
     ylab = "Percentage", xlab = "No. of cumulated cell divisions")

# Generate  distribution over multiple chronological time points
generations <- seq(0,60,2)
distro <- sapply(generations, function(x) sapply(0:100, f,x*365,0.0037,10,prs))

# Normalize the distribution
normlzd <- t(distro)/apply(distro, 2, sum) *100

# save the data for the 3D plot in pythong
write.csv(normlzd, "3DplotData.csv")

normlzdf <- as.data.frame(t(normlzd))
colnames(normlzdf) <- 20+ generations
normlzdf$Divisions <- 0:100

# Create long data for plotting
plottingdf <- melt(normlzdf ,  id.vars = 'Divisions', variable.name = 'years')

# Provisional 2D Plot. 3D Plot generated in python (code in the repo)
BgkCol <- "white"

p <- ggplot(plottingdf, aes(Divisions,value)) + geom_line(aes(colour = years)) 
p <- p +  theme(
  axis.text=element_text(size=15),
  axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = adjustcolor(BgkCol, alpha.f = 0.1)),
  plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
) + labs(x = "No. of accumulated divisions", y = "Percentage") +
  ggtitle("Distribution of T cell replications  with age") 
p
