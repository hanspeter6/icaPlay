## Here I want to set up some random datasets to play with

library(rmutil) # simulate laplacian
library(tuneR) # simualte wave forms

source("R/kurtosisICA.R")
source("R/maxLikeICA.R")
source("R/negentropyICA.R")
source("R/my_Amari.R")
source("R/adICA.R")

# function for plotting matrix of densities (default) or time series (ts)
plots <- function(obj, signals = s, mixes = x, type = c( density, ts)) {
  
  m <- nrow(signals)
  par(mfrow = c(3,m), mar = c(1,1,1,1))
  layout(matrix(c(1:(3*m)), ncol = m, byrow = FALSE))
  
  for(i in 1:m) {
    
    plot(type(signals[i,]), main = paste("Source", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
    plot(type(mixes[i,]), main = paste("Mix", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
    plot(type(obj$S[i,]), main = paste("IC", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
    
  }
}
# # function for plotting time series
# pl_ts <- function(obj, signals = s, mixes = x) {
#   
#   m <- nrow(signals)
#   par(mfrow = c(3,m), mar = c(1,1,1,1))
#   layout(matrix(c(1:(3*m)), ncol = m, byrow = FALSE))
#   
#   for(i in 1:m) {
#     
#     plot(ts(signals[i,]), main = paste("Source", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
#     plot(ts(mixes[i,]), main = paste("Mix", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
#     plot(ts(obj$S[i,]), main = paste("IC", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
#     
#   }
# }

## mix 1:
m1 <- 2
n1 <- 1000
s1 <- rbind(runif(n1), rlaplace(n1))
a1 <- matrix(runif(m1^2), byrow = TRUE, nrow = m1)
x1 <- a1 %*% s1

## mix 2:
m2 <- 3
n2 <- 1000
s2 <- rbind(runif(n2), rlaplace(n2), rnorm(n2))
a2 <- matrix(runif(m2^2), byrow = TRUE, nrow = m2)
x2 <- a2 %*% s2

## mix 3:
m3 <- 2
n3 <- 1000
s3 <- rbind(c(rlaplace(500, m = -2), runif(500)), c(runif(1000)))
a3 <- matrix(runif(m3^2), byrow = TRUE, nrow = m3)
x3 <- a3 %*% s3

## mix 4:
m4 <- 4
n4 <- 1000
s4 <- rbind(sine(440, 1000)@left,
            noise(duration = 1000)@left,
            pulse(220, 1000)@left,
            sine(550, 1000)@left * pulse(600, 1000)@left)
a4 <- matrix(runif(m4^2), byrow = TRUE, nrow = m4)
x4 <- a4 %*% s4 

###### TRIALS
# first id which mix (1-4)
a <- a2
s <- s2
x <- x2
m <- m2

#
system.time(k <- kurtosisICA(x)) # kurtosis
system.time(ad <- adICA(x, max.iters = 100)) # anderson-darling
system.time(ml <- maxLikeICA(x)) # maximum liklihood
system.time(neg <- negentropyICA(x)) # negentropy

# plotting for densities
plots(ml, signals = s, mixes = x, type = density)
plots(k, signals = s, mixes = x, type = density)
plots(ad, signals = s, mixes = x, type = density)
plots(neg, signals = s, mixes = x, type = density)
dev.off()

#plotting for timeseries
plots(ml, signals = s, mixes = x, type = ts)
plots(k, signals = s, mixes = x, type = ts)
plots(ad, signals = s, mixes = x, type = ts)
plots(neg, signals = s, mixes = x, type = ts)
dev.off()


# Amari Error comparisons
my_Amari(solve(a), k$W)
my_Amari(solve(a), ad$W)
my_Amari(solve(a), ml$W)
my_Amari(solve(a), neg$W)

# more on stuff here....like signal to noise ...

# plots to consider convergences...

ad$iters
k$iters
