## Here I want to set up some random datasets to play with

library(rmutil)

# function for plotting
pl_den <- function(obj, signals = s, mixes = x) {
  
  m <- nrow(signals)
  par(mfrow = c(3,m), mar = c(1,1,1,1))
  layout(matrix(c(1:(3*m)), ncol = m, byrow = FALSE))
  
  for(i in 1:m) {
    
    plot(density(signals[i,]), main = paste("Source", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
    plot(density(mixes[i,]), main = paste("Mix", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
    plot(density(obj$S[i,]), main = paste("IC", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
    
  }
}

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

###### TRIALS
# first id which mix 
a <- a3
s <- s3
x <- x3
m <- m3

#
k <- kurtosisICA(x) # kurtosis
ml <- maxLikeICA(x) # maximum liklihood
neg <- negentropyICA(x) # negentropy

# plotting
pl_den(ml, signals = s, mixes = x)
pl_den(k, signals = s, mixes = x)
pl_den(neg, signals = s, mixes = x)
dev.off()

# Amari Error comparisons
my_Amari(solve(a), k$W)
my_Amari(solve(a), ml$W)
my_Amari(solve(a), neg$W)

# more on stuff here....