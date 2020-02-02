## Here I want to set up some random datasets to play with

library(rmutil)

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
# first id which mix, 
a <- a2
s <- s2
x <- x2
m <- m2

k <- kurtosisICA(x) # kurtosis

par(mfrow = c(3,m), mar = c(1,1,1,1))
layout(matrix(c(1:(3*m)), ncol = m, byrow = FALSE))
for(i in 1:m) {
  plot(density(s[i,]), main = paste("Source", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
  plot(density(x[i,]), main = paste("Mix", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
  plot(density(k$S[i,]), main = paste("IC", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
  
}
dev.off()

ml <- maxLikeICA(x) # maximum liklihood
par(mfrow = c(3,m), mar = c(1,1,1,1))
layout(matrix(c(1:(3*m)), ncol = m, byrow = FALSE))
for(i in 1:m) {
  plot(density(s[i,]), main = paste("Source", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
  plot(density(x[i,]), main = paste("Mix", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
  plot(density(ml$S[i,]), main = paste("IC", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
  
}
dev.off()

neg <- negentropyICA(x) # negentropy
par(mfrow = c(3,m), mar = c(1,1,1,1))
layout(matrix(c(1:(3*m)), ncol = m, byrow = FALSE))
for(i in 1:m) {
  plot(density(s[i,]), main = paste("Source", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
  plot(density(x[i,]), main = paste("Mix", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
  plot(density(neg$S[i,]), main = paste("IC", i), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
  
}
dev.off()

# Amari Error comparisons
my_Amari(solve(a), k$W)
my_Amari(solve(a), ml$W)
my_Amari(solve(a), neg$W)

# more on stuff here....