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

###### TRIALS

k <- kurtosisICA(x1) # kurtosis
par(mfrow = c(3,2), mar = c(1,1,1,1))

plot(density(s1[1,]), main = "Source 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "", cex = 0.1) 
plot(density(s1[2,]), main = "Source 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(x1[1,]), main = "Mix 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(x1[2,]), main = "Mix 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(k$S[1,]), main = "IC 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(k$S[2,]), main = "IC 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
dev.off()

ml <- maxLikeICA(x1) # maximum liklihood
par(mfrow = c(3,2), mar = c(1,1,1,1))

plot(density(s1[1,]), main = "Source 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "", cex = 0.1) 
plot(density(s1[2,]), main = "Source 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(x1[1,]), main = "Mix 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(x1[2,]), main = "Mix 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(ml$S[1,]), main = "IC 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(ml$S[2,]), main = "IC 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
dev.off()


neg <- negentropyICA(x1) # negentropy
par(mfrow = c(3,2), mar = c(1,1,1,1))

plot(density(s1[1,]), main = "Source 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "", cex = 0.1) 
plot(density(s1[2,]), main = "Source 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(x1[1,]), main = "Mix 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(x1[2,]), main = "Mix 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(ml$S[1,]), main = "IC 1", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
plot(density(ml$S[2,]), main = "IC 2", axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "") 
dev.off()


# Amari Error comparisons
my_Amari(solve(a1), k$W)
my_Amari(solve(a1), ml$W)
my_Amari(solve(a1), neg$W)
