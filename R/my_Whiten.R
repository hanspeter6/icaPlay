# basic ZCA whitening function: X, datamatrix with rows observations and columns variables
# assume square and positive definite covariance matrix

my_Whiten <- function(x) {

  #center the variables
  x_centered <- apply(x, 2, function(x) {x - mean(x)})

  # determine the covariance matrix
  cvx <- cov(x_centered)

  # determine eigen vectors and values
  eigens <- eigen(cvx)

  # determine the components of Eigen Decomposition
  E <- eigens$vectors
  D_invsqrt <- diag(1/sqrt(eigens$values))

  # determine the whitening matrix
  myBasicWhiteningMatrix <- E %*% D_invsqrt %*% t(E) # compares with ZCA method in package

  # the whitened dataset
  z <- myBasicWhiteningMatrix %*% t(x_centered)

  # attach the whitening matrix
  list(v = myBasicWhiteningMatrix, z = z)

}

