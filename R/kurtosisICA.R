## Kurtosis-based Fixed Iteration algorithm to extract signals by deflation:

# Using deflationary orthogonolisation as in Section 8.4.2 (p194)
# with iteration for one-unit algorithm based on equation 8.20 (p178)

kurtosisICA <- function(x, max.iters = 10, init.w = diag(m), tol = 1e-4) {
  
  # source additional functions
  # basic ZCA whitening function: X, datamatrix with rows observations and columns variables
  # assume square and positive definite covariance matrix
  
  whiten <- function(x) {
    
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
    myBasicWhiteningMatrix %*% t(x_centered)
    
  }
  
  # a function to normalise a vector
  norm_vec <- function(x) {sqrt(sum(x^2))}
  
  # Estimating the gradient of absolute value of the kurtosis p177 equation 8 .14
  kurt_grad <- function(wi, z) {
    
    m <- nrow(z)
    y <- wi %*% z
    y3 <- y^3
    yy3 <- matrix(rep(y3,m), nrow = m, byrow = TRUE)
    rowMeans(z * yy3) - (3 * wi)
    
  }

  # initialsing variables
  m <- nrow(x)                         # number of components/sources/mixtures (for now the same thing)
  n <- ncol(x)                         # number of observations
  z <- whiten(t(x))                 # whiten my mixtures
  w <- init.w/norm_vec(init.w)         # initial unmixing matrix
  iters <- max.iters                   # maximum iterations
  
  # initialising containers
  ws <- vector("list", m)
  thetas <- vector("list", m)
  ks <- vector("list", m)
  is <- vector()
  
  # deflationary algorithm, extracting componentwise and maintaining orthogonolisation
  for(p in 1: m) {
    
    wp <- w[p,]
    theta_vector <- vector()
    ws_vector <- vector("list")
    ks_vector <- vector()
    
    # iterations to convergence for each component
    i <- 1
    repeat {
      
      y <- wp %*% z             # estimated signal vector
      
      k <- mean(y^4) - 3        # measuring kurtosis
      
      wp <- kurt_grad(wp, z)    # update with the estimated gradient
      
      # GS Orthogonalization
      if(p > 1) {
        
        terms <- matrix(rep(0, (p-1) * m), ncol = m)
        
        for(j in 1:(p-1)) {
          
          terms[j,] <- t(wp %*% w[j,]) %*% w[j,]
        }
        
        wp <- wp - colSums(terms)
        
      }
      
      # normalizing
      wp <- wp/norm_vec(wp)
      
      # angle between previous and current wp and collect them in a vector
      val <- sum(wp*w[p,]) / ( sqrt(sum(wp * wp)) * sqrt(sum(w[p,] * w[p,])) )
      theta <- acos(round(val,6)) # needed to round to deal with NaN produced by floats when arg is -1
      theta_vector[i] <- theta   # collecting thetas by 
      ks_vector[i] <- k           # note. Not abs()
      ws_vector[[i]] <- wp

      # updating W and other containers
      w[p,] <- wp
      is[p] <- i
      
      # checking convergence conditions based on "tol" = tolerance
      if(i == iters | abs(theta) < tol | abs(theta - pi) < tol) { break } # 
      
      # update iteration counter
      i <- i + 1
      
    }
    
    # collecting final W (unmixing matrix) for each component. Last item is final W of the algorithm
    ws[[p]] <- ws_vector
    # collecting thetas
    thetas[[p]] <- theta_vector
    ks[[p]] <- ks_vector
    
  }
  
  #identify a list of desired output objects
  out <- list(ws = ws,              # the W matrices after each component has been identified
              W = w,          # the final estimated unmixing matrix W
              S = w %*% z,    # the signals in mxn matrix
              thetas = thetas,      # angles
              ks = ks,             # kurtosis
              iters = is)           # iterations per deflation
  
}


