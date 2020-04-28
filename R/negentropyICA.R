## Negentropy-based Fixed Iteration algorithm to extract signals by deflation

# Using deflationary orthogonolisation as in hyvarinen's book Section 8.4.2 (p194)
# with iteration for one-unit algorithm based on the same book's equation 8.43 (p189)

negentropyICA <- function(x, g = g2, max.iters = 10, init.w = diag(m), tol = 1e-4) {
  
  ## additional functions
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
  
  # a function based on negentropy
  # the minimising function: fast fixed point algorithm 8.43 in hyvarinen's book
  negent <- function(wi, z, g, a1 = 1) {
    
    m <- nrow(z) # identify number of componets/sources
    
    y <- wi %*% z # an estimate of a single source
    
    f <- function(y, g, a1) {eval( g[[1]] )} # identify the necessary g as a function from expression
    
    gy <- f(y, g, a1) # g(w'z) or g(y)
    
    dg <- D(g, 'y') # g'(y)
    
    # the expected values of the two terms in the equation:
    as.vector( (rowMeans( z * matrix(rep(gy,m), byrow = TRUE, nrow = m)) - 
                  rowMeans(matrix(rep(eval(dg), m), byrow = TRUE, nrow = m)) * t(wi) ) )
    
  }
  
  # expressions of gs
  g1 <- expression(tanh(a1 * y))
  g2 <- expression(y * exp(-y^2/2))
  g3 <- expression(y^3)
  
  # initialsing variables
  m <- nrow(x)                         # number of components/sources/mixtures (for now the same thing)
  n <- ncol(x)                         # number of observations
  z <- whiten(t(x))                    # whiten my mixtures
  w <- init.w/norm_vec(init.w)         # initial unmixing matrix normalised
  iters <- max.iters                   # maximum iterations
  
  # initialising containers
  ws <- vector("list", m)
  thetas <- vector("list", m)
  is <- vector()
  
  # deflationary algorithm, extracting componentwise and maintaining orthogonolisation
  for(p in 1: m) {
    
    wp <- w[p,]
    theta_vector <- vector()
    ws_vector <- vector("list")
    
    # iterations to convergence for each component
    i <- 1
    repeat {
      
      y <- wp %*% z                    # estimated signal vector
      
      wp <- negent(wp, z, g2, a1 = 1)  # estimate update of unmixing vector
      
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
      val <- sum(wp*w[p,]) / ( sqrt(sum(wp * wp)) * sqrt(sum(w[p,] * w[p,])) )  # expression inside acos
      theta <- acos(round(val,6))                # needed to round to deal with NaN produced by floats when arg is -1
      theta_vector[i] <- theta                   # collecting thetas by 
      ws_vector[[i]] <- wp
      
      # updating W and iteration container
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
    
  }
  
  #identify a list of desired output objects
  out <- list(ws = ws,              # the W matrices after each component has been identified
              W = w,          # the final estimated unmixing matrix W
              S = w %*% z,       # the signals in mxn matrix
              thetas = thetas,      # angles
              iters = is)           # iterations per deflation
  
}


