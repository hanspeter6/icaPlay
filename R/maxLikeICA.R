# Maximum-Likelikhood- based fixed point ICA (can only be symmetrical...) 
# implementation of the algorithm Table 9.2 p(212). Prewhitening version...
# also, for now only g = tanh()

maxLikeICA <- function(x, g = g, max.iters = 12, init.w = diag(m), tol = 1e-4) { # get back to vary g????
  
  # source additional functions
  source("my_Whiten.R")
  
  # expressions of gs
  g <- expression(tanh(y))
  
  # a function to normalise a vector
  norm_vec <- function(x) {sqrt(sum(x^2))}
  
  # a function to calc Expected Value of the dot product of two vectors # conformable of form E(xy') (m x n)(( n x m))
  # NBNB double check with Bosco...paper
  exp_xty <- function(x,y) { 
    
    n <- ncol(x)
    m <- nrow(x)
    
    mat <- matrix(rep(0, m * m), nrow = m)
    
    for(i in 1:n) {
      
      mat <- mat + matrix(x[,i], nrow = m) %*% matrix(t(y)[i,], nrow = 1)
      
    }
    
    return(mat/n)
    
  }
  
  # initialsing variables
  m <- nrow(x)                              # number of components/sources/mixtures (for now the same thing)
  n <- ncol(x)                              # number of observations
  z <- my_Whiten(t(x))$z                      # whiten the mixtures
  w <- init.w / apply(init.w, 1, norm_vec)  # initial unmixing matrix (nb do norming!!!)
  iters <- max.iters                        # maximum iterations, normalised
   
  dg <- D(g,'y') # derivative of g (tanh)

  # # initialising containers
  ws <- vector("list")
  ws[[1]] <- w
  distances <- vector()
  distances[1] <- 0
  
  # initialise iteration counter
  i <- 1
  
  repeat {
    
    # Step 3 in algorithm
    y <- w %*% z                                     # estimated signal vector
    beta_i <- rowMeans( y * eval(g) )
    alpha_i <- -1 / ( beta_i + rowMeans(eval(dg)) )
    
    # Step 4
    w <- w + diag(alpha_i) %*% (diag(beta_i) + exp_xty( eval(g), y )) %*% w
    
    # Step 5 orthogonolise and normalise

    ## calculate inverse sq root of WW'
    # determine eigen vectors and values
    eigens <- eigen(w %*% t(w))
    
    # determine components of Eigen Decomposition
    E <- eigens$vectors
    D_invsqrt <- diag(1/sqrt(eigens$values))
    
    # the inverse sqrt of WW'
    ww_inSqrt <- E %*% D_invsqrt %*% t(E)
    
    # and update w
    w <- ww_inSqrt %*% w

    # collect iteration's w
    ws[[i + 1]] <- w
    
    # # absolute distance between previous and current w and collect them
    distances[i + 1] <- abs( norm(ws[[i]] - ws[[i+1]]) )
    
    # checking convergence conditions based on "tol" = tolerance
    if(i == iters | abs(distances[i] - distances[i+1]) < tol) { break } #
    
    # update iteration counter
    i <- i + 1
  
  }
  
  # #identify a list of desired output objects
  out <- list(ws = ws,
              W = w,
              S = w %*% z,
              iters = i,
              distances = distances)
}
