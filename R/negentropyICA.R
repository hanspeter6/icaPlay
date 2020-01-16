## Negentropy-based Fixed Iteration algorithm to extract signals by deflation

# Using deflationary orthogonolisation as in Section 8.4.2 (p194)
# with iteration for one-unit algorithm based on equation 8.43 (p189)

negentropyICA <- function(x, g = g2, max.iters = 10, init.w = diag(m), tol = 1e-4) {
  
  # source additional functions
  source("my_Whiten.R")
  
  # expressions of gs
  g1 <- expression(tanh(a1 * y))
  g2 <- expression(y * exp(-y^2/2))
  g3 <- expression(y^3)
  
  # a function to normalise a vector
  norm_vec <- function(x) {sqrt(sum(x^2))}
  
  # a function based on negentropy
  # the minimising function: fast fixed point algorithm 8.43
  negent <- function(wi, z, g, a1 = 1) {
    
    m <- nrow(z) # identify number of componets/sources
    
    y <- wi %*% z # an estimate of a single source
    
    f <- function(y, g, a1) {eval( g[[1]] )} # identify the necessary g as a function from expression
    
    gy <- f(y, g, a1) # g(w'z) / g(y)
    
    dg <- D(g, 'y') # g'(y)
    
    # the expected values of the two terms in the equation:
    as.vector( (rowMeans( z * matrix(rep(gy,m), byrow = TRUE, nrow = m)) - 
                  rowMeans(matrix(rep(eval(dg), m), byrow = TRUE, nrow = m)) * t(wi) ) )
    
  }
  
  # initialsing variables
  m <- nrow(x)                         # number of components/sources/mixtures (for now the same thing)
  n <- ncol(x)                         # number of observations
  z <- my_Whiten(t(x))$z                  # whiten my mixtures
  w <- init.w/norm_vec(init.w)         # initial unmixing matrix
  iters <- max.iters                   # maximum iterations
  
  # initialising containers
  ws <- vector("list", m)
  thetas <- vector("list", m)
  is <- vector()
  
  # deflationary algorithm, extracting componentwise and maintaining orthogonolisation
  for(p in 1: m) {
    
    wp <- w[p,]
    theta_vector <- vector()
    
    # iterations to convergence for each component
    i <- 1
    repeat {
      
      y <- wp %*% z             # estimated signal vector
      
      wp <- negent(wp, z, g2, a1 = 1) # estimate update of unmixing vector
      
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
      val <- wp %*% w[p,]/(length(wp) * length(w[p,]))
      theta <- acos(round(val,6)) # needed to round to deal with NaN produced by floats when arg is -1
      theta_vector[i] <- theta
      
      # updating W and iteration container
      w[p,] <- wp
      is[p] <- i
      
      # checking convergence conditions based on "tol" = tolerance
      if(i == iters | abs(theta) < tol | abs(theta - pi) < tol) { break } # 
      
      # update iteration counter
      i <- i + 1
      
    }
    
    # collecting final W (unmixing matrix) for each component. Last item is final W of the algorithm
    ws[[p]] <- w
    # collecting thetas
    thetas[[p]] <- theta_vector
    
  }
  
  #identify a list of desired output objects
  out <- list(ws = ws,              # the W matrices after each component has been identified
              W = ws[[m]],          # the final estimated unmixing matrix W
              S = ws[[m]] %*% z,       # the signals in mxn matrix
              thetas = thetas,      # angles
              iters = is)           # iterations per deflation
  
}
