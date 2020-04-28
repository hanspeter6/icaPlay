## Anderson-Darling-based (Fixed Iteration) algorithm to extract signals by deflation:
# Source: Razali et al, 2011
# Using deflationary orthogonolisation as in Section 8.4.2 (p194)

adICA <- function(x, max.iters = 10, init.w = diag(m), tol = 1e-4) {
  
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
  
  # My experimental algorithm
  ad_grad <- function(wi, z) {
    
    require(dplyr)
    
    m <- nrow(z)
    
    y <- as.vector(wi %*% z)
    
    yden <- dnorm(y) # stanard normal density of y
    ycum <- pnorm(y) # standard normal distribution of y
    alp <- yden/ycum # 
    bet <- yden/(1 - ycum) #
    lhs <- matrix(rep(alp, m), byrow = TRUE, nrow = m) * z
    rhs <- matrix(rep(bet, m), byrow = TRUE, nrow = m) * z
    fin <- lhs - t(apply(rhs, 1, rev))
    rowMeans(fin)
    # 
    # 
    # 
    # 
    # mat <- as_tibble(t(rbind(y = y, z1 = z[1,], z2 = z[2,], z3 = z[3,] )))
    # 
    # mat_fin <- mat %>%
    #   arrange(y)  %>%
    #   mutate(y_rev = rev(y)) %>%
    #   mutate(z1_rev = rev(z1)) %>%
    #   mutate(z2_rev = rev(z2)) %>%
    #   mutate(z3_rev = rev(z3)) %>%
    #   mutate(den_y = dnorm(y)) %>%
    #   mutate(dist_y = pnorm(y)) %>%
    #   mutate(alpha = den_y/dist_y) %>%
    #   mutate(den_y_rev = dnorm(y_rev)) %>%
    #   mutate(dist_y_rev = pnorm(y_rev)) %>%
    #   mutate(beta = den_y_rev/(1-dist_y_rev)) %>%
    #   mutate(lhs1 = alpha * z1) %>%
    #   mutate(lhs2 = alpha * z2) %>%
    #   mutate(lhs3 = alpha * z3) %>%
    #   mutate(rhs1 = beta * z1_rev) %>%
    #   mutate(rhs2 = beta * z2_rev) %>%
    #   mutate(rhs3 = beta * z3_rev) %>%
    #   mutate(fin1 = lhs1 - rhs1) %>%
    #   mutate(fin2 = lhs2 - rhs2) %>%
    #   mutate(fin3 = lhs3 - rhs3)
    # 
    # with(mat_fin, c(mean(fin1), mean(fin2), mean(fin3)))
    
  }
  
  # initialsing variables
  m <- nrow(x)                         # number of components/sources/mixtures (for now the same thing)
  n <- ncol(x)                         # number of observations
  z <- whiten(t(x))                    # whiten my mixtures
  w <- init.w/norm_vec(init.w)         # initial unmixing matrix normalised
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
      
      wp <- wp + 0.5 * ad_grad(wp, z)    # update with the estimated gradient with learning rate alpha = 1 for now
      
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
      ks_vector[i] <- k                          # note. Not abs(k)
      ws_vector[[i]] <- wp                       # collecing each iteration w
      
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
    
    # collecting kurtosis
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
