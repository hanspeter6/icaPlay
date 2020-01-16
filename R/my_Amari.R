# define the amari error (see Amari etal 1996. A new learning...)
# I got it from Sorensen
# input: inverse mixing matrix A and unmixing matrix W
# (NB assumption is symmetrical and equal dimension matrices)

my_Amari <- function(A,W) {

  # initialising
  AInv <- A                   # inverse of A (should be entered correctly)
  WInv <- solve(W)            # inverse of W
  m <- nrow(A)                # dimensions
  AInvWInv <- AInv %*% WInv  # setting up bij matrix

  # define the two key terms in the expression
  terms1 <- rowSums(abs(AInvWInv)) / apply(abs(AInvWInv), 1, max) - 1
  terms2 <- colSums(abs(AInvWInv)) / apply(abs(AInvWInv), 2, max) - 1

  # final calculation for output
  return( 1/(2 * m) * ( sum(terms1) + sum(terms2) ) )

}
