\name{maxLikeICA}
\alias{maxLikeICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Likelihood Independent Component Analysis
%%  ~~function to do ... ~~
}
\description{Symmetrical Estimate of Independent Components Using Maximum Likelihood.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maxLikeICA(x, g = g, max.iters = 12, init.w = diag(m), tol = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A mixture set with variables in rows and observations in columns.
%%     ~~Describe \code{x} here~~
}
  \item{g}{Current version only uses tanh()
%%     ~~Describe \code{g} here~~
}
  \item{max.iters}{The maximum number of iterations. Default = 10.
%%     ~~Describe \code{max.iters} here~~
}
  \item{init.w}{The initial unmixing matrix. Default is identity.
%%     ~~Describe \code{init.w} here~~
}
  \item{tol}{The tolerance for iteration differences...(?? come back here..). Default is 1e-4.
%%     ~~Describe \code{tol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, g = g, max.iters = 12, init.w = diag(m), tol = 1e-04) 
{
    source("my_Whiten.R")
    g <- expression(tanh(y))
    norm_vec <- function(x) {
        sqrt(sum(x^2))
    }
    exp_xty <- function(x, y) {
        n <- ncol(x)
        m <- nrow(x)
        mat <- matrix(rep(0, m * m), nrow = m)
        for (i in 1:n) {
            mat <- mat + matrix(x[, i], nrow = m) \%*\% matrix(t(y)[i, 
                ], nrow = 1)
        }
        return(mat/n)
    }
    m <- nrow(x)
    n <- ncol(x)
    z <- my_Whiten(t(x))$z
    w <- init.w/apply(init.w, 1, norm_vec)
    iters <- max.iters
    dg <- D(g, "y")
    ws <- vector("list")
    ws[[1]] <- w
    distances <- vector()
    distances[1] <- 0
    i <- 1
    repeat {
        y <- w \%*\% z
        beta_i <- rowMeans(y * eval(g))
        alpha_i <- -1/(beta_i + rowMeans(eval(dg)))
        w <- w + diag(alpha_i) \%*\% (diag(beta_i) + exp_xty(eval(g), 
            y)) \%*\% w
        eigens <- eigen(w \%*\% t(w))
        E <- eigens$vectors
        D_invsqrt <- diag(1/sqrt(eigens$values))
        ww_inSqrt <- E \%*\% D_invsqrt \%*\% t(E)
        w <- ww_inSqrt \%*\% w
        ws[[i + 1]] <- w
        distances[i + 1] <- abs(norm(ws[[i]] - ws[[i + 1]]))
        if (i == iters | abs(distances[i] - distances[i + 1]) < 
            tol) {
            break
        }
        i <- i + 1
    }
    out <- list(ws = ws, W = w, S = w \%*\% z, iters = i, distances = distances)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
