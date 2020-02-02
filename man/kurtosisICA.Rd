\name{kurtosisICA}
\alias{kurtosisICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fixed-point algorithm using kurtosis
%%  ~~function to do ... ~~
}
\description{A deflationary fixed-point algorithm to estimate an unmixing matrix using kurtosis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{kurtosisICA(x, max.iters = 10, init.w = diag(m), tol = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A mixture set with variables in rows and observations in columns.
%%     ~~Describe \code{x} here~~
}
  \item{max.iters}{The maximum number of iterations. Default = 10.
%%     ~~Describe \code{max.iters} here~~
}
  \item{init.w}{The initial unmixing matrix. Default is identity.
%%     ~~Describe \code{init.w} here~~
}
  \item{tol}{The tolerance for theta (the angle between subsequent estimates of the relevant W vector). Default is 1e-4.
%%     ~~Describe \code{tol} here~~
}
}
\details{This fixed-point algorithm using the absolute value of kurtosis is an implementation of the deflationary algorithm described in section 8.3.5 in \cite{Hyvarinen, Karunen and Oja. 2001. Independent Component Analysis}.

More specifically this is the implementation of iteration equation 8.20 on page 178 and the algorithm indicated in section 8.4.2 page 194.

The unmixing matrix is normalised throughout and the data is whitened in a pre-processing step.}
%%  ~~ If necessary, more details than the description above ~~
\value{
\item{ws}{A list per component of unmixing matrices by iteration}
\item{W}{The final converged estimated unmixing matrix}
\item{S}{The estimated independent components by row with observations by column}
\item{thetas}{A list per component of angles in radians between subsequent estimates of W}
\item{ks}{A list per component of kurtosis per iteration}
\item{iters}{The number of iterations to convergence for each component}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{See details.
%% ~put references to the literature/web site here ~
}
\author{Hans-Peter Bakker
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
library(rmutil)
s <- rbind(runif(1000), rlaplace(1000))
a <- matrix(runif(4), nrow = 2)
x <- a %*% s

out <- kurtosisICA(x)

par(mfrow = c(3,2))
plot(density(s[1,]), main = "signal 1")
plot(density(s[2,]), main = "signal 2")
plot(density(x[1,]), main = "mixture 1")
plot(density(x[2,]), main = "mixture 2")
plot(density(out$S[1,]), main = "Indep. Comp. 1")
plot(density(out$S[2,]), main = "Indep. Comp. 2")

## The function is currently defined as
function (x, max.iters = 10, init.w = diag(m), tol = 1e-04) 
{
    source("my_Whiten.R")
    norm_vec <- function(x) {
        sqrt(sum(x^2))
    }
    kurt_grad <- function(wi, z) {
        m <- nrow(z)
        y <- wi \%*\% z
        y3 <- y^3
        yy3 <- matrix(rep(y3, m), nrow = m, byrow = TRUE)
        rowMeans(z * yy3)
    }
    m <- nrow(x)
    n <- ncol(x)
    z <- my_Whiten(t(x))$z
    w <- init.w/norm_vec(init.w)
    iters <- max.iters
    ws <- vector("list", m)
    thetas <- vector("list", m)
    is <- vector()
    for (p in 1:m) {
        wp <- w[p, ]
        theta_vector <- vector()
        i <- 1
        repeat {
            y <- wp \%*\% z
            k <- mean(y^4) - 3
            g <- kurt_grad(wp, z)
            wp <- g - 3 * wp
            if (p > 1) {
                terms <- matrix(rep(0, (p - 1) * m), ncol = m)
                for (j in 1:(p - 1)) {
                  terms[j, ] <- t(wp \%*\% w[j, ]) \%*\% w[j, ]
                }
                wp <- wp - colSums(terms)
            }
            wp <- wp/norm_vec(wp)
            val <- wp \%*\% w[p, ]/(length(wp) * length(w[p, ]))
            theta <- acos(round(val, 6))
            theta_vector[i] <- theta
            w[p, ] <- wp
            is[p] <- i
            if (i == iters | abs(theta) < tol | abs(theta - pi) < 
                tol) {
                break
            }
            i <- i + 1
        }
        ws[[p]] <- w
        thetas[[p]] <- theta_vector
    }
    out <- list(ws = ws, W = ws[[m]], S = ws[[m]] \%*\% z, thetas = thetas, 
        iters = is)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ica }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kurtosis }% __ONLY ONE__ keyword per line
