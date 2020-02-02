\name{my_Amari}
\alias{my_Amari}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Amari Error
%%  ~~function to do ... ~~
}
\description{Measures a distance between a mixing matrix inverse of A and an estimate of an unmixing matrix W.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_Amari(A, W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ The inverse of the mixing matrix
%%     ~~Describe \code{A} here~~
}
  \item{W}{ The estimated unmixing matrix
%%     ~~Describe \code{W} here~~
}
}
\details{ This is an implementation of the algorithm found in \cite{Sorensen, J. 2008. An Assessment of Hermite Function-based Approximations of Mutual Information Applied to Independent Component Analysis. Entropy. Vol. 10. pp745-756}

\cite{Sorensen, 2008} cites \cite{Amari et al. 1996} for the expression.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ The Amari Error/Distance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ See details
%% ~put references to the literature/web site here ~
}
\author{ Hans-Peter Bakker
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
set.seed(1)
A <- matrix(runif(4), nrow = 2)
set.seed(2)
W <- matrix(runif(4), nrow = 2)

my_Amari(solve(A), W)
[1] 0.6755254


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (A, W) 
{
    AInv <- A
    WInv <- solve(W)
    m <- nrow(A)
    AInvWInv <- AInv \%*\% WInv
    terms1 <- rowSums(abs(AInvWInv))/apply(abs(AInvWInv), 1, 
        max) - 1
    terms2 <- colSums(abs(AInvWInv))/apply(abs(AInvWInv), 2, 
        max) - 1
    return(1/(2 * m) * (sum(terms1) + sum(terms2)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
