\name{my_Amari}
\alias{my_Amari}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to measure the Amari error
%%  ~~function to do ... ~~
}
\description{ This function measures the Amari Error. A measure of the distance between an ideal unmixing matrix (the inverse of the mixing matrix) and an estimated unmixing matrix.
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
\details{

Amari's Error provides a measure of the distance between the inverse of the mixing matrix A and the estimated unmixing matrix W.

%This function implements the description found in
%\cite{Sorensen, J. 2008. An Assessment of Hermite Function Based Approximations of Mutual Information Applied to Independent Component Analysis. Entropy Volume 10 pp 745-756}

%This description references \cite{Amari S.,Chichochi A., Yank H.H. 1996. A new learning algorithm for blind signal separation. Advances in Neural Information Processing Systems. Volume 8, pp 757-763}

}
%%  ~~ If necessary, more details than the description above ~~
}
\value{

The function outputs only the computed Amari's Error.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ See Details
%% ~put references to the literature/web site here ~
}
\author{

Hans-Peter Bakker
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
set.seed(1)
W <- matrix(runif(4), nrow = 2)

AInv <- solve(A)

my_Amari(A,Anv)
#[1] 0.5349698

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
