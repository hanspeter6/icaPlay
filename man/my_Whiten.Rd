\name{my_Whiten}
\alias{my_Whiten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Helper Function to Whiten Datasets
%%  ~~function to do ... ~~
}
\description{ This function takes in a single dataset X with
$$X^2$$ and outputs blah

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_Whiten(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A real-valued matrix with observations.. and variables ...
%%     ~~Describe \code{x} here~~
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
\item{v }{The whitening matrix}
\item{z }{The dataset in vars and observations, whitened}

}
\references{
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

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
    x_centered <- apply(x, 2, function(x) {
        x - mean(x)
    })
    cvx <- cov(x_centered)
    eigens <- eigen(cvx)
    E <- eigens$vectors
    D_invsqrt <- diag(1/sqrt(eigens$values))
    myBasicWhiteningMatrix <- E \%*\% D_invsqrt \%*\% t(E)
    z <- myBasicWhiteningMatrix \%*\% t(x_centered)
    list(v = myBasicWhiteningMatrix, z = z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~whiten }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
