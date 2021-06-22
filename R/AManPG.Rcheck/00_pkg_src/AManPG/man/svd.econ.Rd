\name{svd.econ}
\alias{svd.econ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Economy-Size Singular Value Decomposition}
\description{Produces an economy-size decomposition of \eqn{m}-by-\eqn{n} matrix \eqn{x}, similar to the functionality in MATLAB.}
\usage{
svd.econ(x, type=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Matrix to be decomposed into \eqn{UDV}}
  \item{type}{If 0, will only consider the \code{m > n} case. If 1, will consider \code{m > n} and \code{m < n}}
}
\details{
If \code{m > n}, only the first \eqn{n} columns of \eqn{U} are computed, and \eqn{D} is of length \eqn{n}.

If \code{m < n}, only the first \eqn{m} columns of \eqn{V} are computed, and \eqn{D} is of length \eqn{m}.

If \code{m == n}, the base \code{svd} function will be executed.

Note that for \code{type = 0}, only the condition \code{m > n} will perform an economy-size decomposition. For \code{m <= n}, the base \code{svd} function is run.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
https://www.mathworks.com/help/matlab/ref/double.svd.html
}
\author{Justin Huang}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{svd}}
}
\examples{
set.seed(10)
x <- matrix(rnorm(3 * 4), 3, 4)
udv <- svd.econ(x)

udv.0 <- svd.econ(x, 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
