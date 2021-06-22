\name{spca_amanpg}
\alias{spca_amanpg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Alternativing Manifold Proximal Gradient Method for Sparse PCA}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Executes sparse PCA using the AManPG algorithm.
}
\usage{
spca_amanpg(B, mu, lambda, n, type, maxiter, tol, X0, Y0, F_palm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{Data or covariance matrix to be analyzed}
  \item{mu}{Good question!}
  \item{n}{Number of columns in x, y}
  \item{x0}{Initial x value for gradient algorithm}
  \item{y0}{Initial y value for gradient algorithm}
  \item{type}{
    0 indicates \code{b} is data matrix, 1 indicates \code{b} is covariance matrix, default is 0
  }
  \item{maxiter}{Maximum number of iterations prior to termination of algorithm, default is 1e4}
  \item{tol}{Tolerance level for detecting convergence, default is 1e-5}
  \item{f_palm}{Good question! I can see it's a convergence condition tho}
  \item{verbose}{Logical where TRUE displays progress during execution, default is FALSE}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{iter}{Number of iterations used in the algorithm}
  \item{f_amanpg}{Final value returned}
  \item{sparsity}{Ratio of all values in returned matrix with zeroes}
  \item{time}{Seconds used during execution}
  \item{x}{Final x matrix}
  \item{y_man}{Final y matrix of principal components}
}
\references{
Shixiang Chen, Shiqian Ma, Lingzhou Xue and Hui Zou. An Alternating Manifold Proximal Gradient Method for Sparse PCA and Sparse CCA. 2019.
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
function (x)
{
  }
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
