\name{spca.amanpg}
\alias{spca.amanpg}
\title{
Alternating Manifold Proximal Gradient algorithm for Sparse PCA
}
\description{
Performs sparse principal component analysis on the input matrix using an alternating manifold
proximal gradient (AManPG) method
}
\usage{spca.amanpg(z, lambda1, lambda2, f_palm = 1e5, x0 = NULL, y0 = NULL, k = 0, type = 0,
       gamma = 0.5, maxiter = 1e4, tol = 1e-5, normalize = TRUE, verbose = FALSE)
}
\arguments{
  \item{z}{Either the data matrix or sample covariance matrix}
  \item{lambda1}{List of parameters of length n for L1-norm penalty}
  \item{lambda2}{L2-norm penalty term}
  \item{f_palm}{Upper bound for the gradient value to reach convergence, default value is 1e5}
  \item{x0}{Initial x-values for the gradient method, default value is the first n right singular vectors}
  \item{y0}{Initial y-values for the gradient method, default value is the first n right singular vectors}
  \item{k}{Number of principal components desired, default is 0 (returns min(n-1, p) principal components)}
  \item{type}{If 0, b is expected to be a data matrix, and otherwise b is expected to be a covariance matrix; default is 0}
  \item{gamma}{Parameter to control how quickly the step size changes in each iteration, default is 0.5}
  \item{maxiter}{Maximum number of iterations allowed in the gradient method, default is 1e4}
  \item{tol}{Tolerance value required to indicate convergence (calculated as difference between iteration f-values), default is 1e-5}
  \item{normalize}{Center and normalize rows to Euclidean length 1 if True, default is True}
  \item{verbose}{Function prints progress between iterations if True, default is False}

}
\value{
  \item{iter}{total number of iterations executed in the algorithm}
  \item{f_amanpg}{final gradient value}
  \item{sparsity}{Number of sparse loadings (loadings == 0) divided by number of all loadings}
  \item{time}{execution time in seconds}
  \item{x}{corresponding matrix in subproblem to the loadings}
  \item{loadings}{loadings of the sparse principal components}
}
\references{
Chen, S., Ma, S., Xue, L., and Zou, H. (2020) "An Alternating Manifold Proximal Gradient Method for Sparse Principal Component Analysis and Sparse Canonical Correlation Analysis" *INFORMS Journal on Optimization* 2:3, 192-208
}
\author{
Shixiang Chen, Justin Huang, Benjamin Jochem, Shiqian Ma, Lingzhou Xue and Hui Zou
}

\examples{
#see SPCA.R for a more in-depth example
      d <- 500  # dimension
      m <- 1000 # sample size
      a <- normalize(matrix(rnorm(m * d), m, d))
      lambda1 <- 0.1 * matrix(data=1, nrow=4, ncol=1)
      x0 <- svd(a, nv=4)$v
      sprout <- spca.amanpg(a, lambda1, lambda2=Inf, f_palm=1e5, x0=x0, y0=x0, k=4, type=0, 
                            gamma=0.5, maxiter=1e4, tol=1e-5, normalize = FALSE, verbose=FALSE)
      print(paste(sprout$iter, "iterations,", sprout$sparsity, "sparsity,", sprout$time))

      #extract loadings
      #print(sprout$loadings)

}
