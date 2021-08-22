\name{prox.l1}
\alias{prox.l1}
\title{
Proximal L1 Mapping
}
\description{
Calculates the proximal L1 mapping for the given input matrix
}
\usage{
prox.l1(z, lambda, r)
}
\arguments{
  \item{z}{input matrix}
  \item{lambda}{parameters for calculating proximal L1 mapping}
  \item{r}{number of columns used in matrix}
}

\value{
  \item{x_prox}{proximal L1 Mapping}
}
\references{
Chen, S., Ma, S., Xue, L., and Zou, H. (2020) "An Alternating Manifold Proximal Gradient Method for Sparse Principal Component Analysis and Sparse Canonical Correlation Analysis" *INFORMS Journal on Optimization* 2:3, 192-208
}
\author{
Shixiang Chen, Justin Huang, Benjamin Jochem, Shiqian Ma, Lingzhou Xue and Hui Zou
}
