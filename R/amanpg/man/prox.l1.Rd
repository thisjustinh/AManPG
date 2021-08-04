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
Shixiang Chen, Shiqian Ma, Lingzhou Xue and Hui Zou. An Alternating Manifold Proximal Gradient Method for Sparse PCA and Sparse CCA. 2019.
}
\author{
Justin Huang, Benjamin Jochem, Shiqian Ma, and Lingzhou Xue
}
