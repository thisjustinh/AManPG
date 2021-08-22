\name{normalize}
\alias{normalize}
\title{
Matrix Normalization
}
\description{
Center the input matrix to mean 0 and scale to Euclidean length 1
}
\usage{
normalize(x, center=TRUE, scale=TRUE)
}
\arguments{
  \item{x}{matrix to be normalized}
  \item{center}{centers the input matrix to mean 0 if TRUE, default if TRUE}
  \item{scale}{scales the input matrix to Euclidean length 1 if TRUE, default is TRUE}
}
\value{
\item{x}{normalized matrix}
}
\author{
Shixiang Chen, Justin Huang, Benjamin Jochem, Shiqian Ma, Lingzhou Xue and Hui Zou
}
