# Alternating Manifold Proximal Gradient Method (A-ManPG)

![pypi version](https://img.shields.io/pypi/v/sparsepca.svg) ![python version](https://img.shields.io/pypi/pyversions/sparsepca.svg) ![pypi downloads](https://img.shields.io/pypi/dm/sparsepca)

![CRAN version](https://www.r-pkg.org/badges/version/amanpg) ![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/amanpg)

- [Introduction](#introduction)
- [Installation](#installation)
- [Documentation](#documentation)
	- [Usage](#usage)
	- [Arguments](#arguments)
	- [Values](#values)
- [Quick Start](#quick-start)
	- [Python Quick Start](#python-quick-start)
	- [R Quick Start](#r-quick-start)
- [References](#references)

## Introduction

`sparsepca` and `amanpg` find sparse loadings in principal component analysis (PCA) via an alternating manifold proximal gradient method (A-ManPG). Seeking a sparse basis allows the leading principal components to be easier to interpret when modeling with high-dimensional data. Due to the nonsmoothness and nonconvexity numerical difficulties, A-ManPG is implemented to guarantee convergence. 

The package provides a function for performing sparse PCA and a function for normalizing data.

The authors of A-ManPG are Shixiang Chen, Shiqian Ma, Lingzhou Xue, and Hui Zou. The Python and R packages are maintained by Justin Huang and Benjamin Jochem. A MATLAB implementation is maintained by Shixiang Chen.

## Installation

To install the Python package, use `pip` to obtain `sparsepca` from PyPI:

```python
pip3 install sparsepca
```

To install the R package, install `amanpg` directly from CRAN:

```r
install.packages("amanpg")
```

## Documentation

### Usage

#### Python

```python
spca(z, lambda1, lambda2, 
     x0=None, y0=None, k=0, gamma=0.5, type=0, 
     maxiter=1e4, tol=1e-5, f_palm=1e5,
	 normalize=True, verbose=False):
```

#### R

```r
spca.amanpg(z, lambda1, lambda2, 
			f_palm = 1e5, x0 = NULL, y0 = NULL, k = 0, type = 0, 
			gamma = 0.5, maxiter = 1e4, tol = 1e-5, 
			normalize = TRUE, verbose = FALSE)
```

### Arguments

| Name | Python Type | R Type | Description |
| --- | --- | --- |
| `z` | numpy.ndarray | matrix | Either the data matrix or sample covariance matrix |
| `lambda1` | float list | numeric vector | List of parameters of length n for L1-norm penalty |
| `lambda2` | float or numpy.inf | numeric or Inf | L2-norm penalty term |
| `x0` | numpy.ndarray | matrix | Initial x-values for the gradient method, default value is the first n right singular vectors |
| `y0` | numpy.ndarray | matrix | Initial y-values for the gradient method, default value is the first n right singular vectors |
| `k` | int | int | Number of principal components desired, default is 0 (returns min(n-1, p) principal components) |
| `gamma` | float | numeric | Parameter to control how quickly the step size changes in each iteration, default is 0.5 |
| `type` | int | int | If 0, b is expcted to be a data matrix, and otherwise b is expected to be a covariance matrix; default is 0 |
| `maxiter` | int | int | Maximum number of iterations allowed in the gradient method, default is 1e4 |
| `tol` | float | numeric | Tolerance value required to indicate convergence (calculated as difference between iteration f-values), default is 1e-5 |
| `f_palm` | float | numeric | Upper bound for the F-value to reach convergence, default is 1e5 |
| `normalize` | bool | logical | Center and normalize rows to Euclidean length 1 if True, default is True |
| `verbose` | bool | logical | Function prints progress between iterations if True, default is False |e

### Values

Python returns a dictionary with the following key-value pairs, while R returns a list with the following elements:

| Key | Python Value Type | R Value Type | Value |
| --- | --- | --- |
| `loadings` | numpy.ndarray | matrix | Loadings of the sparse principal components |
| `f_manpg` | float | numeric | Final F-value |
| `x` | numpy.ndarray | matirx | Corresponding ndarray in subproblem to the loadings |
| `iter` | int | numeric | Total number of iterations executed |
| `sparsity` | float | numeric | Number of sparse loadings (loadings == 0) divided by number of all loadings |
| `time` | float | numeric | Execution time in seconds |

## Quick Start

### Python Quick Start

Note that the Python package depends on numpy.

In the following example, the package function is imported first. The appropriate parameters are defined&mdash;in this case, we want four sparse principal components (rank-`k` loadings)&mdash;from a 1000x500 data matrix. The L1-penalty terms are set to 0.1, and the L2-penalty term is set to 1. Note that any positive value can be used for the L2-penalty term, up to `np.inf`. 

A random 1000x500 matrix is generated from the normal distribution, and then the function is called through `spca()`. A printout of the results follows, along with observing the loadings. 

The second example keeps the same parameters except switching `lambda2` with infinity. Again, the results are printed out and the loadings are observed.

```python
import numpy as np
from sparsepca import spca

k = 4  # columns
d = 500  # dimensions
m = 1000  # sample size
lambda1 = 0.1 * np.ones((k, 1))
lambda2 = 1

np.random.seed(10)
a = np.random.normal(0, 1, size=(m, d))  # generate random normal 1000x500 matrix
fin_sprout = spca(a, lambda1, lambda2, k=k)
print(f"Finite: {fin_sprout['iter']} iterations with final value 
		{fin_sprout['f_manpg']}, sparsity {fin_sprout['sparsity']}, 
		timediff {fin_sprout['time']}.")

fin_sprout['loadings']

inf_sprout = spca_amanpg(a, lambda1, np.inf, k=4)
print(f"Infinite: {inf_sprout['iter']} iterations with final value 
		{inf_sprout['f_manpg']}, sparsity {inf_sprout['sparsity']}, 
		timediff {inf_sprout['time']}.")

inf_sprout['loadings']
```

### R Quick Start

In the following example, we load the library using `library(amanpg)` and then define a 1000x500 randomly-generated matrix from the normal distribution. We set the L1-penalty term to 0.1 and L2-penalty term to infinity, and seek the first four principal components.

The default initial point are the `k` right singular vectors from SVD, which we can see manually broken down here. In the function call, we pass the parameters in and output our list sprout. 

The results are printed out, and then we view the loadings.

```r
d <- 500  # dimension
m <- 1000 # sample size
a <- normalize(matrix(rnorm(m * d), m, d))
lambda1 <- 0.1 * matrix(data=1, nrow=4, ncol=1)
x0 <- svd(a, nv=4)$v
sprout <- spca.amanpg(a, lambda1, lambda2=Inf, x0=x0, y0=x0, k=4) 
print(paste(sprout$iter, "iterations,", sprout$sparsity, "sparsity,", sprout$time))

# extract loadings
View(sprout$loadings)
```

## References

Chen, S., Ma, S., Xue, L., and Zou, H. (2020) "An Alternating Manifold Proximal Gradient Method for Sparse Principal Component Analysis and Sparse Canonical Correlation Analysis" INFORMS Journal on Optimization 2:3, 192-208 <[doi:10.1287/ijoo.2019.0032](https://doi.org/10.1287%2Fijoo.2019.0032)>.

Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse principal component analysis. Journal of Computational and Graphical Statistics, 15(2), 265-286 <[doi:10.1198/106186006X113430](https://doi.org/10.1198%2F106186006X113430)>.

Zou, H., & Xue, L. (2018). A selective overview of sparse principal component analysis. Proceedings of the IEEE, 106(8), 1311-1320 <[doi:10.1109/JPROC.2018.2846588](https://doi.org/10.1109%2FJPROC.2018.2846588)>.


