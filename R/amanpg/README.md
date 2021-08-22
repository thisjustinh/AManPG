# amanpg

## Description

Uses an alternating manifold proximal gradient (AManPG) method to find sparse principal components from the given data or covariance matrix. 

Only base R is required to be installed.

## Usage
```R
spca.amanpg(z, lambda1, lambda2, f_palm = 1e5, x0 = NULL, y0 = NULL, k = 0, type = 0, gamma = 0.5,
                        maxiter = 1e4, tol = 1e-5, normalize = TRUE, verbose = FALSE)
```


## Arguments

| Name | Type | Description |
| --- | --- | --- |
| `z` | matrix | Either the data matrix or sample covariance matrix |
| `lambda1` | matrix | List of parameters of length n for L1-norm penalty |
| `lambda2` | double| L2-norm penalty term |
| `f_palm` | double | Upper bound for the gradient value to reach convergence, default value is 1e5 |
| `x0` | matrix | Initial x-values for the gradient method, default value is the first n right singular vectors |
| `y0` | matrix | Initial y-values for the gradient method, default value is the first n right singular vectors |
| `k` | integer | Number of principal components desired, default is 0 (returns min(n-1, p) principal components) |
| `type` | integer | If 0, b is expected to be a data matrix, and otherwise b is expected to be a covariance matrix; default is 0 |
| `gamma` | double | Parameter to control how quickly the step size changes in each iteration, default is 0.5 |
| `maxiter` | integer | Maximum number of iterations allowed in the gradient method, default is 1e4 |
| `tol` | double | Tolerance value required to indicate convergence (calculated as difference between iteration f-values), default is 1e-5 |
| `normalize` | logical | Center and normalize rows to Euclidean length 1 if True, default is True |
| `verbose` | logical | Function prints progress between iterations if True, default is False |

## Value

Returns a dictionary with the following key-value pairs:

| Key | Value Type | Value |
| --- | --- | --- |
| `iter` | integer | Total number of iterations executed |
| `f_manpg` | double | Final gradient value |
| `sparsity` | float | Number of sparse loadings (loadings == 0) divided by number of all loadings |
| `time` | double | Number of seconds for execution |
| `x` | matrix | Corresponding ndarray in subproblem to the loadings |
| `loadings` | matrix | Loadings of the sparse principal components |



## Authors
 
Shixiang Chen, Justin Huang, Benjamin Jochem, Shiqian Ma, Lingzhou Xue and Hui Zou

## References

Chen, S., Ma, S., Xue, L., and Zou, H. (2020) "An Alternating Manifold Proximal Gradient Method for Sparse Principal Component Analysis and Sparse Canonical Correlation Analysis" *INFORMS Journal on Optimization* 2:3, 192-208

Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse principal component analysis. 
Journal of Computational and Graphical Statistics, 15(2), 265-286.

Zou, H., & Xue, L. (2018). A selective overview of sparse principal component analysis. 
Proceedings of the IEEE, 106(8), 1311-1320.

## Example

See `SPCA.R` for a more in-depth example.

```R
library('SPCA')

#see SPCA.R for a more in-depth example
      d <- 500  # dimension
      m <- 1000 # sample size
      set.seed(10)
      a <- normalize(matrix(rnorm(m * d), m, d))
      lambda1 <- 0.1 * matrix(data=1, nrow=4, ncol=1)
      x0 <- svd(a, nv=4)$v
      sprout <- spca.amanpg(a, lambda1, lambda2=Inf, f_palm=1e5, x0=x0, y0=x0, k=4, type=0, gamma=0.5,
                            maxiter=1e4, tol=1e-5, normalize = FALSE, verbose=FALSE)
      print(paste(sprout$iter, "iterations,", sprout$sparsity, "sparsity,", sprout$time))

      #extract loadings
      #print(sprout$loadings)
```
