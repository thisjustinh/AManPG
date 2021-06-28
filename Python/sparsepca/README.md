# SparsePCA

## Description

Uses an alternating manifold proximal gradient (AManPG) method to find sparse principal components from the given data.

## Usage

```python
spca_amanpg(b, mu, lamb, n, f_palm, 
            x0=None, y0=None, gamma=0.5, type=0, 
            maxiter=1e4, tol=1e-5, verbose=False):
```

## Arguments

| Name | Type | Description |
| --- | --- | --- |
| b | numpy.ndarray | Either the data matrix or sample covariance matrix |
| mu | float list | List of parameters of length n for L1-norm penalty |
| lamb | float or numpy.inf | L2-norm penalty term |
| n | int | Number of principal components desired |
| f_palm | float | Upper bound for the gradient value to reach convergence |
| x0 | numpy.ndarray | Initial x-values for the gradient method, default value is the first n right singular vectors |
| y0 | numpy.ndarray | Initial y-values for the gradient method, default value is the first n right singular vectors |
| gamma | float | Parameter to control how quickly the step size changes in each iteration, default is 0.5 |
| type | int | If 0, b is expected to be a data matrix, and otherwise b is expected to be a covariance matrix; default is 0 |
| maxiter | int | Maximum number of iterations allowed in the gradient method, default is 1e4 |
| tol | float | Tolerance value required to indicate convergence (calculated as difference between iteration f-values), default is 1e-5 |
| verbose | bool | Function prints progress between iterations if True, default is False |

## Value

Returns a dictionary with the following key-value pairs:

| Key | Value Type | Value |
| --- | --- | --- |
| loadings | numpy.ndarray | Loadings of the sparse principal components |
| f_manpg | float | Final gradient value |
| x | numpy.ndarray | Corresponding ndarray in subproblem to the loadings |
| iter | int | Total number of iterations executed |
| sparsity | float | Number of sparse loadings (loadings == 0) divided by number of all loadings |
| time | Number of seconds for execution |

## Authors
 
Justin Huang and Benjamin Jochem

## References

Chen, S., Ma, S., Xue, L., and Zou, H. (2020) "An Alternating Manifold Proximal Gradient Method for Sparse Principal Component Analysis and Sparse Canonical Correlation Analysis" *INFORMS Journal on Optimization* 2:3, 192-208

## Example

