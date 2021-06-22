# mupar <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
# TODO: Ask why demo file doesn't have y0 and F_palm initialized

maxiter <- 1e4
tol <- 1e-5
n <- 4  # columns
d <- 500  # dimension
m <- 1000  # sample size
mu <- 0.1 * matrix(data=1, nrow=n, ncol=1)
type <- 0
# lambda <- Inf
f_palm <- 1e5  # what is this?

set.seed(10)
a <- normalize(matrix(rnorm(m * d), m, d))
x0 <- svd(a, nv=n)$v

sprout <- spca.amanpg(a, mu, lambda=1, n, x0, x0, type, maxiter, tol, f_palm, verbose=TRUE)
sprout <- spca.amanpg(a, mu, lambda=Inf, n, x0, x0, type, maxiter, tol, f_palm, verbose=TRUE)
