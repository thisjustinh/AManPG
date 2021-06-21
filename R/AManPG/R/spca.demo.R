# mupar <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
# TODO: Ask why demo file doesn't have y0 and F_palm initialized

maxiter <- 1e4
tol <- 1e-5
n <- 4  # columns
d <- 500  # dimension
m <- 1000  # sample size
mu <- 0.1 * matrix(data=1, nrow=n, ncol=1)
type <- 0
lambda <- 1

set.seed(10)
a <- matrix(rnorm(m * d), m, d)
a <- scale(a)
x0 <- svd(a, nv=n)$v

sprout <- spca_amanpg(a, mu, lambda, n, type, maxiter, tol, x0, x0, 0)
