scca.amanpg <- function(a, b, n, p, q, b1, b2, maxiter, tol, inner_tol, sp_type) {
  start <- Sys.time()

  # initial value assignment
  m <- dim(A)[1]  # TODO: pretty sure this should be number of rows since "sample size" but should prob check
  tau1 <- b1 * sqrt((n + log(p)) / m)
  tau2 <- b2 * sqrt((n + log(p)) / m)


  # define anonymous functions
  row.norms <- function(x) apply(x, 1, function(y) sqrt(sum(y ^ 2)))

  if (sp_type <- 'l1') {  # uses 1-norm using abs
    h1 <- function(x) tau1 * sum(abs(x))
    h2 <- function(x) tau2 * sum(abs(x))
    prox.func <- function(b, lambda, r) prox.l1(b, lambda, r)
  } else if (sp_type <- 'l21') {  # uses 2-norm
    h1 <- function(x) tau1 * sum(row.norms(x))  # TODO: this is almost certainly wrong lol
    h2 <- function(x) tau2 * sum(row.norms(x))
    prox.func <- function(b, lambda, r) prox.l21(b, lambda, r)
  }

  inner_flag1 <- 0
  # TODO: Finish this
}
