prox_l1 <- function(b, lambda, r) {
  a <- abs(b) - lambda
  if (r < 15) act_set <- as.numeric(a > 0) else act_set = a > 0
  x_prox <- act_set * sign(b) * a
  inact_set <- a <= 0
  return(x_prox)
}
