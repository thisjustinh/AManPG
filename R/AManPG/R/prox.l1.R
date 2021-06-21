prox.l1 <- function(b, lambda, r) {
  #print(b)
  #print(lambda)
  #print(r)
  a <- abs(b) - lambda

  # TODO: Check what the heck the difference between these is
  if (r < 15) act_set <- as.numeric(a > 0) else act_set = a > 0
  x_prox <- act_set * sign(b) * a
  inact_set <- a <= 0
  return(x_prox)
}
