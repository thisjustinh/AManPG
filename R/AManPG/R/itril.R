itril <- function(sz, k=0) {
  if (is.atomic(sz) && length(sz) == 1L) {
    sz = c(sz, sz)
  }

  m = sz[1]
  n = sz[2]

  nc = min(n, m + k)  # number of columns of the triangular part
  lo = max((1:nc) - k, 1)  # lower row index for each column
  hi = m + matrix(0, nc, 1)  # upper row index for each column

  if (!length(lo)) {
    i = matrix(0, 0, 1)  # zeros(0, 1)
    j = matrix(0, 0, 1)
  } else {
    c = cumsum((0:hi-lo) + 1)
    # TODO: Finish this
  }
}
