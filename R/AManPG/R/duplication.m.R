duplication.m <- function(n, option=0) {
  if (is.atomic(n) && length(n) == 1L) {
    n = c(n, n)
  }

  if (!option) {  # lower
    # TODO: Finish this
  } else {  # upper

  }
}
