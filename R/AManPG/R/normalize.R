normalize <- function(x, center=TRUE, scale=TRUE) {
  if (center) {
    x <- scale(x, scale=FALSE)
  }

  if (scale) {
    dims <- dim(x)
    rows <- sqrt(rowSums(x^2))
    x <- x / matrix(rows, nrow=dims[1], ncol=dims[2])
  }

  return(x)
}
