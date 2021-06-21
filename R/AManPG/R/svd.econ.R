svd.econ <- function(x, type=1) {
  dims <- dim(x)

  if (dims[1] > dims[2]) {  # m > n
    return(svd(x, nu=dims[2]))
  } else if (dims[1] < dims[2] && type) {  # m < n
    return(svd(x, nv=dims[1]))
  } else {  # m == n
    return(svd(x))
  }
}
