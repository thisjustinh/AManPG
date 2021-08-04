spca.amanpg <- function(z, lambda1, lambda2, f_palm = 1e5, x0 = NULL, y0 = NULL, k = 0, type = 0, gamma = 0.5,
                        maxiter = 1e4, tol = 1e-5, normalize = TRUE, verbose = FALSE) {

  start <- Sys.time()

  if (normalize){
    z <- normalize(z)
  }

  dims <- dim(z);
  m <- dims[1];
  d <- dims[2];
  # anonymous function that gets sum of matrix X times lambda1
  h <- function(x) sum(lambda1 * colSums(abs(x)))

  if (d < m * 2) {
    z <- Conj(t(z)) %*% z
    type <- 1
  }

  svds1 <- svd(z)$d[1]
  if (!type) ly <- 2 * svds1^2 + 2 * lambda2 else ly <- 2 * svds1 + 2 * lambda2

  ### Initialization ###
  if (k==0){
    k <- d
  }
  if (is.null(x0) && is.null(y0)) {  # assign initial point if not provided
    x0 <- svd(z, nv=k)$v
    y0 <- x0
  }

  x <- x0
  y <- y0
  total_linesearch <- 0
  linesearch_flag <- 1
  min_step <- 0
  linesearch_flag_y <- 1
  t <- 1 / ly
  tau <- 100 / d

  if (!type) {
    ay <- Conj(t(z)) %*% (z %*% y)
    ax <- Conj(t(z)) %*% (z %*% x)
  } else {
    ay <- z %*% y
    ax <- z %*% x
  }

  ### Main Loop ###
  if (is.finite(lambda2)) {
    fx <- -2 * sum(x * ay)
    fy <- sum(y * ay) + lambda2 * norm(y, 'F')^2 + h(y)
    f_rgd <- fx + fy

    for (iter in 2:maxiter) {
      if (verbose) {
        iter_start <- Sys.time()
        print("=========================")
        print(paste("On iteration", iter))
      }

      ### Update y ###
      if (!linesearch_flag_y) t <- t * 1.01 else t <- max(1/ly, t/1.01)
      linesearch_flag_y <- 0

      y_t <- prox.l1(y - t*2 * (ay - ax + lambda2 * y), lambda1 * t, k)

      if (!type) ayt <- Conj(t(z)) %*% (z %*% y_t) else ayt <- z %*% y_t
      f_ytrial <- -2 * sum(x * ayt) + sum(y_t * ayt) + lambda2 * norm(y_t, 'F')^2 + h(y_t)
      normpg <- norm(y_t - y, 'F')^2 / t^2

      # adjust step size to be appropriate and recalculate values
      while (f_ytrial > f_rgd[iter - 1] - 1e-3 * t * normpg) {
        t <- gamma * t
        if (t < 1e-5 / d) {
          break
        }

        y_t <- prox.l1(y - t * 2 * (ay - ax + lambda2 * y), lambda1 * t, k)
        if (!type) ayt <- Conj(t(z)) %*% (z*y_t) else ayt <- z %*% y_t
        f_ytrial <- -2 * sum(x * ayt) + sum(y_t * ayt) + lambda2 * norm(y_t, 'F')^2 + h(y_t)
        linesearch_flag_y <- 1
      }

      # assign updated values from loop
      y <- y_t
      ay <- ayt

      ### update x ###
      if (!linesearch_flag) tau <- tau*1.1
      if (min_step == 1) tau <- 1/d
      linesearch_flag <- 0
      min_step <- 0
      gx <- -2 * ay
      xgx <- Conj(t(gx)) %*% x
      rgx <- gx - 0.5 * x %*% (xgx + Conj(t(xgx)))  # Projected gradient

      tx <- x - tau * rgx
      eigendecomp <- eigen(Conj(t(tx)) %*% tx)
      u <- eigendecomp$vectors
      sigma <- eigendecomp$values
      j <- u %*% diag(sqrt(1 / sigma)) %*% Conj(t(u))
      x_trial <- tx %*% j
      f_xtrial <- -2 * sum(x_trial * ay)
      fxval <- -2 * sum(x * ay)
      normpg <- norm(rgx, 'F')^2

      while (f_xtrial > fxval - 1e-3 * tau * normpg) {
        tau <- gamma * tau
        if (tau < 1e-5 / d) {
          min_step <- 1
          break
        }

        tx <- x - tau * rgx
        eigendecomp <- eigen(Conj(t(tx)) %*% tx)
        u <- eigendecomp$vectors
        sigma <- eigendecomp$values
        j <- u %*% diag(sqrt(1 / sigma)) %*% Conj(t(u))
        x_trial <- tx %*% j
        total_linesearch <- total_linesearch + 1
        linesearch_flag <- 1
        f_xtrial <- -2 * sum(x_trial * ay)
      }

      x <- x_trial
      if (!type) ax <- Conj(t(z)) %*% (z %*% x) else ax <- z %*% x
      fx <- f_xtrial
      fy <- sum(y * ay) + lambda2 * norm(y, 'F')^2 + h(y)
      f_rgd <- c(f_rgd, fx + fy)

      if (verbose) {
        print(paste("fx:", fx))
        print(paste("fy:", fy))
        print(paste("Finished with value", f_rgd[iter], "and difference", abs(f_rgd[iter]-f_rgd[iter-1])))
        print(difftime(Sys.time(), iter_start))
      }

      if (abs(f_rgd[iter] - f_rgd[iter - 1]) < tol && f_rgd[iter] < f_palm || abs(f_rgd[iter] - f_rgd[iter - 1]) < 1e-12) {
        if (verbose) print(paste("Difference of", abs(f_rgd[iter] - f_rgd[iter-1])))
        break
      }
    }
  } else {  # lambda2 is Inf
    fx <- -2 * sum(x * ay)
    fy <- norm(y, 'F')^2 + h(y)
    f_rgd <- fx + fy

    for (iter in 2:maxiter) {
      if (verbose) {
        iter_start <- Sys.time()
        print("=========================")
        print(paste("On iteration", iter))
      }

      if (!linesearch_flag) tau <- tau * 1.1
      if (min_step == 1) tau <- 1/d

      min_step <- 0
      linesearch_flag <- 0

      ### Update y ###
      t <- 1/2
      y <- prox.l1(y - 2 * t * (-ax + y), lambda1 * t, k)
      if (!type) ay <- Conj(t(z)) %*% (z %*% y) else ay <- z %*% y

      ### update X ###
      gx <- -2 * ay
      xgx <- Conj(t(gx)) %*% x
      rgx <- gx - x %*% xgx  # Canonical Riemannian gradient
      tx <- x - tau * rgx
      tx.dims <- dim(tx)
      tmp <- svd(tx, nu=min(tx.dims[1], tx.dims[2]), nv=min(tx.dims[1], tx.dims[2]))
      x_trial <- tmp$u %*% Conj(t(tmp$v))
      f_xtrial <- -2 * sum(x_trial * ay)
      fxval <- -2 * sum(x * ay)
      normpg <- norm(rgx, 'F')^2

      while (f_xtrial > fxval - 1e-3 * tau * normpg) {
        tau <- gamma * tau
        if (tau < 1e-3 / d) {
          min_step <- 1
          break
        }

        tx <- x - tau * rgx
        tmp <- svd(tx, nu=min(tx.dims[1], tx.dims[2]), nv=min(tx.dims[1], tx.dims[2]))
        x_trial <- tmp$u %*% Conj(t(tmp$v))
        total_linesearch <- total_linesearch + 1
        linesearch_flag <- 1
        f_xtrial <- -2 * sum(x_trial * ay)
      }

      #assign updated values from loop
      x <- x_trial
      fx <- f_xtrial
      if (!type) ax <- Conj(t(z)) %*% (z %*% x) else ax <- z %*% x
      fy <- norm(y, 'F')^2 + h(y)
      f_rgd <- c(f_rgd, fx + fy)

      if (verbose) {
        print(paste("fx:", fx))
        print(paste("fy:", fy))
        print(paste("Finished with value", f_rgd[iter]))
        print(difftime(Sys.time(), iter_start))
      }

      # if normDsquared < tol^2, as in the e-stationary point
      if (abs(f_rgd[iter] - f_rgd[iter - 1]) < tol) {
        if (verbose) {
          print(paste("Final difference of", abs(f_rgd[iter] - f_rgd[iter - 1])))
        }
        break
      }
    }

  }

  ### Process return list ###
  y_norm <- sqrt(colSums(y^2))
  y_norm[y_norm == 0] <- 1

  results <- list(
    iter=iter,
    f_amanpg=f_rgd[iter],
    sparsity=sum(y == 0) / (d * k),
    time=difftime(Sys.time(), start),
    x=x,
    loadings=y / (matrix(1, d, 1) %*% y_norm)
  )

  return(results)
}


prox.l1 <- function(z, lambda, r) {
  a <- abs(z) - as.vector(lambda)

  if (r < 15) act_set <- as.numeric(a > 0) else act_set <- a > 0
  x_prox <- act_set * sign(z) * a
  inact_set <- a <= 0
  return(x_prox)
}

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


if (sys.nframe() == 0) {
  maxiter <- 1e4
  tol <- 1e-5
  k <- 4  # columns
  d <- 500  # dimension
  m <- 1000  # sample size
  lambda1 <- 0.1 * matrix(data=1, nrow=k, ncol=1)
  type <- 0
  gamma <- 0.5
  lambda2 <- Inf
  f_palm <- 1e5

  test_finite <- FALSE

  if (!test_finite){
    # test inf
    for (i in 1:10) {
      set.seed(i)
      a <- normalize(matrix(rnorm(m * d), m, d))
      x0 <- svd(a, nv=k)$v
      sprout <- spca.amanpg(a, lambda1, lambda2, f_palm, x0, x0, k, type, gamma, maxiter, tol, normalize = FALSE, verbose=FALSE)
      print(paste(sprout$iter, "iterations,", sprout$sparsity, "sparsity,", sprout$time))

      #extract loadings
      #print(sprout$loadings)
    }
  }
  else{
    # test finite
    path <- system.file('extdata', 'a.RDS', package = 'amanpg')
    a <- readRDS(path)
    x0 <- svd(a, nv=k)$v
    sprout <- spca.amanpg(a, lambda1, lambda2=1, f_palm, x0, x0, k, type, gamma, maxiter, tol, normalize = TRUE, verbose=FALSE)
    print(paste(sprout$iter, "iterations,", sprout$sparsity, "sparsity,", sprout$time))

    #extract loadings
    #print(sprout$loadings)
  }
}
