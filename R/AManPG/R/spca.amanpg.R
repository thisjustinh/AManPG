# if option.lambda < inf,  min Tr(Y'*B'B*Y)- 2Tr(X'*B'*B*Y) + lambda*norm(Y,'fro')^2  + mu*norm(Y,1) s.t. X'*X=In.
# if option.lambda = inf,  min - 2Tr(X'*B'*B*Y)+ norm(Y,'fro')^2  + mu*norm(Y,1) s.t. X'*X=In.
# TODO: Check if conjugate transpose is necessary. If not, get rid of Conj() calls.
# TODO: Consider custom gamma value?

spca.amanpg <- function(b, mu, lambda, n, x0, y0, f_palm, type=0, gamma=0.5,
                        maxiter=1e4, tol=1e-5, verbose=FALSE) {

  start <- Sys.time()

  dims <- dim(b);
  m <- dims[1];
  d <- dims[2];
  # anonymous function that gets sum of matrix X times mu
  h <- function(x) sum(mu * colSums(abs(x)))

  if (d < m * 2) {
    b <- Conj(t(b)) %*% b  # B^T B
    type <- 1
  }

  svds1 <- svd(b)$d[1]
  if (!type) ly <- 2 * svds1^2 + 2 * lambda else ly <- 2 * svds1 + 2 * lambda

  # initialization
  if (missing(x0) && missing(y0)) {  # assign initial point if not provided
    x0 <- svd(b, nv=n)$v
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
    ay <- Conj(t(b)) %*% (b %*% y)
    ax <- Conj(t(b)) %*% (b %*% x)
  } else {
    ay <- b %*% y
    ax <- b %*% x
  }

  ### Main Loop ###
  if (is.finite(lambda)) {  # elastic net parameter isn't Inf
    fx <- -2 * sum(x * ay)
    fy <- sum(y * ay) + lambda * norm(y, 'F')^2 + h(y)
    f_rgd <- fx + fy  # will become vector in the loop

    for (iter in 2:maxiter) {
      if (verbose) {
        iter_start <- Sys.time()
        print("=========================")
        print(paste("On iteration", iter))
      }

      ### update y ###
      if (!linesearch_flag_y) t <- t * 1.01 else t <- max(1/ly, t/1.01)
      linesearch_flag_y <- 0

      y_t <- prox.l1(y - t*2 * (ay - ax + lambda * y), mu * t, n)

      if (!type) ayt <- Conj(t(b)) %*% (b %*% y_t) else ayt <- b %*% y_t
      f_ytrial <- -2 * sum(x * ayt) + sum(y_t * ayt) + lambda * norm(y_t, 'F')^2 + h(y_t)
      normpg <- norm(y_t - y, 'F')^2 / t^2

      # adjust step size to be appropriate
      while (f_ytrial > f_rgd[iter - 1] - 1e-3 * t * normpg) {
        t <- 0.5 * t  # gamma = 0.5
        if (t < 1e-5 / d) {
          break
        }

        y_t <- prox.l1(y - t * 2 * (ay - ax + lambda * y), mu * t, n)
        if (!type) ayt <- Conj(t(b)) %*% (b*y_t) else ayt <- b %*% y_t
        f_ytrial <- -2 * sum(x * ayt) + sum(y_t * ayt) + lambda * norm(y_t, 'F')^2 + h(y_t)
        linesearch_flag_y <- 1
      }

      y <- y_t
      ay <- ayt

      ### update x ###
      if (!linesearch_flag) tau <- tau*1.1
      if (min_step == 1) tau <- 1/d
      linesearch_flag <- 0
      min_step <- 0
      gx <- -2 * ay
      xgx <- Conj(t(gx)) %*% x
      # RGX <- gx - X %*% xgx  # Canonical Riemannian gradient
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
        tau <- 0.5 * tau
        if (tau < 1e-5 / d) {
          min_step <- 1
          break
        }

        tx <- x - tau * rgx
        eigendecomp <- eigen(Conj(t(tx)) %*% tx)
        u <- eigendecomp$vectors
        sigma <- eigendecomp$values
        j <- u %*% diag(sqrt(1 / sigma)) %*% Conj(t(u))  # note use of abs
        x_trial <- tx %*% j
        total_linesearch <- total_linesearch + 1
        linesearch_flag <- 1
        f_xtrial <- -2 * sum(x_trial * ay)
      }

      x <- x_trial
      if (!type) ax <- Conj(t(b)) %*% (b %*% x) else ax <- b %*% x
      fx <- f_xtrial
      fy <- sum(y * ay) + lambda * norm(y, 'F')^2 + h(y)
      f_rgd <- c(f_rgd, fx + fy)  # concatenate vector

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
  } else {  # lambda is Inf
    fx <- -2 * sum(x * ay)
    fy <- norm(y, 'F')^2 + h(y)
    f_rgd <- fx + fy  # will become vector

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

      ### update y ###
      t <- 1/2
      y <- prox.l1(y - 2 * t * (-ax + y), mu * t, n)
      if (!type) ay <- Conj(t(b)) %*% (b %*% y) else ay <- b %*% y

      ### update X ###
      gx <- -2 * ay
      xgx <- Conj(t(gx)) %*% x
      rgx <- gx - x %*% xgx  # Canonical Riemannian gradient
      tx <- x - tau * rgx
      # tmp <- svd.econ(tx, type=0)
      tx.dims <- dim(tx)
      tmp <- svd(tx, nu=min(tx.dims[1], tx.dims[2]), nv=min(tx.dims[1], tx.dims[2]))
      x_trial <- tmp$u %*% Conj(t(tmp$v))
      f_xtrial <- -2 * sum(x_trial * ay)
      fxval <- -2 * sum(x * ay)
      normpg <- norm(rgx, 'F')^2

      while (f_xtrial > fxval - 1e-3 * tau * normpg) {
        tau <- 0.5 * tau
        if (tau < 1e-3 / d) {
          min_step <- 1
          break
        }

        tx <- x - tau * rgx
        # tmp <- svd.econ(tx, type=0)
        tmp <- svd(tx, nu=min(tx.dims[1], tx.dims[2]), nv=min(tx.dims[1], tx.dims[2]))
        x_trial <- tmp$u %*% Conj(t(tmp$v))
        total_linesearch <- total_linesearch + 1
        linesearch_flag <- 1
        f_xtrial <- -2 * sum(x_trial * ay)
      }

      x <- x_trial
      if (!type) ax <- Conj(t(b)) %*% (b %*% x) else ax <- b %*% x
      fx <- f_xtrial
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
    sparsity=sum(y == 0) / (d * n),
    time=difftime(Sys.time(), start),  # TODO: Check if acceptable value of time
    x=x,
    y_man=y / (matrix(1, d, 1) %*% y_norm)
  )

  return(results)
}
