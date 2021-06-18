# if option.lambda < inf,  min Tr(Y'*B'B*Y)- 2Tr(X'*B'*B*Y) + lambda*norm(Y,'fro')^2  + mu*norm(Y,1) s.t. X'*X=In.
# if option.lambda = inf,  min - 2Tr(X'*B'*B*Y)+ norm(Y,'fro')^2  + mu*norm(Y,1) s.t. X'*X=In.

# TODO: Check that matrix multiplication is correct. Use crossprod when possible.
# TODO: Check if conjugate transpose is necessary. If not, get rid of Conj() calls.

spca_amanpg = function(B, mu, lambda, n, type, maxiter, tol, X0, Y0, F_palm) {
  start <- Sys.time()

  dims <- dim(B);
  m <- dims[1];
  b <- dims[2];
  # anonymous function that gets sum of matrix X times mu
  # TODO: Check that this matches MATLAB functionality (should colSums be used?)
  h <- function(x) mu * sum(x)

  if (d < m * 2) {
    # TODO: Is conjugate transpose necessary?
    B <- Conj(t(B)) %*% B  # B^T B
    type <- 1
  }

  # TODO: Check that the svds replacement works
  svds1 <- svd(B)$d[1]
  if (!type) LY <- 2 * svds1^2 + 2 * lambda else LY <- 2 * svds1 + 2 * lambda

  # initial point
  X <- X0
  Y <- Y0
  total_linesearch <- 0
  linesearch_flag <- 1
  min_step <- 0
  linesearch_flag_Y <- 1
  t <- 1 / LY
  tau <- 100 / d

  if (!type) {
    AY <- Conj(t(B)) %*% (B %*% Y)
    AX <- Conj(t(B)) %*% (B %*% X)
  } else {
    AY <- B %*% Y
    AX <- B %*% X
  }

  if (is.finite(lambda)) {
    fx <- -2 * sum(X * AY)
    fy <- sum(Y * AY) + lambda * norm(Y, 'F')^2 + h(Y)
    F_rgd <- fx + fy  # will become vector in the loop

    for (iter in 2:maxiter) {
      ### update Y ###
      if (!linesearch_flag_Y) t <- t * 1.01 else t <- max(1/LY, t/1.01)
      linesearch_flag_Y <- 0
      Y_t <- prox_l1(Y - t*2 * (AY - AX + lambda * Y), mu * t, n)
      if (!type) AYt <- Conj(t(B)) %*% (B %*% Y_t) else AYt <- B %*% Y_t
      f_ytrial <- -2 * sum(X * AYt) + sum(Y_t * AYt) + lambda * norm(Y_t, 'F')^2 + h(Y_t)
      normpg <- norm(Y_t - Y, 'F')^2 / t^2

      while (f_ytrial > F_rgd(iter - 1) - 1e-3 *t * normpg) {
        # TODO: Consider custom gamma value?
        t <- 0.5 * t  # gamma = 0.5
        if (t < 1e-5 / d) {
          break
        }
        Y_t <- prox_l1(Y - t * 2 * (AY - AX + lambda * Y), mu * t, n)
        if (!type) AYt <- Conj(t(B)) %*% (B*Y_t) else AYt <- B %*% Y_t
        f_ytrial <- -2 * sum(X * AYt) + sum(Y_t * AYt) + lambda * norm(Y_t, 'F')^2 + h(Y_t)
        linesearch_flag_Y <- 1
      }

      Y <- Y_t
      AY <- AYt

      ### update X ###
      if (!linesearch_flag) tau <- tau*1.1 else tau <- 1/d
      linesearch_flag <- 0
      min_step <- 0
      gx <- -2 * AY
      xgx <- Conj(t(gx)) %*% X
      # RGX <- gx - X %*% xgx  # Canonical Riemannian gradient
      RGX <- gx - 0.5 * X %*% (xgx + Conj(t(xgx)))  # Projected gradient
      TX <- X - tau %*% RGX
      eigendecomp <- eigen(Conj(t(TX)) %*% TX)
      U <- diag(x = eigendecomp$values)
      SIGMA <- diag(SIGMA)
      J <- U %*% diag(sqrt(1 / SIGMA)) %*% Conj(t(U))
      X_trial <- TX %*% J
      f_xtrial <- -2 * sum(X_trial * AY)
      fXval <- -2 * sum(X * AY)
      normpg <- norm(RGX, 'F')^2

      while (f_xtrial > fXval - 1e-3 * tau * normpg) {
        # TODO: Configurable gamma parameter or nah?
        tau <- 0.5 * tau
        if (tau < 1e-5 / d) {
          min_step <- 1
          break
        }

        TX <- X - tau %*% RGX
        eigendecomp <- eigen(Conj(t(TX)) %*% TX)
        U <- diag(x = eigendecomp$values)
        SIGMA <- diag(SIGMA)
        J <- U %*% diag(sqrt(1 / SIGMA)) %*% Conj(t(U))
        X_trial <- TX %*% J
        total_linesearch <- total_linesearch + 1
        linesearch_flag <- 1
        f_xtrial <- -2 * sum(X_trial * AY)
      }

      X <- X_trial
      if (!type) AX <- Conj(t(B)) %*% (B %*% X) else AX <- B %*% X
      fx <- f_xtrial
      fy <- sum(Y * AY) + lambda * norm(Y, 'F')^2 + h(Y)
      F_rgd <- c(F_rgd, fx + fy)  # concatenate vector

      if (iter > 1)
        if ((abs(F_rgd[iter]) - F_rgd[iter - 1]) < tol && F_rgd[iter] < F_palm || abs(F_rgd[iter] - F_rgd[iter - 1]) < 1e-12) break
    }
  } else {  # lambda is Inf
    fx <- -2 * sum(X * AY)
    fy <- norm(Y, 'F')^2 + h(Y)
    F_rgd <- fx + fy  # will become vector

    for (iter in 2:maxiter) {
      if (!linesearch_flag) tau <- tau * 1.1
      if (min_step == 1) tau <- 1/d

      min_step <- 0
      linesearch_flag <- 0

      ### update Y ###
      t <- 1/2
      Y <- prox_l1(Y - 2 * 2 * (-AX + Y), mu * t, n)
      if (!type) AY <- Conj(t(B)) %*% (B %*% Y) else AY <- B %*% Y

      ### update X ###
      gx <- -2 * AY
      xgx <- Conj(t(gx)) %*% X
      RGX <- gx - X %*% xgx  # Canonical Riemannian gradient
      TX <- X - tau * RGX
      tmp <- svd(TX)  # TODO!!!!!: Figure out economy-sized decomposition
      X_trial <- tmp$u %>% tmp$v
      f_xtrial <- -2 * sum(X * AY)
      fXval <- -2 * sum(X * AY)
      normpg <- norm(RGX, 'F')^2
      while (f_xtrial > fXval - 1e-3 * tau * normpg) {
        # TODO: Configure gamma?
        tau <- 0.5 * tau
        if (tau < 1e-3 / d) {
          min_step <- 1
          break
        }

        TX <- X - tau * RGX
        tmp <- svd(TX)  # TODO!!!!!: Figure out economy-sized decomposition
        X_trial <- tmp$u %>% tmp$v
        total_linesearch <- total_linesearch + 1
        linesearch_flag <- 1
        f_xtrial <- -2 * sum(X * AY)
      }

      X <- X_trial
      if (!type) AX <- Conj(t(B)) %*% (B %*% X) else AX <- B %*% X
      fx <- f_xtrial
      fy <- norm(Y, 'F')^2 + h(Y)
      F_rgd <- c(F_rgd, fx + fy)

      # if normDsquared < tol^2, as in the e-stationary point
      if (iter > 1)
        if (abs(F_rgd[iter] - F_rgd[iter - 1]) < tol)
          break
    }
  }

  ### Process return list ###
  Y_norm <- sqrt(colSums(Y^2))  # TODO: Check colSums?
  Y_norm[Y_norm == 0] <- 1

  results <- list(
    iter=iter,
    F_amanpg=F_rgd[iter],
    sparsity=sum(Y == 0) / (d * n),
    time=difftime(Sys.time(), start),  # TODO: Check if acceptable value of time
    X=X,
    Y_man=Y / (matrix(0, d, 1) %*% Y_norm)
  )

  return(results)
}
