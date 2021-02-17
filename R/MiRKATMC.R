#' The Microbiome Regression-based Kernel Association Test for multi-categorical
#' (nominal or ordinal) data
#'
#'  \code{MiRKATMC} is able to test for the association between multi-categorical
#'  outcomes (nominal and ordinal) and microbiome, for both independent and
#'  clustered/longitudinal data.
#'
#' @import stats PearsonDS
#' @importFrom Matrix Diagonal Matrix
#' @param formula The model formula with fixed effects only. The response should be
#' a factor.
#' @param data.type Either 'nominal' or 'ordinal'.
#' @param Ks A kernel matrix or a list of kernels.
#' @param data An optional data frame (default is NULL).
#' @param random An optional formula for the random effect (default is NULL).
#' Currently, only models with random intercepts or with random intercepts and
#' slopes are supported.
#'
#' @return A p value or a list of p values corresponding to multiple kernels. If
#' multiple kernels provided, the omnibus test to combine multiple p values
#' will be conducted and returns the p value as \code{omnibus.HMP}.
#' @export
#'
#' @author
#' Zhiwen Jiang
#'
#' @examples
#' test.data <- data.frame(outcome = as.factor(sample(4, 100, replace = TRUE)),
#' ID = gl(20, 5), time = rep(1:5, 20), age = rnorm(n = 100, mean = 30, sd = 5),
#' sex = rbinom(100, 1, 1/2))
#' D <- matrix(rbinom(10000, 2, 0.05), 100, 100)
#' K <- crossprod(D) # kernel matrix
#'
#' # independent nominal data
#' MiRKATMC(formula = outcome ~ sex + age, data = test.data, data.type = 'nominal',
#' Ks = K)
#'
#' # clustered nominal data
#' MiRKATMC(formula = outcome ~ sex + age, data = test.data, data.type = 'nominal',
#' Ks = K, random = ~ 1 | ID)
#'
#' # longitudinal nominal data
#' MiRKATMC(formula = outcome ~ sex + age, data = test.data, data.type = 'nominal',
#' Ks = K, random = ~ 1 + time | ID)
MiRKATMC <- function(formula, data.type, Ks, random = NULL, data = NULL){
  # browser()

  if (is.matrix(Ks)){
    Ks <- list(Ks)
  } else if (!is.list(Ks)){
    stop("Please enter either a single kernel matrix or a list of kernels for Ks.")
  }

  # if (is.null(names(Ks)) & length(Ks) != 1){
  #   warning('No kernel names provided in the list of kernels')
  # }

  if (data.type == 'nominal'){
    if (is.null(random)){
      KY <- inner.GLM.null(formula, data)
    } else{
      KY <- inner.GLM.random.null(formula, random, data)
    }
  } else if (data.type == 'ordinal'){
    if (is.null(random)){
      KY <- inner.POM.null(formula, data)
    } else{
      KY <- inner.POM.random.null(formula, random, data)
    }
  } else{
    stop('data.type should be either \"nominal\" or \"ordinal\".')
  }

  pvs <- mapply(function(x){DKAT(x, KY)}, Ks)
  # browser()
  if (length(pvs) > 1){
    p.value.combined <- harmonicmeanp::p.hmp(p = pvs, w = NULL, L = length(pvs),
                                  w.sum.tolerance = 1e-6, multilevel = F)[1]
    pvs['omnibus.HMP'] <- p.value.combined
  }
  return(pvs)
}




inner.GLM.null <- function(formula, data = NULL){
  # browser()
  # debug(mclogit::mblogit)
  invisible(utils::capture.output(fit <- mclogit::mblogit(formula, data, na.action = 'na.fail')))
  X <- model.matrix(formula, fit$model)
  nSam <- dim(X)[1]
  y.matrix <- matrix(fit$y, nrow = nSam, byrow = T)
  J <- ncol(y.matrix)

  mu.hat <- matrix(fit$fitted.values, ncol = J, byrow = T)
  D <- c((mu.hat[, 2:J] + mu.hat[, 1]) / (mu.hat[, 2:J] * mu.hat[, 1]))
  mu.hat <- c(mu.hat[, 2:J])
  V <- Matrix(0, nrow = nSam * (J-1), ncol = nSam)
  V[cbind(c(1:(nSam * (J-1))), rep(1:nSam, J-1))] <- mu.hat
  V <- Diagonal(x = mu.hat) - Matrix::tcrossprod(V)
  VI <- chol2inv(chol(V))
  Y <- matrix(1/D * VI %*% (c(y.matrix[, 2:J]) - mu.hat), ncol = J - 1)
  KY <- tcrossprod(Y)
  return(KY)
}




inner.GLM.random.null <- function(formula, random, data = NULL){
  # browser()
  arg <- list(formula = formula, data = data, random = random, na.action = 'na.fail',
              method = 'MQL')
  # debug(mclogit::mblogit)
  cat('Fitting the null model...\n')
  fit <- do.call(mclogit::mblogit, arg)
  # browser()
  if (is.null(data)) data <- fit$model
  X <- model.matrix(formula, data)
  nSam <- dim(X)[1]
  y.matrix <- matrix(fit$y, nrow = nSam, byrow = T)
  J <- ncol(y.matrix)

  random <- setupRandomFormula(random)
  rt <- terms(random$formula)
  Z <- model.matrix(rt, data)
  Z <- get_Z(Z, data[, random$groups], nSam)

  if (as.character(random$formula)[[2]] == '1'){
    Vb <- fit$VarCov$`1`
    ZVbZ <- kronecker(Vb, Matrix::tcrossprod(Z))
  }
  else {
    Vb <- fit$VarCov$`1`[order(rep(seq(1, J-1), 2)), order(rep(seq(1, J-1), 2))]
    to.split <- kronecker(matrix(1:(J-1)^2, J-1, byrow = TRUE), matrix(1, 2, 2))
    Vb.list <- lapply(split(Vb, to.split), matrix, nr = 2)
    nclusters <- length(unique(data[, random$groups]))
    Vb.list <- mapply(function(x){kronecker(Diagonal(nclusters), x)}, Vb.list, SIMPLIFY = F)
    sqrt.len <- sqrt(length(Vb.list))
    Vb.m <- do.call(cbind, Vb.list[1: sqrt.len])
    Vb.list <- Vb.list[-(1: sqrt.len)]
    while(length(Vb.list) > 0){
      Vb.mi <- do.call(cbind, Vb.list[1: sqrt.len])
      Vb.m <- rbind(Vb.m, Vb.mi)
      Vb.list <- Vb.list[-(1: sqrt.len)]
    }
    Z.k <- kronecker(Diagonal(J-1), Z)
    ZVbZ <- Z.k %*% Vb.m %*% Matrix::t(Z.k)
  }
  # rf <- c(matrix(unlist(fit$random.effects), ncol = J - 1, byrow = T))
  # Zb.hat <- kronecker(diag(1, J - 1), Z) %*% rf

  # X.matrix <- kronecker(Diagonal(J-1), X)
  mu.hat <- matrix(fit$fitted.values, ncol = J, byrow = T)
  # eta.hat <- matrix(fit$linear.predictors, ncol = J, byrow = T)
  D <- c((mu.hat[, 2:J] + mu.hat[, 1]) / (mu.hat[, 2:J] * mu.hat[, 1]))
  mu.hat <- c(mu.hat[, 2:J])
  V <- Matrix(0, nrow = nSam * (J-1), ncol = nSam)
  V[cbind(c(1:(nSam * (J-1))), rep(1:nSam, J-1))] <- mu.hat
  V <- Diagonal(x = mu.hat) - Matrix::tcrossprod(V)
  res <- D * (c(y.matrix[, 2:J]) - mu.hat) # MQL
  # res <- D * c(y.matrix[, 2:J] - mu.hat[, 2:J]) + c(eta.hat[, 2:J]) -
  #   X.matrix %*% c(matrix(fit$coefficients, ncol = J - 1, byrow = T))
  WI <- Matrix::t(D * Matrix::t(D * V))
  Sigma <- WI + ZVbZ
  SigmaI <- chol2inv(chol(Sigma))
  Y <- matrix(SigmaI %*% res, ncol = J - 1)
  KY <- tcrossprod(Y)
  return(KY)
}






inner.POM.null <- function(formula, data = NULL){
  # browser()
  fit <- MASS::polr(formula, data, method = 'logistic', na.action = 'na.fail')
  if (is.null(data)) data <- fit$model
  X <- model.matrix(formula, data)
  nSam <- dim(X)[1]
  y <- as.formula(paste0('~ 0 +', as.character(formula)[[2]]))
  y.matrix <- model.matrix(y, data)
  J <- ncol(y.matrix)

  pi.hat <- fit$fitted.values
  A <- matrix(1, nrow = J-1, ncol = J-1)
  A[upper.tri(A)] <- 0
  A <- kronecker(A, Diagonal(nSam, 1))

  mu.hat <- matrix(A %*% c(pi.hat[, 1: (J-1)]), ncol = J-1)
  y.tilde <- matrix(A %*% c(y.matrix[, 1: (J-1)]), ncol = J-1)
  KY <- tcrossprod(y.tilde - mu.hat)

  return(KY)
}


inner.POM.random.null <- function(formula, random, data = NULL){
  # browser()
  rf <- random
  rf[[1]] <- NULL
  rf <- sprintf('(%s)', rf)
  rf <- paste(c(".~.", rf), collapse="+")
  rf <- as.formula(rf)
  formula.model <- update(formula, rf)
  fit <- ordinal::clmm(formula.model, data, link = 'logit', threshold = 'flexible',
                       na.action = 'na.fail')

  if (is.null(data)) data <- fit$model
  X <- model.matrix(formula, data)
  nSam <- dim(X)[1]
  y <- as.formula(paste0('~ 0 +', as.character(formula)[[2]]))
  y.matrix <- model.matrix(y, data)
  J <- ncol(y.matrix)

  # if (is.null(fit$Zt)){
  #   random <- setupRandomFormula(random)
  #   rt <- terms(random$formula)
  #   Z <- model.matrix(rt, data)
  #   Z <- get_Z(Z, data[, random$groups], nSam)
  # } else {
  #   Z <- Matrix::t(fit$Zt)
  # }

    random <- setupRandomFormula(random)
    rt <- terms(random$formula)
    Z <- model.matrix(rt, data)
    Z <- get_Z(Z, data[, random$groups], nSam)
    nlev <- fit$dims$nlev.gf

  # browser()
  ranef <- c(t(matrix(fit$ranef, nrow = nlev)))
  # Zb.hat <- rep(Z %*% fit$ranef, J - 1)
  Zb.hat <- rep(Z %*% ranef, J - 1)
  if(is.null(dim(fit$condVar))){
    Vb <- Diagonal(nlev, fit$ST[[1]]^2)
  } else{
    # Vb <- fit$condVar
    ST <- fit$ST[[1]]
    T.matrix <- diag(2)
    T.matrix[2,1] <- ST[2,1]
    S.matrix <- diag(diag(ST))
    Vb <- kronecker(diag(nlev), tcrossprod(T.matrix %*% S.matrix))
  }
  ZVbZ <- kronecker(Diagonal(J - 1), Z %*% Vb %*% Matrix::t(Z))

  X.new <- cbind(kronecker(diag(J - 1), rep(1, nSam)),
                 kronecker(rep(1, J - 1), -X[, -1]))
  eta.hat <- X.new %*% fit$coefficients - Zb.hat
  mu.hat <- matrix(exp(eta.hat)/(1 + exp(eta.hat)), ncol = J-1)

  # browser()
  # start1 <- Sys.time()
  # V <- get_V_POM(mu.hat, J = J, nSam)
  # D <- quickInverse(V, dimen = J, nSam)
  # D.matrix <- get_matrix_W(D, J = J)
  # end1 <- Sys.time()
  # end1 - start1
  # browser()

  A <- matrix(1, nrow = J-1, ncol = J-1)
  A[upper.tri(A)] <- 0
  A <- kronecker(A, Diagonal(nSam, 1))
  y.tilde <- matrix(A %*% c(y.matrix[, 1: (J-1)]), ncol = J-1)

  # browser()
  # start1 <- Sys.time()
  pi.hat <- c(mu.hat[, 1], c(mu.hat[, 2:(J - 1)] - mu.hat[, 1:(J - 2)]))
  V <- Matrix(0, nrow = nSam * (J-1), ncol = nSam)
  V[cbind(c(1:(nSam * (J-1))), rep(1:nSam, J-1))] <- pi.hat
  V <- Diagonal(x = pi.hat) - Matrix::tcrossprod(V)
  V <- A %*% V %*% Matrix::t(A)
  VI <- chol2inv(chol(V))
  # end1 <- Sys.time()
  # end1 - start1
  # browser()

  res <- VI %*% c(y.tilde - mu.hat) - Zb.hat
  Sigma <- VI + ZVbZ
  SigmaI <- chol2inv(chol(Sigma))
  Y <- matrix(SigmaI %*% res, ncol = J - 1)
  KY <- tcrossprod(Y)

  return(KY)
}



