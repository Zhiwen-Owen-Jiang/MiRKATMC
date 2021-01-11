#' The Microbiome Regression-based Kernel Association Test for multi-categorical
#' (nominal or ordinal) data
#'
#'  MiRKATMC is based on the baseline-category logit model for nominal data or
#'  the proportional odds model for ordinal data, which is able to deal with
#'  independent and clustered/longitudinal data.
#' @import stats
#' @importFrom Matrix Diagonal Matrix
#' @param formula the model formula with fixed effects only. The response should be
#' a factor.
#' @param data an optional data frame, ...
#' @param data.type either 'nominal' or 'ordinal'
#' @param Ks a kernel or a list of kernels
#' @param id an optional vector for subjects or clusters. Default is NULL for
#' independent response. It is required when slope is not null.
#' @param slope an optional vector for the random slope (e.g. time points).
#' Default si NULL.
#'
#' @return p value
#' @export
#'
#' @examples
MiRKATMC <- function(formula, data, data.type, Ks, id = NULL, slope = NULL){
  if (missing(data)) {
    data <- environment(formula)
  }

  if (any(is.na(data))){
    stop('There are missing values in the data. Please remove before proceeding')
  } else {
    data <- as.data.frame(data)
  }

  if (is.matrix(Ks)){
    Ks <- list(Ks)
  } else if (!is.list(Ks)){
    stop("Please enter either a single kernel matrix or a list of kernels for Ks.")
  }

  if (is.null(names(Ks)) & length(Ks) != 1){
    warning('No kernel names provided in the list of kernels')
  }

  if (data.type == 'nominal'){
    if (is.null(id) & is.null(slope)){
      KY <- inner.GLM.null(formula, data)
    } else if (id & is.null(slope) | id & slope){
      KY <- inner.GLM.random.null(formula, id, slope, data)
    } else {
      stop('id is required when random slope is provided.')
    }
  } else if (data.type == 'ordinal'){
    if (is.null(id) & is.null(slope)){
      KY <- inner.POM.null(formula, data)
    } else if (id & is.null(slope) | id & slope){
      KY <- inner.POM.random.null(formula, id, slope, data)
    } else {
      stop('id is required when random slope is provided.')
    }
  } else {
    stop('data.type should be either \"nominal\" or \"ordinal\".')
  }

  pvs <- mapply(DKAT, KY, Ks)
  if (length(pvs) > 1){
    pvs$p.value.combined <- harmonicmeanp::p.hmp(p = unlist(pvs), w = NULL, L = length(pvs),
                                  w.sum.tolerance = 1e-6, multilevel = F)[1]
  }
  return(pvs)
}




inner.GLM.null <- function(formula, data){
  # y <- eval(as.character(formula)[2], data)
  X <- model.matrix(formula, data)
  nSam <- dim(X)[1]

  # if (is.factor(y)){
  #   J <- length(levels(y))
  #   y <- factor(y, levels = c(1: J))
  #   y.matrix <- matrix(0, nrow = nSam, ncol = J)
  #   y.matrix[cbind(1:nSam, y)] <- 1
  #   # model.matrix( ~ 0 + y, test.data)
  # }

  fit <- mclogit::mblogit(formula, data)
  y <- fit$y
  if (fit$response.type == 'factor'){
    J <- nlevels(y)
    I <- diag(J)
    dimnames(I) <- list(levels(y), levels(y))
    y.matrix <- t(I[, y])
  } else if (fit$response.type == 'matrix'){
    y.matrix <- fit$y
    J <- ncol(y.matrix)
  }

  mu.hat <- matrix(fit$fitted.values, ncol = J, byrow = T)[, c(2:J, 1)]
  D <- get_D_GLM(mu.hat, J, nSam)
  V <- Matrix(0, nrow = nSam * (J-1), ncol = nSam)
  V[cbind(c(1:(nSam * (J-1))), rep(1:nSam, J-1))] <- c(mu.hat[, 1: (J-1)])
  V <- Diagonal(x = c(mu.hat[, 1: (J-1)])) - tcrossprod(V)
  VI <- chol2inv(chol(V))
  Y <- matrix(1/D * VI %*% c(y.matrix[, 2:J] - mu.hat[, 1:(J-1)]), ncol = J - 1)
  KY <- tcrossprod(Y)
    # Vb <- fit$VarCov$`1`
    # Z
    # ZVbZ <- kronecker(Vb, tcrossprod(Z))
    # rf <- c(matrix(unlist(fit$random.effects), ncol = J - 1, byrow = T))
    # Zb.hat <- kronecker(diag(1, J - 1), Z) %*% rf
    # res <- D * c(y.matrix[, 2:J] - mu.hat[, 1:(J-1)]) + Zb.hat
    # WI <- t(D * t(D * V))
    # Sigma <- WI + ZVbZ
    # SigmaI <- chol2inv(chol(Sigma))
    # Y <- matrix(SigmaI %*% res, ncol = J - 1)
  return(KY)
}




inner.GLM.random.null <- function(formula, id, slope, data){
  if (is.null(slope)){
    random <- ~ 1 | id
  } else {
    random <- ~ 1 + slope | id
  }
  fit <- mclogit::mblogit(formula, data, random)

  X <- model.matrix(formula, data)
  nSam <- dim(X)[1]
  y <- fit$y

  random <- setupRandomFormula(random)
  rt <- terms(random$formula)
  Z <- model.matrix(rt, data)
  Z <- get_Z(Z, data[, rt$groups], nSam)

  if (fit$response.type == 'factor'){
    J <- nlevels(y)
    I <- diag(J)
    dimnames(I) <- list(levels(y), levels(y))
    y.matrix <- t(I[, y])
  } else if (fit$response.type == 'matrix'){
    y.matrix <- fit$y
    J <- ncol(y.matrix)
  }

  if (is.null(slope)){
    Vb <- fit$VarCov$`1`
    ZVbZ <- kronecker(Vb, tcrossprod(Z))
  }
  else {
    Vb <- fit$VarCov$`1`[order(rep(seq(1, J-1), 2)), order(rep(seq(1, J-1), 2))]
    to.split <- kronecker(matrix(1:(J-1)^2, J-1, byrow = TRUE), matrix(1, 2, 2))
    Vb.list <- lapply(split(Vb, to.split), matrix, nr = 2)
    nclusters <- length(unique(data[, rt$groups]))
    Vb.list <- mapply(function(x){kronecker(diag(nclusters), x)}, Vb.list, SIMPLIFY = F)
    sqrt.len <- sqrt(length(Vb.list))
    Vb.m <- do.call(cbind, Vb.list[1: sqrt.len])
    Vb.list <- Vb.list[-(1: sqrt.len)]
    while(length(Vb.list) > 0){
      Vb.mi <- do.call(cbind, Vb.list[1: sqrt.len])
      Vb.m <- rbind(Vb.m, Vb.mi)
      Vb.list <- Vb.list[-(1: sqrt.len)]
    }
    Z.k <- kronecker(diag(J-1), Z)
    ZVbZ <- Z.k %*% Vb.m %*% t(Z.k)
  }
  rf <- c(matrix(unlist(fit$random.effects), ncol = J - 1, byrow = T))
  Zb.hat <- kronecker(diag(1, J - 1), Z) %*% rf

  mu.hat <- matrix(fit$fitted.values, ncol = J, byrow = T)[, c(2:J, 1)]
  D <- get_D_GLM(mu.hat, J, nSam)
  V <- Matrix(0, nrow = nSam * (J-1), ncol = nSam)
  V[cbind(c(1:(nSam * (J-1))), rep(1:nSam, J-1))] <- c(mu.hat[, 1: (J-1)])
  V <- Diagonal(x = c(mu.hat[, 1: (J-1)])) - tcrossprod(V)
  res <- D * c(y.matrix[, 2:J] - mu.hat[, 1:(J-1)]) + Zb.hat
  WI <- t(D * t(D * V))
  Sigma <- WI + ZVbZ
  SigmaI <- chol2inv(chol(Sigma))
  Y <- matrix(SigmaI %*% res, ncol = J - 1)
  KY <- tcrossprod(Y)
  return(KY)
}






inner.POM.null <- function(formula, data){
  X <- model.matrix(formula, data)
  nSam <- dim(X)[1]
  fit <- MASS::polr(formula, data, method = 'logistic')
  y <- eval(as.character(formula)[2], data)

  J <- nlevels(y)
  I <- diag(J)
  dimnames(I) <- list(levels(y), levels(y))
  y.matrix <- t(I[, y])

  pi.hat <- fit$fitted.values
  A <- matrix(1, nrow = J-1, ncol = J-1)
  A[upper.tri(A)] <- 0
  A <- kronecker(A, Diagonal(nSam, 1))

  mu.hat <- matrix(A %*% c(pi.hat[, 1: (J-1)]), ncol = J-1)
  y.tilde <- matrix(A %*% c(y.matrix[, 1: (J-1)]), ncol = J-1)
  KY <- tcrossprod(y.tilde - mu.hat)

  return(KY)
}


inner.POM.random.null <- function(formula, id, slope, data){
  X <- model.matrix(formula, data)
  nSam <- dim(X)[1]
  y <- eval(as.character(formula)[2], data)

  J <- nlevels(y)
  I <- diag(J)
  dimnames(I) <- list(levels(y), levels(y))
  y.matrix <- t(I[, y])

  if (is.null(slope)){
    rf <- sprintf('(1 | %s)',
                  as.character(substitute(id)))
  } else {
    rf <- sprintf('(1 + %s | %s)', as.character(substitute(slope)),
                  as.character(substitute(id)))
  }
  # rf <- sprintf('(%s)', as.character(random)[[2]])
  rf <- paste(c(".~.", rf), collapse="+")
  rf <- as.formula(rf)
  formula <- update(formula,rf)
  fit <- ordinal::clmm(formula, data, link = 'logit', threshold = 'flexible')

  if (is.null(slope)){
    random <- ~ 1 | id
    random <- setupRandomFormula(random)
    rt <- terms(random$formula)
    Z <- model.matrix(rt, data)
    Z <- get_Z(Z, data[, rt$groups], nSam)
  } else {
    Z <- t(fit$Zt)
  }

  Zb.hat <- rep(Z %*% fit$ranef, J - 1)
  ZVbZ <- kronecker(diag(J - 1), Z %*% fit$condVar %*% t(Z))

  X.new <- cbind(kronecker(diag(J - 1), rep(1, nSam)),
                 kronecker(rep(1, J - 1), -X[, -1]))
  eta.hat <- X.new %*% fit$coefficients - Zb.hat
  mu.hat <- matrix(exp(eta.hat)/(1 + exp(eta.hat)), ncol = J-1)
  V <- get_V_POM(mu.hat, J = J, nSam)
  D <- quickInverse(V, dimen = J, nSam)
  D.matrix <- get_matrix_W(D, J = J)

  A <- matrix(1, nrow = J-1, ncol = J-1)
  A[upper.tri(A)] <- 0
  A <- kronecker(A, Diagonal(nSam, 1))
  y.tilde <- matrix(A %*% c(y.matrix[, 1: (J-1)]), ncol = J-1)

  res <- D.matrix %*% c(y.tilde - mu.hat) - Zb.hat
  Sigma <- D.matrix + ZVbZ
  SigmaI <- chol2inv(chol(Sigma))
  Y <- matrix(SigmaI %*% res, ncol = J - 1)
  KY <- tcrossprod(Y)

  return(KY)
}



