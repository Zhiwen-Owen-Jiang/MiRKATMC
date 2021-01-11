get_D_GLM <- function(mu, J){
  ## reference level should be put in the last column
  d <- rep(0, (J-1)*nSam)
  n <- 1
  for (j in c(1:(J-1))){
    for (i in c(1:nSam)) {
      # num <- 1/(mu[i,j]*(1-mu[i,j]))
      # num <- 1/(mu[i,j])
      num <- (mu[i,j] + mu[i,J])/(mu[i,j]*mu[i,J])
      d[n] <- num
      n = n + 1
    }
  }
  return(d)
}



get_matrix_W <- function(W, J){
  W.matrix <- c()
  for (i in 1: (J - 1)){
    W.i <- diag(W[[1 + sum(0:(i-1))]])
    for (j in 2: (J - 1)){
      if (i == 1 | i == 2){
        W.i <- cbind(W.i, diag(W[[i + sum(0:(j-1))]]))
      } else {
        if (j <= i){
          W.i <- cbind(W.i, diag(W[[sum(0:(i-1)) + j]]))
        } else {
          W.i <- cbind(W.i, diag(W[[sum(0:(i-1)) + i + sum(i:(j - 1))]]))
        }
      }
    }
    W.matrix <- rbind(W.matrix, W.i)
  }
  return(W.matrix)
}


get_V_POM <- function(mu, J){
  ## reference level should be put in the last column
  vv <- rep(0, (J - 1) * nSam)
  vc <- rep(0, sum(1: (J - 2)) * nSam)
  n <- 1
  m <- 1

  for (j in c(1:(J-1))){
    for (k in c(1:(J-1))){
      for (i in c(1:nSam)) {
        if (j == k){
          num <- mu[i,j]*(1-mu[i,j])
          vv[n] <- num
          n = n + 1
        }
        else if (j > k){
          num <- mu[i,k]*(1 - mu[i,j])
          vc[m] <- num
          m = m + 1
        }
      }
    }
  }

  V <- list()
  k = 1
  for (i in c(1: (J-1))){
    for (j in c(1 : i)){
      if (i == j) {
        V[[paste0('w',i,j)]] = vv[((j - 1) * nSam + 1): (j * nSam)]
      }else{
        V[[paste0('w',j,i)]] = vc[((k - 1) * nSam + 1): (k * nSam)]
        k = k + 1
      }
    }
  }

  return(V)
}




blockdiag.prod <- function(A, B){
  ## A and B should be block diagonal matrix, stored in form of list containing
  ## only upper triangle (because here we are dealing with symmetric matrix)
  ## A should be a row vector, B should be square matrix
  out <- list()
  lenA <- length(A)
  # lenB <- length(B)
  for (i in 1:lenA){
    out[[paste0('w', i)]] <- rep(0, nSam)
    for (j in 1:lenA){
      if (i == 1 | i == 2){
        out[[paste0('w', i)]] <- out[[paste0('w', i)]] + A[[j]] * B[[i + sum(0:(j-1))]]
      } else{
        if (j <= i){
          out[[paste0('w', i)]] <- out[[paste0('w', i)]] + A[[j]] * B[[sum(0:(i-1)) + j]]
        } else{
          out[[paste0('w', i)]] <- out[[paste0('w', i)]] + A[[j]] * B[[sum(0:(i-1)) + i + sum(i:(j - 1))]]
        }
      }
    }
  }
  return(out)
}

matrixlist.prod <- function(A, B, type){
  ## product of two 'vectors' in form of list
  ## type 1 represents a row (1*n) times a column (n*1)
  ## type 2 represents a column (n*1) times a row (1*n)
  ## type 3 represents a number (1*1) times a row (1*n) or a column (n*1) times
  ## a number (1*1)
  lenA <- length(A)
  lenB <- length(B)
  if (lenA == lenB & type == 1){
    out <- rep(0, nSam)
    for (i in 1:lenA){
      out <- out + A[[i]] * B[[i]]
    }
    return(list(w1 = out))
  } else if (lenA == lenB & type == 2){
    out <- list()
    for (i in 1:lenA){
      for (j in 1:i){
        if (i == j){
          out[[paste0('w', i, j)]] <- A[[i]] * B[[j]]
        } else{
          out[[paste0('w', j, i)]] <- A[[i]] * B[[j]]
        }
      }
    }
    return(out)
  }


  if (lenA == 1 & type == 3){
    out <- list()
    for (i in 1:lenB){
      out[[paste0('w', i, lenB)]] <- A[[1]] * B[[i]]
    }
    return(out)
  }

  if (lenB == 1 & type == 3){
    out <- list()
    for (i in 1:lenA){
      out[[paste0('w', i, lenA)]] <- A[[i]] * B[[1]]
    }
    return(out)
  }
}

quickInverse <- function(M, dimen){
  ## suppose M is a square matrix
  length.M <- length(M)

  if (length.M != 3){
    Ai <- quickInverse.1011(M[c(1: sum(1:(dimen - 2)))], dimen - 1)
  }else {
    # Ai <- diag(1/diag(A))
    Ai <- list(w11 = 1/unlist(M[1]))
  }

  B <- M[c(sum(1: (dimen - 2)) + 1): c(sum(1: (dimen - 2)) + dimen - 2)]
  D <- unlist(M[length.M])

  # DI <- 1/(D - B * Ai * B)
  # BI <- -Ai * B * DI
  # AI <- Ai - BI * B * Ai
  # MI <- list(AI = AI, BI = BI, DI = DI)
  BAi <- blockdiag.prod(B, Ai)
  DI <- list()
  DI[[paste0('w', dimen, dimen)]] <- 1/(D - unlist(matrixlist.prod(BAi, B, type = 1)))
  BI <- lapply(matrixlist.prod(BAi, DI, type = 3), function(x) -1 * x)
  AI <- mapply('-', Ai, matrixlist.prod(BI, BAi, type = 2), SIMPLIFY = F)
  MI <- c(AI, BI, DI)
  return(MI)
}


setupRandomFormula <- function(formula){
  trms <- terms(formula)
  fo <- delete.response(trms)

  attributes(fo) <- NULL
  if(length(fo[[2]]) < 2 || as.character(fo[[2]][1])!="|")
    stop("missing '|' operator")

  groups <- fo
  fo[2] <- fo[[2]][2]
  groups[2] <- groups[[2]][3]

  list(
    formula=structure(fo,class="formula"),
    groups=all.vars(groups)
  )
}


get_Z <- function(Z, id){
  id <- factor(id)
  levels.id <- levels(id)
  nclusters <- length(levels.id)
  Z.m <- Matrix(0, nrow = nSam, ncol = dim(Z)[2] * nclusters)
  colnames(Z.m) <- rep(levels.id, dim(Z)[2])
  loc <- unlist(sapply(id, function(x){which(levels.id == x)}))
  Z.m[cbind(rep(1: nSam, dim(Z)[2]), rep(loc, dim(Z)[2]))] <- c(Z)
  return(Z.m)
}
