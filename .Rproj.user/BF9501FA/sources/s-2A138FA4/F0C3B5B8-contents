library(devtools)
library(roxygen2)

#' Solves a linear system using Gauss-Seidel or Jacobi method
#'
#' @param A The coefficient matrix for solving x in b = Ax
#' @param b The outcome matrix for solving x in linear system b = Ax
#' @param v The correct outcome matrix for comparing relative error per iteration
#' @param iteration Optional Parameter for number of iterations desired for the algorithm; 10000 default
#' @param method Optional Parameter for method, choice of "gauss_seidel" or "jacobi"
#' @param no_core Optional Parameter for number of cores if using jacobi parallel
#'
#' @return A list of predicted x, relative error if given v, and time for each iteration
#'
#' @import doParallel
#' @import MASS
#' @export
#'
#' @examples A = matrix(0, nrow = 100, ncol=100)
#' diag(A) = rep(3, 100)
#' index = seq(1, 99)
#' A[cbind(index+1, index)] = rep(-1, 99)
#' A[cbind(index, index+1)] = rep(-1, 99)
#' v = matrix(rep(1,100), nrow=100)
#' b = A %*% v
#' solve_ols(A, b, v, method="gauss_seidel")
#' solve_ols(A, b, v, method="jacobi")
solve_ols = function(A, b, v=NULL, iteration=10000, method="gauss_seidel", no_core=1){
  if (method == "gauss_seidel"){
    gauss_seidel(A, b, v)
  } else if (method == "jacobi"){
    if (no_core == 1){
      jacobi(A, b, v)
    } else {
      jacobi_par(A, b, v, no_core)
    }
  } else{
    stop("Incorrect Method Name")
  }

}

gauss_seidel = function(A, b, v, iteration=10000){
  x = matrix(0, nrow=length(b))
  rel_err = c()
  time = c()
  timer0 = proc.time()[3]
  for (iter in 1:iteration) {
    for (i in 1:length(b)){
      temp = 0
      for (j in 1:length(b)){
        if (j != i) {
          temp = temp + A[i, j]*x[j]
        }
      }
      x[i] = 1/A[i,i]*(b[i]-temp)
    }
    if (!is.null(v)){
      rel_err[iter] = norm(x-v, type="2")/norm(v, type="2")
    }
    time[iter] = round(proc.time()[3] - timer0, digits=3)
  }
  list(x = x, rel_err = rel_err, time = time)
}

jacobi = function(A, b, v, iteration=10000){
  x = matrix(0, nrow=length(b))
  rel_err = c()
  time = c()
  timer0 = proc.time()[3]
  for (iter in 1:iteration) {
    temp_x = x
    for (i in 1:length(b)){
      temp = 0
      for (j in 1:length(b)){
        if (j != i) {
          temp = temp + A[i, j]*temp_x[j]
        }
      }
      x[i] = 1/A[i,i]*(b[i]-temp)
    }
    if (!is.null(v)){
      rel_err[iter] = norm(x-v, type="2")/norm(v, type="2")
    }
    time[iter] = round(proc.time()[3] - timer0, digits=3)
  }
  list(x = x, rel_err = rel_err, time = time)
}

jacobi_par = function(A, b, v, no_core = 1, iteration=10000){
  registerDoParallel(cores=7)
  x = matrix(0, nrow=length(b))
  rel_err = c()
  time = c()
  timer0 = proc.time()[3]
  for (iter in 1:iteration) {
    print(iter)
    temp_x = x
    x = foreach(i = 1:length(b), .combine=c) %dopar% {
      temp = 0
      for (j in 1:length(b)){
        if (j != i) {
          temp = temp + A[i, j]*temp_x[j]
        }
      }
      1/A[i,i]*(b[i]-temp)
    }
    if (!is.null(v)){
      rel_err[iter] = norm(x-v, type="2")/norm(v, type="2")
    }
    time[iter] = round(proc.time()[3] - timer0, digits=3)
  }
  list(x = x, rel_err = rel_err, time = time)
}

#' Implements algorithmic leveraging for linear regression using uniform and leverage score based subsampling of rows
#'
#' @param y the response variable matrix
#' @param X the predictor variable matrix
#' @param r the number of subsamples to choose from
#' @param plot if true, a plot of points selected will be shown
#' @param method either "unif" or "lev", determines wheter uniform or leverage score will be used
#'
#' @return the linear regression coefficient on the subsample
#'
#' @import MASS
#' @export
#'
#' @examples X = matrix(rt(500, 6), nrow=500)
#' ep = matrix(rnorm(500, 0, 1), nrow=500)
#' y = -X + ep
#' r = 100
#' algo_leverage(y, X, r, method="unif")
#' algo_leverage(y, X, r, method="lev")
algo_leverage = function(y, X, r, plot=F, method="unif"){
  if (method == "unif"){
    unif_subs(y, X, r, plot)
  } else if (method == "lev"){
    lev_subs(y, X, r, plot)
  } else{
    stop("Incorrect Method Name")
  }
}

unif_subs = function(y, X, r, plot=F) {
  prob = rep(1/(length(y)), length(y))
  samp = sample(seq(1, length(y)), r, prob = prob)
  data = data.frame(y = y[samp], X = X[samp])
  model = lm(y ~ X, data = data, weights = 1/prob[samp])
  if (plot) {
    plot(X, y, main="Uniform")
    points(X[samp], y[samp], col = "red", lwd=3)
    abline(a = model$coefficients[1], b = model$coefficients[2], col="red")
  }
  model$coefficients[2]
}

lev_subs = function(y, X, r, plot=F) {
  h = diag(X%*%solve(t(X)%*%X)%*%t(X))
  prob = h/sum(h)
  samp = sample(seq(1, length(y)), r, prob = prob)
  data = data.frame(y = y[samp], X = X[samp])
  model = lm(y ~ X, data = data, weights = 1/prob[samp])
  if (plot) {
    plot(X, y, main="Leverage")
    points(X[samp], y[samp], col = "red", lwd=3)
    abline(a = model$coefficients[1], b = model$coefficients[2], col="red")
  }
  model$coefficients[2]
}

#' Fits elastic net to data using coordinate descent algorithm
#'
#' @param X The predictor variable matrix
#' @param y The response variable matrix
#' @param alpha A value between 0 and 1 indicating weight toward L1 vs L2 shrinkage
#' @param lambda The weight parameter for determining significance of the shrinkage parameters
#' @param eps Stopping criteria
#'
#' @return The estimated coefficient vector beta
#'
#' @import MASS
#' @export
#'
#' @examples n = 100
#' alpha = 1
#' eps = matrix(rnorm(n, 0, 1), nrow=n)
#' X = mvrnorm(n, mu=rep(0, 20), Sigma = diag(rep(1, 20)))
#' X = scale(X)
#' beta = matrix(c(2,0,-2,0,1,0,-1,0,rep(0, 12)), nrow=20)
#' y = X%*%beta+eps
#' elnet_coord(X, y, alpha, 0.05)
elnet_coord = function(X, y, alpha, lambda, eps=0.000001){
  beta = matrix(rep(0,dim(X)[2]))
  nrow = dim(X)[2]
  n = length(y)
  for (itr in 1:1000){
    beta_old = beta
    for (j in 1:20){
      temp = 0
      for (i in 1:n){
        temp = temp + 1/n*X[i,j]*(y[i]-X[i,]%*%beta+X[i,j]*beta[j])
      }
      beta[j] = sign(temp)*max(abs(temp)-lambda*alpha, 0)/(1+2*lambda*(1-alpha))
    }
    if (sum(abs(beta-beta_old))<eps){
      break
    }
  }
  beta
}

