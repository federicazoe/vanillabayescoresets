
#' Create a list data formatted in the way required by LaplaceApproximation()
#' @param x a n x d matrix of coefficients
#' @param y a n x 1 vector of binary observations (1 or -1)
#' @return a list of data satisfying the minimum requirements to be considered 
#'   as a valid "Data" input of the function LaplaceApproximation()
#' @import LaplacesDemon
#' @importFrom stats rnorm
#' @noRd

format_data_laplace <- function(x, y) {
  
  n <- length(y)
  d <- dim(x)[2]
  
  mon_names <- "LP"
  parm_names <- as.parm.names(list(theta = rep(0, d)))
  pos_theta <- grep("theta", parm_names)
  
  PGF <- function(Data) {
    theta <- rnorm(Data$D)
    return(theta)
  }
  
  data_laplace_format <- list(D = d,
                              N = n,
                              PGF = PGF,
                              x = x,
                              mon.names = mon_names,
                              parm.names = parm_names,
                              pos.theta = pos_theta,
                              y = y)
  
  data_laplace_format
}



#' Generate a list model in the way required by LaplaceApproximation()
#' @param parm a d x 1 vector of parameters
#' @param Data a list of data satisfying the minimum requirements to be 
#'   considered a LaplaceDemon Data object
#' @return a list of model specifications satisfying the minimum requirements
#'    to be considered  as a valid "Model" input of the function 
#'    LaplaceApproximation()
#' @import LaplacesDemon
#' @noRd

Model <- function(parm, Data){
  
  ### Parameters
  theta <- parm
  d <- length(theta)
  ### Log-Prior
  log_prior <- - ((d / 2) * log(2 * pi)) - (sum(theta^2) / 2)
  ### Log-Likelihood
  mu <- tcrossprod(Data$x, t(theta))
  log_likelihood <- log_logit_likelihood(x = Data$x, y = Data$y, theta = theta)
  ### Log-Posterior
  log_posterior <- log_likelihood + log_prior
  Modelout <- list(LP = log_posterior, 
                   Dev = - (2 * log_likelihood), 
                   Monitor = log_posterior,
                   yhat = (rbinom(length(mu), 
                                  size = 1, 
                                  prob = (1 / (1 + exp(-mu)))) 
                           * 2 - 1),
                   parm = parm)
  return(Modelout)
}