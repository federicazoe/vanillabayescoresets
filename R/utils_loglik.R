
#' Compute the log likelihood of theta given a sample from the logistic model 
#'   Y ~ Bernoulli(1/(1+exp(-X*theta)), where Y is coded either 1 or -1
#' @param y a n x 1 vector of binary observations (1 or -1)
#' @param x a n x d matrix of coefficients
#' @param theta a d x 1 vector of parameters at which the likelihood is
#'   evaluated
#' @return a scalar
#' @noRd
log_logit_likelihood <- function(x, y, theta){
  
  y_x_theta <- y*x %*% theta
  
  summands <- -y_x_theta
  # log(1-exp(n)) = (virtually) n, for n>=100
  less_100 <- (summands < 100) 
  summands[less_100] <- -log(1 + exp(summands[less_100]))
  summands[!less_100] <- -summands[!less_100]
  
  # We now obtain the loglikelihood by summing over summands
  sum(summands)
}

#' Compute the gradient of the log_likelihood of a single observation
#'   evaluated at the index_j element of the j-th parameter vector
#'    sampled from the approximate posterior
#' @param theta_for_projection_j a scalar at which the gradient is
#'   evaluated
#' @param index_j an index from 1 to d, that indicates which parameter index
#'   the evaluated theta_for_projection_j corresponds
#' @param y a single observation (1 or -1)
#' @param x a d x 1 vector containing the covariates of this observation
#' @return a scalar grad(loglik_i(mu_j))_{d_{index_1}} for this observation
#' @noRd

loglik_grad_datapt_elements <- function(theta_for_projection_j, 
                                        index_j, 
                                        y, 
                                        x){
  
  delta_loglik_obs <- - ((y * x) %*% theta_for_projection_j)
  lower_100 <- (delta_loglik_obs < 100)
  delta_loglik_obs[lower_100] <- exp(delta_loglik_obs[lower_100]) / 
    (1 + exp(delta_loglik_obs[lower_100]))
  delta_loglik_obs[!lower_100] <- 1
  delta_loglik_obs <- (y * x[index_j]) * c(delta_loglik_obs)
  
  delta_loglik_obs
  
}

#' Compute the gradient of the log_likelihood for each observation
#'   evaluated at the index_j-th element of the j-th parameter vector
#'    sampled from the approximate posterior
#' @param a n x 1 vector of binary observations (1 or -1)
#' @param x a n x d matrix of coefficients
#' @param theta_for_projection_j a scalar at which the gradient is
#'   evaluated
#' @param index_j an index from 1 to d, that indicates which parameter index
#'   the evaluated theta_for_projection_j corresponds
#' @return a n x 1 vector with delta(L_n(mu_j))_d_j for all n
#' @noRd

loglik_j_grad_elements <- function(theta_for_projection_j, 
                                   index_j, 
                                   y, 
                                   x){
  indexed_x <- x[, index_j]
  
  delta_loglik_j <- - ((y * x) %*% theta_for_projection_j)
  lower_100 <- (delta_loglik_j < 100)
  delta_loglik_j[lower_100] <- exp(delta_loglik_j[lower_100]) / 
    (1 + exp(delta_loglik_j[lower_100]))
  delta_loglik_j[!lower_100] <- 1
  delta_loglik_j <- (y * indexed_x) * c(delta_loglik_j)
  
  # Returns a n x 1 vector with delta(L_n(mu_j))_d_j for all n
  delta_loglik_j
  
}

#' Compute the gradient of the log_likelihood for each observation
#'   evaluated at each of the j-th index element of the j-th parameter vector
#'    sampled from the approximate posterior
#' @param a n x 1 vector of binary observations (1 or -1)
#' @param x a n x d matrix of covariates
#' @param thetas_for_projection a num_projections x 1 vector at which 
#'   the gradients are evaluated
#' @param indices a num_projections x 1 vector with indices from 1 to d, 
#'   that indicates for each element in the vector thetas_for_projection 
#'   which is the corresponding parameter index
#' @return n x num_projections matrix whose i-th row is  ell_hat_i, the gradient 
#'    of the likelihood's projection for observation i multiplied by the 
#'    square root of d/num_projections
#' @importFrom purrr pmap
#' @noRd

get_all_hat_ell_n <- function(thetas_for_projection, indices, y, x){
  
  # map all rows of thetas_for_projection, with the corresponding index,
  # to elements_grad_loglik_j() 
  
  thetas_projection_row_list <- as.list(as_tibble(t(thetas_for_projection),
                                                  .name_repair = "minimal"))
  indices_list <- as.list(indices)
  
  all_projections <- pmap(list(theta_for_projection_j = thetas_projection_row_list,
                               index_j = indices_list), 
                          loglik_j_grad_elements,
                          y = y,
                          x = x) # This produces a list whose elements are 
  # delta(L_n(mu_j))_d_j for all n
  
  # Convert the list to a matrix whose rows are
  # (grad(loglik_i(mu_1))_{d_{index_1}}, ..., 
  #   grad(loglik_i(mu_{num_projections}))_{d_{index_{num_projections}}})
  
  all_projections <- as.matrix(as_tibble(all_projections,
                                         .name_repair = "minimal"))
  
  # Multiply all_projections by sqrt(d/num_projections) to obtain a matrix whose
  # row is an ell_hat_i for i = 1, ..., n
  d <- dim(x)[2]
  num_projections <- dim(all_projections)[2]
  
  ell_n_hat_matrix <- sqrt(d/num_projections) * all_projections
  ell_n_hat_matrix
  
}
