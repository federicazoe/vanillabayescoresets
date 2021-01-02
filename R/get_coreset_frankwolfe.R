#' Function for implementing Bayesian Coresets with the Hilbert Frank Wolfe 
#' method
#' 
#' This function implements Algorithm 2 and 3 described in the paper
#'   \href{https://arxiv.org/abs/1710.05053}{Automated Scalable Bayesian Inference 
#'   via Hilbert Coresets} for the case of logistic regression.
#'   It is a simplified version of the algorithm written in Python by the paper's
#'   first author, Trevor Campbell, which is available
#'   \href{https://github.com/trevorcampbell/bayesian-coresets}{here}. For
#'   obtaining the approximate posterior required in order to compute 
#'   (approximate) sensitivities and norms it uses the function `LaplaceApproximation()`
#'   from the package `LaplacesDemon`.
#'   
#'
#' @param x a (n x d) covariate matrix.
#' @param y a (n x 1) vector of {-1, 1} observed data.
#' @param m (optional) a scalar, indicating the maximum number of data points that 
#'   will be included in the returned Bayesian coreset. If not specified, this 
#'   value will be set to n/10 rounded down to the nearest integer.
#' @param num_projections a scalar indicating the dimension of finite-dimensional 
#'   approximation of log-likelihood to be used.
#' @param seed (optional) a value to set the seed. Since the uniform
#'   coreset method is not deterministic, given all other function inputs remain
#'   the same, this value can be set to ensure that the same coreset is returned.
#' @param verbose (optional) default is FALSE. If set to TRUE, status messages
#'   will be printed while the function runs.
#'
#' @return a list containing the following elements:
#'   \itemize{
#'     \item datapoints_selected - the indices for the data points selected to
#'           be part of the coreset.
#'     \item weights - the computed weights for each point in the coreset.
#'     \item xc - the feature matrix containing only the coreset data points.
#'     \item yc - the reponse vector containing only the responses for the
#'           observations in the coreset.
#'   }
#' @export
#' 
#' @importFrom LaplacesDemon LaplaceApproximation
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
#' # Using the defaults
#'
#' ## simulate binary data
#' data <- simulate_logit_data()
#' get_coreset_frankwolfe(data$x, data$y)
#'
#' # Customizations
#'
#' ## choose number of projections
#' get_coreset_frankwolfe(data$x, data$y, num_projections = 1000)
#'
#' ## change coreset size
#' get_coreset_frankwolfe(data$x, data$y, m = 100)
#'
#' ## set the seed
#' get_coreset_frankwolfe(data$x, data$y, seed = 5678)


get_coreset_frankwolfe <- function(x,
                                   y,
                                   m = NA,
                                   num_projections = 500,
                                   seed = 1234,
                                   verbose = FALSE) {
  
  # Outline:
  # 1) Obtain random projections of individual data points' log likelihoods:
  #    1.1) Sample feature points mu_j ~ Laplace approx. of posterior ; j = 1, ..., num_projections
  #    1.2) Sample gradient dimension d_j ~ Uniform(1, ..., d) ; j = 1, ..., num_projections
  #    1.3) Compute grad(loglik_i(mu_j))_d_j for each i = 1, ..., n and j = 1, ..., num_projections
  #    1.4) Construct Random Projection ell_hat_i using grad(loglik_i(mu_j))_d_j)
  # 2) Draw Frank-Wolfe coreset
  #    2.1) Compute approximate norms sigma_i using ell_hat_i
  #    2.2) Initialize:
  #         2.2.1) Compute ell_hat = sum(ell_hat_i)
  #         2.2.2) Find first_point = argmax (1 / sigma_i) * ell_hat^T * ell_hat_i

  set.seed(seed)
  
  # setting and error checking input parameters
  
  # get num obs and num features
  n <- dim(x)[1]
  d <- dim(x)[2]
  
  if (length(y) != n) {
    stop("x and y dimensions are incompatible.")
  }
  
  if (!is.na(m)) {
    if (m >= n) {
      stop("Coreset size m should be less than number of observations.")
    }
  } else {
    m <- as.integer(n / 10)
  }
  
  # Create a list with the structure of a Data object as required by 
  #   LaplaceApproximation()
  mydata <- format_data_laplace(x = x, y = y)
  
  # Initialize LaplaceApproximation() search
  initial_theta <- rep(0, d)
  
  print_verbose('Finding mode and covariance of Laplace posterior approximation',
                verbose)
  
  # Find approximate posterior mode and covariance through the Laplace approximation
  laplace_approx <- LaplaceApproximation(Model,
                                         parm = initial_theta,
                                         Data = mydata)
  
  approx_mode <- laplace_approx$Summary1[, "Mode"]
  approx_cov <- laplace_approx$Covar
  
  print_verbose(paste0('Sampling ', num_projections, ' sets of parameters for projections'),
                verbose)
  
  # Sample J thetas from N(approx_mode, approx_cov)
  thetas_for_projection <- rmvnorm(num_projections, 
                                   mean = approx_mode,
                                   sigma = approx_cov)
  
  print_verbose(paste0('Sampling ', num_projections, ' indices of parameters for projections'),
                verbose)
  
  # Sample gradient dimension d_j ~ Uniform(1,...,d), j = 1, ..., num_projections
  indices_for_projection <- sample(x = seq(1, d),
                                   size = num_projections,
                                   replace = TRUE)
  
  print_verbose("Computing gradients' projections, approximate sensitivities and approximate norms",
                verbose)
  
  # Compute gradients projections (hat_ell_n)
  ell_n_hat_matrix <- get_all_hat_ell_n(thetas_for_projection, 
                                        indices_for_projection, 
                                        y, 
                                        x)
  
  # Compute approximate sensitivities
  sigma_n_hat_vec <- sqrt(rowSums(ell_n_hat_matrix^2))
  
  # Compute sigma
  sigma_hat <- sum(sigma_n_hat_vec)
  
  # Compute vector ell_hat
  ell_hat <- colSums(ell_n_hat_matrix)
  
  # Compute (1/sigma_n)*ell_n for all n
  # (to avoid computing it at each iteration)
  ell_n_hat_divided_sigma_n <- ell_n_hat_matrix / sigma_n_hat_vec
  
  print_verbose('Initializing Frank-Wolfe Hilbert coreset',
                verbose)
  
  # Compute, for each n, (1 / sigma_n) * ell_hat^T * ell_hat_n
  norm_ells <-  ell_n_hat_divided_sigma_n %*% ell_hat
  
  # Select first point in coreset and set its weight
  first_pt <- which.max(norm_ells) 
  weights <- rep(0, n)
  weights[first_pt] <- sigma_hat / sigma_n_hat_vec[first_pt]
  
  print_verbose('Constructing Frank-Wolfe Hilbert coreset',
        verbose)
  
  # Running Frank-Wolfe Hilbert coreset contruction
  for (m in 2:m){
    
    ell_hat_w <- colSums(ell_n_hat_matrix * weights)
      
    ell_minus_ell_w <- ell_hat - ell_hat_w
    
    norm_ells <- ell_n_hat_divided_sigma_n %*% ell_minus_ell_w
    new_pt <- which.max(norm_ells) 
    
    gamma <- compute_gamma_frank_wolfe(new_pt,
                                       ell_n_hat_divided_sigma_n,
                                       ell_hat_w,
                                       ell_minus_ell_w,
                                       sigma_hat)
    
    dirac_new_pt <- rep(0, n)
    dirac_new_pt[new_pt] <- 1 
    
    weights <- ((1 - gamma) * weights) + 
      (gamma * (sigma_hat / sigma_n_hat_vec[new_pt]) * dirac_new_pt)
    
  }
  
  datapoints_selected <- seq(1, n)[(weights > 0)]
  coreset <- list(datapoints_selected = datapoints_selected, 
                  weights = weights,
                  xc = x[weights > 0, ],
                  yc = y[weights > 0])  
  
  coreset
}

#' Compute gamma as in Algorithm 2 of Campbell and Broderick (2019)
#' @param new_pt index of the observation selected in the current iteration
#' @param ell_n_hat_divided_sigma_n a n x d matrix corresponding to the object
#'  (1/sigma_f)*L_f, that features both in the left term of the numerator
#'  and in the right term of the denominator 
#' @param ell_hat_w a n x d matrix corresponding to the object L(w), which is 
#'   subtracted from(sigma/sigma_f)*L_f both in the left term of the 
#'   numerator and in the right term of the denominator
#' @param ell_minus_ell_w a n x d matrix corresponding to the object (L - L(w))
#'   featuring on the right term of the numerator 
#' @param sigma_hat a scalar corresponding to the object sigma, which
#'   premultiplies (1/sigma_f)*L_f on the right term of the numerator
#' @return a scalar 
#' @noRd 
compute_gamma_frank_wolfe <- function(new_pt,
                                      ell_n_hat_divided_sigma_n,
                                      ell_hat_w,
                                      ell_minus_ell_w,
                                      sigma_hat) {
  
  ell_f_w <- (sigma_hat * ell_n_hat_divided_sigma_n[new_pt, ]) - ell_hat_w
  
  gamma_num <- t(ell_f_w) %*% ell_minus_ell_w
  gamma_den <- sum(ell_f_w^2)
  
  c(gamma_num / gamma_den)
}
  


