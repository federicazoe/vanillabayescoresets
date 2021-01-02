#' Function for getting uniform Bayesian coresets for logistic regression
#'
#' This function implements the algorithm described in the paper
#'   \href{https://arxiv.org/abs/1605.06423}{Coresets for Scalable Bayesian
#'   Logistic Regression}.
#'   It is a simplified version of the algorithm written in Python by the paper's
#'   first author, Jonathan H. Huggins, which is available
#'   \href{https://bitbucket.org/jhhuggins/lrcoresets/src/master/}{here}.
#'
#' @param x a feature matrix of dimension n x d
#' @param y a response vector of dimension n x 1
#' @param m (optional) a scalar, indicating the maximum number of data points that 
#'   will be included in the returned Bayesian coreset. If not specified, this 
#'   value will be set to n/10 rounded down to the nearest integer.
#' @param num_clusters (optional) the number of clusters for kmeans.
#'   If not specified, this value will be set to 4 for data with fewer than 10,000
#'   observations, and to n/2500 rounded down to the nearest integer for data
#'   with more than 10,000 observations.
#' @param r (optional) the radius of the ball representing the parameter space.
#'   Must be greater than 0. Per Huggins et. al., the default is to set to 4.
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
#' @importFrom stats kmeans dist rmultinom
#'
#' @examples
#' # Using the defaults
#'
#' ## simulate binary data
#' data <- simulate_logit_data()
#' get_coreset_uniform(data$x, data$y)
#'
#' # Customizations
#'
#' ## choose number of kmeans clusters
#' get_coreset_uniform(data$x, data$y, num_clusters = 5)
#'
#' ## change radius of parameter space
#' get_coreset_uniform(data$x, data$y, r = 3)
#'
#' ## change maximum coreset size
#' get_coreset_uniform(data$x, data$y, m = 500)
#'
#' ## set the seed
#' get_coreset_uniform(data$x, data$y, seed = 5678)

get_coreset_uniform <- function(x, y,
                                m = NA,
                                num_clusters = NA,
                                r = 4,
                                seed = 1234,
                                verbose = FALSE) {
  
  print_verbose("Running uniform coresets method", verbose)

  # allow for seed to be set for reproducible results
  set.seed(seed)

  # get num obs and num features
  n <- dim(x)[1]
  d <- dim(x)[2]

  # setting and error checking input parameters
  if (n < 100) {
    stop("Must be at least 100 observations.")
  }

  if (length(y) != n) {
    stop("x and y dimensions are incompatible.")
  }

  if (!is.na(num_clusters)) {
    if (num_clusters > n) {
      stop("Number of kmeans clusters num_clusters is larger than number of observations.")
    }
  } else {
    if (as.integer(n / 2500) < 4) {
      num_clusters <- 4
    } else {
      num_clusters <- as.integer(n / 2500)
    }
  }

  if (r <= 0) {
    stop("Radius r must be greater than 0.")
  }

  if (!is.na(m)) {
    if (m >= n) {
      stop("Coreset size m should be less than number of observations.")
    }
  } else {
    m <- as.integer(n / 10)
  }

  # get z matrix
  z <- matrix(NA, nrow = n, ncol = d + 1)
  z[, 1:d] <- y * x
  
  print_verbose("Starting kmeans clustering", verbose)

  # do kmeans clustering
  # cluster_sizes: number of data points assigned (i.e. closest) to each cluster
  # cluster_means: averages of the data assigned to each cluster
  # cluster_assignments: cluster assignments for each data point
  k_means_results <- kmeans(z[, 1:d], num_clusters)
  cluster_sizes <- k_means_results$size
  cluster_means <- k_means_results$centers
  cluster_assignments <- k_means_results$cluster
  z[, d + 1] <- cluster_assignments
  
  print_verbose("Kmeans clustering completed", verbose)
  
  print_verbose("Calculating sensitivities", verbose)

  # helper function to apply to z matrix in next step
  calc_sensitivity <- function(obs) {
    # get the cluster assignment, mean, and size
    k <- obs[length(obs)] # last col is cluster assignment
    obs <- obs[1:(length(obs) - 1)]
    true_mean <- cluster_means[k, ]
    nk <- cluster_sizes[k]
    log_sizes <- log(cluster_sizes)
    # G_k^(-n) is of size nk - 1 since it excludes this data pt
    log_sizes[k] <- log(nk - 1)

    # update the mean without the data row
    cluster_means[k, ] <- cluster_means[k, ] * nk / (nk - 1)
    cluster_means[k, ] <- cluster_means[k, ] - 1 / (nk - 1) * obs

    # calculate denominator
    dists <- apply(cluster_means, 1, function(x) dist(rbind(x, obs)))
    exp_arg <- log_sizes - r * dists
    denominator <- 1 + sum(exp(exp_arg))

    # restore cluster mean
    cluster_means[k, ] <- true_mean

    # calculate and return sensitivity
    sensitivity <- n / denominator
    sensitivity
  }

  # calculate sensitivity upper bounds and mean
  sensitivities <- apply(z, 1, calc_sensitivity)
  mean_sensitivity <- (1 / n) * sum(sensitivities)
  
  print_verbose("Sampling data for coreset", verbose)

  # calculate importance weights
  p <- sensitivities / (n * mean_sensitivity)

  # sample data for coreset
  samp <- rmultinom(1, m, p)
  
  print_verbose("Calculating weights", verbose)

  # calculate coreset weights and get indices of selected points
  weights <- (1 / m) * samp / p
  datapoints_selected <- seq(1, n)[(weights > 0)]

  # return data with nonzero weights
  coreset <- list(datapoints_selected = datapoints_selected,
                  weights = weights,
                  xc = x[weights > 0, ],
                  yc = y[weights > 0])
  print_verbose("Complete! Returning coreset", verbose)
  coreset
}
