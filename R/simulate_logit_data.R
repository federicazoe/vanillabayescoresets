#' Function for simulating binary data
#'
#' This function simulates binary data for a logistic regression model.
#' It currently supports drawing the feature matrix from a normal or bernoulli
#' distribution. These models are inspired by these two papers on coresets:
#' \itemize{
#'   \item \href{https://arxiv.org/abs/1605.06423}{Coresets for Bayesian
#'   Logistic Regression}.
#'   \item \href{https://arxiv.org/abs/1710.05053}{Automated Scalable Bayesian
#'   Inference via Hilbert Coresets}.
#' }
#'
#' @param model a string which provides the distribution according to which the
#'   feature matrix, X, will be drawn. Default value is "normal". Must be one of:
#'   \itemize{
#'     \item "normal".
#'     \item "bernoulli".
#'   }
#' @param params an optional list of model parameters that can be specified.
#'   Supported parameters include:
#'   \itemize{
#'     \item n - an integer specifying the number of observations to simulate.
#'     \item theta - a vector containing the true values of the logistic
#'           regression coefficients.
#'     \item intercept - a boolean for whether or not the model should include an
#'           intercept.
#'     \item sigma - a matrix to be used as the covariance matrix if the model is
#'           normal.
#'     \item px - a vector of probabilities to use for the draws if the model is
#'           bernoulli.
#'   }
#'   Default parameters for both models are n = 10000, theta = c(3, 3, 0),
#'   and intercept = TRUE.
#'   Default model is the normal model with the identity covariance matrix. If
#'   bernoulli model is specified, default probabilities are 0.5 for all features.
#'
#' @return a list with the response vector, feature matrix, and logistic
#'   regression coefficients. Details on the elements of the returned list:
#'   \itemize{
#'     \item y - vector of binary response variables (elements are 1 or -1).
#'     \item x - feature matrix.
#'     \item theta - vector of coefficients for the logistic regression model
#'       logit(p(y = 1)/(1 - p(y = 1))) = x * theta.
#'   }
#'   
#' @export
#'
#' @importFrom stats rbinom rnorm
#'
#' @examples
#' ## use all defaults (gives a normal model)
#' simulate_logit_data()
#'
#' ## use defaults for bernoulli model
#' simulate_logit_data(model = "bernoulli")
#'
#' # Possible customizations for both models
#'
#' ## set number of independent observations to 1000
#' simulate_logit_data(params = list(n = 1000))
#'
#' ## set true logistic regression coefficients to c(1, 1, 1)
#' simulate_logit_data(params = list(theta = c(1, 1, 1)))
#'
#' ## return data without intercept column
#' simulate_logit_data(params = list(intercept = FALSE))
#'
#' # Customizations for normal model
#'
#' ## change covariance matrix
#' simulate_logit_data(params = list(sigma = diag(10, 2)))
#'
#' # Customizations for bernoulli model
#'
#' ## change probability vector
#' simulate_logit_data(params = list(px = c(0.1, 0.2, 0.3)))

simulate_logit_data <- function(model = "normal", params = list()) {

  # extract general parameters
  if ("n" %in% names(params)) {
    n <- params$n
  } else {
    n <- 10000
  }

  if ("theta" %in% names(params)) {
    theta <- params$theta
  } else {
    theta <- c(3, 3, 0)
  }

  if ("intercept" %in% names(params)) {
    if (params$intercept == FALSE) {
      d <- length(theta)
    } else {
      d <- length(theta) - 1
    }
  } else {
    d <- length(theta) - 1
  }

  if (model != "bernoulli" & model != "normal") {
    # check that specified model is supported
    stop("model must be either normal or bernoulli.")
  }

  if (model == "normal") {
    # extract normal model parameters
    if ("sigma" %in% names(params)) {
      sigma <- params$sigma
    } else {
      sigma <- diag(1, d)
    }

    # simulate X
    a <- chol(sigma)
    z <- matrix(rnorm(n * d), nrow = n, ncol = d)
    x <- z %*% a

  } else if (model == "bernoulli") {
    # extract bernoulli parameters
    if ("px" %in% names(params)) {
      px <- params$px
      if (d != length(px)) {
        stop("The lengths of d and px are incompatible.")
      }

      if (!all(px <= 1 & px >= 0)) {
        stop("All probabilities must be between 0 and 1.")
      }
    } else {
      px <- 0.5
    }

    # simulate X
    x <- matrix(rbinom(n * d, 1, px), nrow = n, ncol = d, byrow = TRUE)

  }
  else {
    stop("Model must be normal or bernoulli")
  }

  # add intercept column if desired
  if (d < length(theta)) {
    x <- cbind(x, rep(1, n))
  }

  # simulate Y
  p <- 1 / (1 + exp(-x %*% theta))
  y <- 2 * rbinom(n, 1, prob = p) - 1

  # return as list
  list(y = y, x = x, theta = theta)
}
