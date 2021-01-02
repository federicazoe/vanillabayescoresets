

# Correct computation in toy example --------------------------------------

test_that(
  "Loglikelihood is computed correctly in a toy example",
  {
    y <- c(1, -1)
    x <- matrix(c(0.1, 0.5, -0.2, -0.3), nrow = 2, byrow = TRUE)
    theta <- c(3, 3)
    result_expected <- - 
      (log(1 + exp(- (y[1] * x[1, 1] * theta[1] + y[1] * x[1, 2] * theta[2]))) +
         log(1 + exp(- (y[2] * x[2, 1] * theta[1] + y[2] * x[2, 2] * theta[2]))))
    expect_equal(log_logit_likelihood(x, y, theta), result_expected)
  }
)



# No overflow for very large -Y_X_theta -----------------------------------------

test_that(
  "Computation of the loglikelihood does not incur into problems due to overflow
  when (-y * x * theta) is very large causing 
  log(1 + exp(1 + exp(-y * x * theta)) = Inf  instead of  
  log(1 + exp(-y * x * theta)) = (-y * x * theta), which is virtually  
  true when (-y * x * theta)  is very large ",
  {
    y <- c(1, -1)
    x <- matrix(c(-1000, 0.5, -0.2, -0.3), nrow = 2, byrow = TRUE)
    theta <- c(3, 3)
    result_expected <- - 
      (-(y[1] * x[1, 1] * theta[1] + y[1] * x[1, 2] * theta[2]) +
         log(1 + exp(- (y[2] * x[2, 1] * theta[1] + y[2] * x[2, 2] * theta[2])))) 
    expect_equal(log_logit_likelihood(x, y, theta), result_expected)
  }
)

