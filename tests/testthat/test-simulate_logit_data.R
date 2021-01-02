
# Correct number of observations ------------------------------------------

test_that(
  "default number of observations is 10000", 
  {
    y <- simulate_logit_data()$y
    expect_length(y, 10000)
  }
)

test_that(
  "correct number of observations for y", 
  {
    y <- simulate_logit_data(params = list(n = 10))$y
    expect_length(y, 10)
  }
)

test_that(
  "correct number of observations for x", 
  {
    x <- simulate_logit_data(params = list(n = 10))$x
    expect_equal(dim(x)[1], 10)
  }
)

# Correct number of features ----------------------------------------------

test_that(
  "correct number of features for x w/ intercept", 
  {
    x <- simulate_logit_data(params = list(theta = rep(0.5, 5)))$x
    expect_equal(dim(x)[2], 5)
  }
)

test_that(
  "correct number of features for x w/o intercept", 
  {
    x <- simulate_logit_data(params = list(theta = rep(0.5, 5), 
                                           intercept = FALSE))$x
    expect_equal(dim(x)[2], 5)
  }
)

# Correct feature properties ----------------------------------------------

test_that(
  "intercept column is all 1's", 
  {
    x <- simulate_logit_data(params = list(theta = rep(0.5, 5)))$x
    expect_setequal(x[, ncol(x)], 1)
  }
)

test_that(
  "all x's are 0 or 1 if model is bernoulli", 
  {
    x <- simulate_logit_data(model = "bernoulli")$x
    expect_true(all(as.vector(x) %in% c(0, 1)))
  }
)

# Correct response properties ---------------------------------------------
test_that(
  "all y's are -1 or 1", 
  {
    y <- simulate_logit_data()$y
    expect_true(all(as.vector(y) %in% c(-1, 1)))
  }
)

