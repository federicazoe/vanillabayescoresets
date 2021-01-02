
# Correct coreset size --------------------------------------

test_that(
  "size of obtained coreset matches required size", 
  {
    data <- simulate_logit_data(param = list(n = 1000))
    datapoints_selected <- get_coreset_frankwolfe(x = data$x,
                                               y = data$y,
                                               m = 100)$datapoints_selected
    maximum_size <- (length(datapoints_selected) <= 100)
    expect_true(maximum_size)
  }
)