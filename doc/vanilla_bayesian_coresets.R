## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(vanillabayescoresets)

## -----------------------------------------------------------------------------
data <- simulate_logit_data()
names(data)

## -----------------------------------------------------------------------------
length(data$y)
dim(data$x)[1]
dim(data$x)[2]
length(data$theta)
head(data$y)

## -----------------------------------------------------------------------------
data <- simulate_logit_data(params = list(n = 5000))
length(data$y)
dim(data$x)[1]

## -----------------------------------------------------------------------------
data <- simulate_logit_data(params = list(n = 5000, theta = c(1, 2, 3, 4)))
length(data$theta)
dim(data$x)[2]

## -----------------------------------------------------------------------------
head(data$x)

## -----------------------------------------------------------------------------
dependent_sigma <- matrix(c(1, 0.1, 0.2, 1), byrow = TRUE, ncol = 2, nrow = 2)
dependent_sigma
data <- simulate_logit_data(params = list(sigma = dependent_sigma))

## -----------------------------------------------------------------------------
data <- simulate_logit_data(model = "bernoulli")
head(data$x)

## -----------------------------------------------------------------------------
data <- simulate_logit_data(model = "bernoulli", params = list(px = c(0.2, 1)))
head(data$x)

## -----------------------------------------------------------------------------
data <- simulate_logit_data(model = "bernoulli", 
                            params = list(theta = c(3, 3), 
                                          px = c(0.2, 1), 
                                          intercept = FALSE))
head(data$x)

## -----------------------------------------------------------------------------
data <- simulate_logit_data()
coreset <- get_coreset_uniform(data$x, data$y)
length(coreset$datapoints_selected)

## -----------------------------------------------------------------------------
coreset <- get_coreset_uniform(data$x, data$y, m = 500)
length(coreset$datapoints_selected)

## -----------------------------------------------------------------------------
coreset <- get_coreset_uniform(data$x, data$y, num_clusters = 5)

## ---- message = FALSE---------------------------------------------------------
data <- simulate_logit_data()
dplyr::glimpse(data$x)
dplyr::glimpse(data$y)
coreset <- get_coreset_frankwolfe(x = data$x, y = data$y)

## ---- message = FALSE---------------------------------------------------------
coreset <- get_coreset_frankwolfe(x = data$x,
                                  y = data$y,
                                  m = as.integer(length(data$y) / 10),
                                  num_projections = 500,
                                  seed = 1234,
                                  verbose = FALSE)

## ---- fig.align="center"------------------------------------------------------
# Plotting uniform coresets
coreset_uniform <- get_coreset_uniform(data$x, data$y, num_clusters = 5)
visualize_coresets(coreset_uniform, data$x, data$y, intercept_col = 3)

## ---- , fig.align="center"----------------------------------------------------
# Plotting Hilbert Frank-Wolfe coresets
coreset_frankwolfe <- get_coreset_frankwolfe(data$x, data$y)
visualize_coresets(coreset_frankwolfe, data$x, data$y, intercept_col = 3)

## ---- fig.align="center"------------------------------------------------------
# Making optional changes
y_glm <- (data$y + 1) /2 # recode for glm call
glm_estimates <- glm(y_glm ~ data$x[, 1] + data$x[, 2], family = "binomial")
coeffs <- glm_estimates$coefficients
intercept <- - coeffs[1] / coeffs[3]
slope <- - coeffs[2] / coeffs[3]

visualize_coresets(coreset_frankwolfe, 
                   data$x, 
                   data$y, 
                   legend_true_label = FALSE,
                   legend_weights = FALSE,
                   name_variables = c("Variable 1", "Variable 2"),
                   equation_line = c(intercept, slope))

