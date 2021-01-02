#' Function for visualizing Bayesian coresets
#'
#' This function visualizes coresets for binary response data. Note that these
#'   visualizations are limited to models with two covariates
#'   (not counting the intercept). Looking at these visualizations helps to
#'   understand the behavior of the two methods of obtaining coresets.
#'
#' @param coreset_results a list returned by either
#'   \code{\link{get_coreset_uniform}} or \code{\link{get_coreset_frankwolfe}}.
#' @param y the vector of binary response variables.
#' @param x the original, full feature matrix.
#' @param legend_true_label (optional) boolean, set to FALSE to omit plot's
#'   legend for the true class label. Default is TRUE.
#' @param legend_weights (optional) boolean, set to FALSE to omit the plot's
#'   legend for the coreset weights. Default is TRUE.
#' @param name_variables (optional) a vector of length 2 containing the
#'   desired axes labels. Default is c(expression(x[1]), expression(x[2])).
#' @param equation_line (optional) a vector of length 2 specifying the
#'   intercept and slope of the linear classification boundary to draw on the
#'   plot. The first element of the vector is the line's intercept, and the second
#'   element is the line's slope.
#' @param intercept_col (optional) an integer specifying the column index for
#'   the intercept column. Set to NA if model does not have an intercept.
#'   Default is the last column, i.e. x[, 3], which is the default of
#'   \code{\link{simulate_logit_data}}.
#'
#' @return a ggplot which visualizes the points selected for the coreset from
#'   the original data
#'   
#' @export
#' 
#' @import dplyr
#' @import ggplot2
#' @importFrom stats glm
#'
#' @examples
#' # Setup data and sample coreset
#' data <- simulate_logit_data()
#' coreset_results <- get_coreset_uniform(data$x, data$y)
#'
#' # Using all defaults
#' visualize_coresets(coreset_results, data$x, data$y)
#'
#' # Omit class label from legend
#' visualize_coresets(coreset_results, data$x, data$y,
#' legend_true_label = FALSE)
#'
#' # Omit coreset weights from legend
#' visualize_coresets(coreset_results, data$x, data$y, legend_weights = FALSE)
#'
#' # Customize axes labels
#' visualize_coresets(coreset_results, data$x, data$y,
#' name_variables = c("first variable", "second variable"))
#'
#' # Add glm equation results to plot
#' y_glm <- (data$y + 1) /2 # recode for glm call
#' glm_estimates <- glm(y_glm ~ data$x[, 1] + data$x[, 2], family = "binomial")
#' coeffs <- glm_estimates$coefficients
#' intercept <- - coeffs[1] / coeffs[3]
#' slope <- - coeffs[2] / coeffs[3]
#' visualize_coresets(coreset_results, data$x, data$y, equation_line = c(intercept, slope))
#'
#' # No intercept example
#' data <- simulate_logit_data(params = list(intercept = FALSE,
#' theta = c(3, 3)))
#' coreset_results <- get_coreset_uniform(data$x, data$y)
#' visualize_coresets(coreset_results, data$x, data$y, intercept_col = NA)

visualize_coresets <- function(coreset_results, x, y,
                               legend_true_label = TRUE,
                               legend_weights = TRUE,
                               name_variables = NA,
                               equation_line = NA,
                               intercept_col = 3) {

  # error check intercept_col
  if (!is.na(intercept_col)) {
    if (intercept_col < 1) {
      stop("Intercept column must be at least 1.")
    } else if (intercept_col > dim(x)[2]) {
      stop("Intercept column index is larger than the number of columns in x.")
    } else if (sum(x[, intercept_col]) != dim(x)[1]) {
      stop("Intercept column should be all 1's.")
    }
    
    # intercept column is not needed for plots
    x <- x[, -intercept_col]
  }
  
  # get coreset weights
  weights <- coreset_results$weights
  
  # format data
  mydata <- cbind(y, x)
  colnames(mydata) <- c("y", paste0("V", 1:2))
  mydata <- as_tibble(mydata)
  mydata <- mydata %>%
    mutate(weight = weights,
           selected = .data[["weight"]] > 0,
           true_label = case_when(.data[["y"]] == 1 ~ "success",
                                  .data[["y"]] != 1 ~ "failure"))
  
  # base plot
  plot <- ggplot(mydata) +
    aes(x = .data[["V1"]], y = .data[["V2"]]) +
    geom_point(aes(fill = .data[["true_label"]]), shape = 21, color = "blue",
               stroke = 0.25, alpha = 0.25, size = 1.5) +
    scale_fill_brewer(palette = "Set2") +
    geom_point(data = filter(mydata, .data[["selected"]] == TRUE),
               fill = "black",
               shape = 21,
               stroke = 0.01,
               aes(size = .data[["weight"]])) +
    labs(x = expression(x[1]),
         y = expression(x[2])) +
    theme_minimal() +
    theme(axis.title.x.bottom =  element_text(margin = margin(t = 8)),
          axis.title.y.left = element_text(margin = margin(r = 8))) +
    guides(fill = guide_legend(override.aes = list(alpha = 1),
                               title = "true label"))
  
  # set custom plot options
  if (legend_true_label == FALSE) {
    plot <- plot + guides(fill = FALSE)
  }
  
  if (legend_weights == FALSE) {
    plot <- plot + guides(size = FALSE)
  }
  
  if (!is.na(name_variables[1])) {
    plot <- plot + labs(x = name_variables[1],
                        y = name_variables[2])
  }
  
  if (!is.na(equation_line[1])) {
    plot <- plot +
      geom_abline(intercept = equation_line[1], slope = equation_line[2],
                  color = "gray47")
  }
  
  # return plot
  plot
}
