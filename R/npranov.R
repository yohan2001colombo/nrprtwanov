#' Krss: Kernel Ridge Regression and ANOVA-based Analysis
#'
#' This package provides functions to perform kernel ridge regression, permutation-based significance testing,
#' and visualization of main and interaction effects.
#'
#' @docType package
#' @name npranov
NULL

#' Kernel Ridge Regression and ANOVA-based Analysis
#'
#' @param response Numeric vector (dependent variable).
#' @param A Factor variable (independent factor A).
#' @param B Factor variable (independent factor B).
#' @param n_perm Number of permutations for significance testing (default = 1000).
#' @return A dataframe with p-values and significance levels for factor A, factor B, and their interaction.
#' @export
#' @examples
#' data <- data.frame(
#'   Location = factor(rep(c("Olympia", "Ventura"), 6)),
#'   Tribe = factor(rep(c("Jedi", "Sith"), 6)),
#'   Midichlorians = c(10, 4, 12, 5, 15, 4, 15, 9, 15, 11, 18, 12)
#' )
#' Krss(data$Midichlorians, data$Location, data$Tribe)
krss <- function(response, A, B, n_perm = 1000) {
  # Function body
}

#' Generate an Interactive Boxplot
#'
#' @param data Data frame containing the factors and response variable.
#' @param factor Factor variable for grouping.
#' @param response Numeric response variable.
#' @param title Title for the plot.
#' @return An interactive boxplot using ggplot2 and plotly.
#' @export
#' @examples
#' main_effect_boxplot(data, "Location", "Midichlorians", "Main Effect of Location")
main_effect_boxplot <- function(data, factor, response, title) {
  # Function body
}

#' Generate an Interactive Interaction Plot
#'
#' @param data Data frame containing the factors and response variable.
#' @param factor1 First factor variable.
#' @param factor2 Second factor variable.
#' @param response Numeric response variable.
#' @return An interactive interaction plot using ggplot2 and plotly.
#' @export
#' @examples
#' interaction_plot(data, "Location", "Tribe", "Midichlorians")
interaction_plot <- function(data, factor1, factor2, response) {
  # Function body
}

Krss <- function(response, A, B, n_perm = 1000) {
  library(kernlab)
  library(Matrix)
  # For matrix operations

  # Gaussian kernel function
  gaussian_kernel <- function(x, y, sigma = 1) {
    exp(- (x - y)^2 / (2 * sigma^2))
  }

  # Create kernel matrices for two categorical factors
  kernel_matrix <- function(factor_levels, sigma = 1) {
    n <- length(factor_levels)
    K <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- gaussian_kernel(as.numeric(factor_levels[i]), as.numeric(factor_levels[j]), sigma)
      }
    }
    return(K)
  }

  # Function to optimize sigma using cross-validation
  optimize_sigma <- function(response, A, B, sigma_range = seq(0.1, 10, by = 0.1), k_folds = 5) {
    set.seed(123)  # For reproducibility
    fold_indices <- sample(1:k_folds, length(response), replace = TRUE)
    cv_errors <- numeric(length(sigma_range))

    for (s in 1:length(sigma_range)) {
      sigma <- sigma_range[s]
      K1 <- kernel_matrix(A, sigma)
      K2 <- kernel_matrix(B, sigma)
      K12 <- K1 * K2

      fold_errors <- numeric(k_folds)
      for (fold in 1:k_folds) {
        train_indices <- fold_indices != fold
        test_indices <- fold_indices == fold

        # Combine kernel matrices for training data
        K_train <- K1[train_indices, train_indices] + K2[train_indices, train_indices] + K12[train_indices, train_indices]

        # Add regularization (ridge penalty)
        lambda <- 1e-5  # Small regularization parameter
        K_train_reg <- K_train + lambda * diag(nrow(K_train))

        # Solve for kernel ridge regression coefficients
        alpha <- solve(K_train_reg, response[train_indices])

        # Construct test kernel matrix
        K1_test <- K1[test_indices, train_indices, drop = FALSE]
        K2_test <- K2[test_indices, train_indices, drop = FALSE]
        K12_test <- K12[test_indices, train_indices, drop = FALSE]
        K_test <- K1_test + K2_test + K12_test

        # Predict on test data
        predictions <- K_test %*% alpha
        fold_errors[fold] <- mean((response[test_indices] - predictions)^2)  # MSE
      }
      cv_errors[s] <- mean(fold_errors)  # Average MSE across folds
    }

    # Find the sigma with the lowest cross-validation error
    optimal_sigma <- sigma_range[which.min(cv_errors)]
    return(optimal_sigma)
  }

  # Optimize sigma
  optimal_sigma <- optimize_sigma(response, A, B)
  #cat("Optimal sigma:", optimal_sigma, "\n")

  # Compute kernel matrices with optimal sigma
  K1 <- kernel_matrix(A, optimal_sigma)  # Main effect of A
  K2 <- kernel_matrix(B, optimal_sigma)  # Main effect of B
  K12 <- K1 * K2                         # Interaction effect

  # Function to fit kernel ridge regression and compute RSS
  compute_rss <- function(K, response) {
    lambda <- 1e-5  # Small regularization parameter
    K_reg <- K + lambda * diag(nrow(K))
    alpha <- solve(K_reg, response)
    predictions <- K %*% alpha
    rss <- sum((response - predictions)^2)
    return(rss)
  }

  # Permutation test function
  permutation_test <- function(K_full, K_reduced, response, n_perm) {
    # Fit full and reduced models
    rss_full <- compute_rss(K_full, response)
    rss_reduced <- compute_rss(K_reduced, response)
    observed_statistic <- rss_reduced - rss_full  # Test statistic

    # Permutation loop
    perm_statistics <- numeric(n_perm)
    for (i in 1:n_perm) {
      perm_response <- sample(response)
      rss_full_perm <- compute_rss(K_full, perm_response)
      rss_reduced_perm <- compute_rss(K_reduced, perm_response)
      perm_statistics[i] <- rss_reduced_perm - rss_full_perm
    }

    # Compute p-value
    p_value <- mean(perm_statistics >= observed_statistic)
    return(p_value)
  }

  # Compute p-values
  p_value_main1 <- permutation_test(K1 + K2, K2, response, n_perm)  # Test for A
  p_value_main2 <- permutation_test(K1 + K2, K1, response, n_perm)  # Test for B
  p_value_interaction <- permutation_test(K1 + K2 + K12, K1 + K2, response, n_perm)  # Test for interaction

  # Create significance labels
  Significance <- ifelse(c(p_value_main1, p_value_main2, p_value_interaction) <= 0.05, "***", "")

  # Create a data frame similar to the ANOVA output
  cat("Response:", deparse(substitute(response)), "\n")
  anova_results <- data.frame(
    P_value = c(p_value_main1, p_value_main2, p_value_interaction),
    Significance = Significance,
    row.names = c(deparse(substitute(A)), deparse(substitute(B)), "Interaction")
  )

  print(anova_results)
  cat("\nSignificance Level = 0.05 '***'\n")
}


main_effect_boxplot <- function(data, factor, response, title) {
  library(ggplot2)
  ggplot(data, aes_string(x = factor, y = response)) +
    geom_boxplot(fill = "skyblue", color = "black") +
    theme_minimal() +
    labs(title = title, x = factor, y = deparse(substitute(response)))
}

interaction_plot <- function(data, factor1, factor2, response) {
  library(ggplot2)
  ggplot(data, aes_string(x = factor1, y = response, color = factor2, group = factor2)) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    theme_minimal() +
    labs(title = "Interaction Plot", x = factor1, y = deparse(substitute(response)), color = factor2)
}
#set.seed(123)
#soil_type <- factor(rep(c("Sandy", "Clay", "Loamy"), each = 30))  # Soil type (3 levels)
#watering_frequency <- factor(rep(c("Low", "Medium", "High"), times = 30))  # Watering frequency (3 levels)
#plant_height <- rnorm(90, mean = 20, sd = 5) + as.numeric(soil_type) * 2 + as.numeric(watering_frequency) * 3 + rnorm(90, sd = 2)

# Create a data frame
#data <- data.frame(soil_type, watering_frequency, plant_height)
#library(npranv)
#Krss(plant_height,soil_type,watering_frequency)

#interaction_plot(data=data,soil_type,watering_frequency,plant_height)
#main_effect_boxplot(data=data, soil_type,plant_height, "Main Effect of Soil Type")
#main_effect_boxplot(data, watering_frequency,plant_height, "Main Effect of Tribe")










