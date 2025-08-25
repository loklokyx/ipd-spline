library(tidyverse)
library(rms)
library(lme4)
source("data/random_data.r")

myfunc <- "10 + 5 * (1 + exp((.15 * (X1 - 60))))^(-1)"
input <- tribble(
  ~study, ~n, ~noise,     ~Y,         ~X1,
  "1", 100, 0, myfunc, "rnorm(n, 35, 1)",
  "2", 100, 0, myfunc, "rnorm(n, 40, 2)",
  "3", 100, 0, myfunc, "rnorm(n, 50, 3)",
  "4", 100, 0, myfunc, "rnorm(n, 60, 4)",
  "5", 100, 0, myfunc, "rnorm(n, 70, 5)",
  "6", 100, 0, myfunc, "rnorm(n, 80, 5)",
)

set.seed(2025)
ipd_data <- generate_from_table(input)

ipd_data <- ipd_data %>%
  rename(X = X1) %>%
  mutate(Y = abs(Y), Study = as.factor(Study))

# settings
n_sim <- 1
knots_list <- c(3, 4, 5)
noise_sd <- 2

dd <- datadist(ipd_data)
options(datadist = "dd")

# restricted cubic spline formula
make_rcs <- function(x, var = "X", knots = 3, dp = 1) {
  quan <- list(
    '3' = c(0.10, 0.50, 0.90),
    '4' = c(0.05, 0.35, 0.65, 0.95),
    '5' = c(0.05, 0.275, 0.50, 0.725, 0.95)
  )
  if (as.character(knots) %in% names(quan)) {
    knot_values <- quantile(x, probs = quan[[as.character(knots)]])
    sprintf("rcs(%s, c(%s))", var, paste(round(knot_values, dp), collapse = ", "))
  } else {
    sprintf("rcs(%s)", var)
  }
}

decide_knots <- function(n) {
  if (n < 50) {
    return(3)
  } else if (n < 200) {
    return(c(3, 4))
  } else {
    return(c(3, 4, 5))
  }
}

# simulations for each knot setting
for (knots in knots_list) {
  cat("Running simulation for knots =", knots, "\n")
  
  # initialize result storage
  pooled_results <- ipd_data
  meta_results <- ipd_data
  # for bias, variance, MSE calculations
  pooled_preds_matrix <- matrix(NA, nrow = nrow(ipd_data), ncol = n_sim)
  meta_preds_matrix <- matrix(NA, nrow = nrow(ipd_data), ncol = n_sim)
  
  for (i in 1:n_sim) {
    # random shift to true Y (Y = Y + epsilon)
    set.seed(i)
    sim_data <- ipd_data %>% mutate(Y = Y + rnorm(n(), 0, noise_sd))
    
    # FIRST STAGE: fit individual study, get predictions on original data
    studies <- unique(sim_data$Study)
    pred_ori_data <- data.frame()
    
    for (s in studies) {
      study_data <- subset(sim_data, Study == s)
      pred_str <- make_rcs(study_data$X, "X", knots = knots)
      formula_obj <- as.formula(paste0("Y ~ ", pred_str))
      model <- ols(formula_obj, data = study_data)
      temp <- data.frame(
        Study = s,
        X = study_data$X,
        Y = predict(model, newdata = study_data)
      )
      
      pred_ori_data <- rbind(pred_ori_data, temp)
    }

    # SECOND STAGE: Method 1 - Pooled analysis
    combined_data <- pred_ori_data
    combined_data$Study <- 1  # Treat as single study
    combined_pred_str <- make_rcs(combined_data$X, "X", knots = knots)
    pooled_formula <- as.formula(paste0("Y ~ ", combined_pred_str))
    pooled_model <- ols(pooled_formula, data = combined_data)
    # predict on original data
    pooled_pred <- predict(pooled_model, newdata = ipd_data)
    pooled_preds_matrix[, i] <- pooled_pred
    
    # SECOND STAGE: Method 2 - Meta-analysis
    meta_data <- pred_ori_data
    meta_formula <- as.formula(paste0("Y ~ ", combined_pred_str, " + (1 | Study)"))
    meta_model <- lmer(meta_formula, data = meta_data)
    meta_pred <- predict(meta_model, newdata = ipd_data, re.form = NA)
    meta_preds_matrix[, i] <- meta_pred
    
    # add predictions to results
    pooled_results[[paste0("sim_", i)]] <- pooled_pred
    meta_results[[paste0("sim_", i)]] <- meta_pred
  }

  # calculate bias, variance, MSE
  pooled_bias <- rowMeans(pooled_preds_matrix) - ipd_data$Y
  pooled_variance <- apply(pooled_preds_matrix, 1, var)
  pooled_mse <- rowMeans((pooled_preds_matrix - ipd_data$Y)^2)
  
  meta_bias <- rowMeans(meta_preds_matrix) - ipd_data$Y
  meta_variance <- apply(meta_preds_matrix, 1, var)
  meta_mse <- rowMeans((meta_preds_matrix - ipd_data$Y)^2)
  
  # Add bias, variance, MSE to results
  pooled_results$bias <- pooled_bias
  pooled_results$variance <- pooled_variance
  pooled_results$mse <- pooled_mse
  pooled_results$Study <- 1
  
  meta_results$bias <- meta_bias
  meta_results$variance <- meta_variance
  meta_results$mse <- meta_mse
  
  # save in csv
  write_csv(pooled_results, paste0("./data/pooled_knots", knots, ".csv"))
  write_csv(meta_results, paste0("./data/meta_knots", knots, ".csv"))  
}

cat("Simulation completed!!!")