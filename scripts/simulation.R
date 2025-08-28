library(tidyverse)
library(rms)
library(lme4)
source("data/random_data.r")

myfunc <- "10 + 5 * (1 + exp((.15 * (X1 - 60))))^(-1)"
input <- tribble(
  ~study, ~n, ~noise,     ~Y,         ~X1,
  "1", 300, 0, myfunc, "rnorm(n, 35, 1)",
  "2", 300, 0, myfunc, "rnorm(n, 40, 2)",
  "3", 300, 0, myfunc, "rnorm(n, 50, 3)",
  "4", 300, 0, myfunc, "rnorm(n, 60, 4)",
  "5", 300, 0, myfunc, "rnorm(n, 70, 5)",
  "6", 300, 0, myfunc, "rnorm(n, 80, 5)",
)

set.seed(2025)
ipd_data <- generate_from_table(input)

ipd_data <- ipd_data %>%
  rename(Age = X1) %>%
  mutate(Y = abs(Y), Study = as.factor(Study))
write_csv(ipd_data, "data/generated_ipd.csv")

# settings
n_sim <- 1000
knots_list <- c(3, 4, 5)
noise_sd <- 2

dd <- datadist(ipd_data)
options(datadist = "dd")

# restricted cubic spline formula
make_rcs <- function(x, var = "Age", knots = 3, dp = 1) {
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
    return(c(3))
  } else if (n < 200) {
    return(c(3, 4))
  } else {
    return(c(3, 4, 5))
  }
}

# valid knots across all studies
study_sizes <- ipd_data %>% 
  group_by(Study) %>% 
  summarise(n = n(), .groups = 'drop')
study_knots <- lapply(study_sizes$n, decide_knots)
study_valid_knots <- Reduce(intersect, study_knots)

# intersect user-defined list
valid_knots <- intersect(knots_list, study_valid_knots)
cat("Valid knots for all studies:", valid_knots, "\n")

# simulations for each knot setting
for (knots in valid_knots) {
  cat("Running simulation for knots =", knots, "\n")
  temp_data <- ipd_data
  temp_data$knots <- knots
  # initialize result storage
  pooled_results <- temp_data
  meta_results <- temp_data
  meta_re_results <- temp_data
  pooled_results$method <- "pooled"
  meta_results$method <- "meta"
  meta_re_results$method <- "meta_re"
  # for bias, variance, MSE calculations
  pooled_preds_matrix <- matrix(NA, nrow = nrow(ipd_data), ncol = n_sim)
  meta_preds_matrix <- matrix(NA, nrow = nrow(ipd_data), ncol = n_sim)
  meta_preds_re_matrix <- matrix(NA, nrow = nrow(ipd_data), ncol = n_sim)
  
  for (i in 1:n_sim) {
    # random shift to true Y (Y = Y + epsilon)
    set.seed(i)
    sim_data <- ipd_data %>% mutate(Y = Y + rnorm(n(), 0, noise_sd))
    
    # FIRST STAGE: fit individual study, get predictions on original data
    studies <- unique(sim_data$Study)
    pred_ori_data <- data.frame()
    
    for (s in studies) {
      study_data <- subset(sim_data, Study == s)
      pred_str <- make_rcs(study_data$Age, "Age", knots = knots)
      formula_obj <- as.formula(paste0("Y ~ ", pred_str))
      model <- ols(formula_obj, data = study_data)
      temp <- data.frame(
        Study = s,
        Age = study_data$Age,
        Y = predict(model, newdata = study_data)
      )
      pred_ori_data <- rbind(pred_ori_data, temp)
    }

    # SECOND STAGE: Method 1 - Pooled analysis
    combined_data <- pred_ori_data
    combined_data$Study <- 1  # Treat as single study
    combined_pred_str <- make_rcs(combined_data$Age, "Age", knots = knots)
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
    
    # SECOND STAGE: Method 3 - Meta-analysis Random Effects
    meta_pred_re <- predict(meta_model, newdata = ipd_data, re.form = ~(1 | Study))
    meta_preds_re_matrix[, i] <- meta_pred_re
    
    # add predictions to results
    pooled_results[[paste0("sim_", i)]] <- pooled_pred
    meta_results[[paste0("sim_", i)]] <- meta_pred
    meta_re_results[[paste0("sim_", i)]] <- meta_pred_re
  }
  # pooled_results$Study <- 1
  final_results <- rbind(pooled_results,meta_results, meta_re_results)
  # save in csv
  write_csv(final_results, paste0("./data/simul_knots", knots, ".csv"))
}

cat("Simulation completed!!!")
