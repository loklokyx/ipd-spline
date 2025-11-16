library(gtools)
library(tidyverse)
library(rms)
library(lme4)
source("data/random_data.r")
source("scripts/helpers.R")

# Data generation
myfunc <- "10 + 5 * (1 + exp((.15 * (X1 - 60))))^(-1)"

# Configuration - keep all settings together
CONFIG <- list(
  data = NA,
  # n_sim = 1000,
  # knots_list = c(3, 4),
  n_sim = 1,
  knots_list = c(3),
  stage2_knots_list = c(5),
  noise_sd = 1,
  weight_fun = function(se) 1/se^2,  # inverse variance weighting
  overlap_conditions = c("simulated_no", "simulated_some", "simulated_many"),
  setting_dir = "./data/extrapolate/beta/"
)

# Get valid knot combinations
get_valid_combinations <- function(knots_list) {
  study_ids <- unique(CONFIG$data$Study)
  # Get study sizes
  study_max_n <- sapply(study_ids, function(s) sum(CONFIG$data$Study == s))

  valid_knots_per_study <- lapply(study_max_n, function(n) {
    max_knot <- get_max_knots(n)
    all_knots <- 3:max_knot # assume min knots = 3
    intersect(all_knots, knots_list)
  })
  names(valid_knots_per_study) <- study_ids
  
  expand.grid(valid_knots_per_study, KEEP.OUT.ATTRS = FALSE)
}

# Stage 1: Individual study predictions
run_stage1 <- function(sim_data, study_ids, stage1_knots_vec) {
  pred_stage1 <- data.frame()
  stage1_probs <- list()
  
  for (s_idx in seq_along(study_ids)) {
    s <- study_ids[s_idx]
    study_data <- subset(sim_data, Study == s)
    knots_s <- stage1_knots_vec[s_idx]
    
    probs <- make_quantiles(study_data$Age, knots_s)
    stage1_probs[[as.character(s)]] <- probs
    
    pred_str <- make_rcs_formula("Age", probs)
    formula_obj <- as.formula(paste0("Y ~ ", pred_str))
    model <- lm(formula_obj, data = study_data)
    
    # Predict on all data points (for pointwise combining)
    newdata <- data.frame(Age = sim_data$Age)
    pred_se <- predict(model, newdata = newdata, se.fit = TRUE)
    pred_ci <- predict(model, newdata = newdata, interval = "confidence")
    
    temp <- data.frame(
      Study = s,
      Age = sim_data$Age,
      Y = pred_se$fit,
      SE = pred_se$se.fit,
      lower = pred_ci[, "lwr"],
      upper = pred_ci[, "upr"]
    )
    pred_stage1 <- rbind(pred_stage1, temp)
  }
  
  list(data = pred_stage1, probs = stage1_probs)
}

# Pointwise combining function with weights
combine_pointwise <- function(pred_data, weight_fun = NULL) {
  d <- pred_data
  
  if (!is.null(weight_fun)) {
    d <- d %>% mutate(w = weight_fun(SE))
    
    d_summary <- d %>%
      group_by(Age) %>%
      summarise(
        Y_pred = sum(Y * w) / sum(w),
        var_weighted = 1 / sum(w),
        lower = Y_pred - 1.96 * sqrt(var_weighted),
        upper = Y_pred + 1.96 * sqrt(var_weighted),
        .groups = "drop"
      )
  } else {
    d_summary <- d %>%
      group_by(Age) %>%
      summarise(
        Y_pred = mean(Y),
        se = sd(Y) / sqrt(n()),
        lower = Y_pred - 1.96 * se,
        upper = Y_pred + 1.96 * se,
        .groups = "drop"
      )
  }
  
  d_summary %>% 
    rename(Y = Y_pred) %>% 
    slice(match(unique(d$Age), Age)) # reorder by original Age
}

# Stage 2: Pooled method
run_pooled <- function(pred_stage1, weight_fun) {
  combine_pointwise(pred_stage1, weight_fun)
}

# Stage 2: Meta-analysis fixed effects with weights
run_meta_fixed <- function(pred_stage1, stage2_knots, weight_fun) {
  d <- pred_stage1 %>% mutate(weight = weight_fun(SE))
  
  probs2 <- make_quantiles(d$Age, stage2_knots)
  pred_str <- make_rcs_formula("Age", probs2)
  meta_formula <- as.formula(paste0("Y ~ ", pred_str, " + (1 | Study)"))
  meta_model <- lmer(meta_formula, weights = weight, data = d)
  
  # Predict on original data
  newdata <- pred_stage1 %>% select(Age, Study)
  newdata$Y <- predict(meta_model, newdata = newdata, re.form = NA)
  
  # Approximate SE using residual sigma and weights
  sigma_hat <- summary(meta_model)$sigma
  
  # Get weights for prediction points
  pred_with_weights <- d %>%
    group_by(Age) %>%
    summarise(mean_weight = mean(weight), .groups = "drop")
  
  newdata <- newdata %>%
    left_join(pred_with_weights, by = "Age") %>%
    mutate(SE = sigma_hat / sqrt(mean_weight))
  
  # Combine pointwise
  list(results = combine_pointwise(newdata, weight_fun), stage2_probs = probs2)
}

# Single simulation run
run_single_simulation <- function(sim_i, study_ids, stage1_knots_vec, stage2_knots) {
  set.seed(sim_i)
  sim_data <- CONFIG$data %>% 
    mutate(Y = Y + rnorm(n(), 0, CONFIG$noise_sd)) %>%
    group_by(Study) %>%
    mutate(Age = make_unique_by_precision(Age, digits = 5)) %>%
    ungroup()
  
  stage1_results <- run_stage1(sim_data, study_ids, stage1_knots_vec)
  # Stage 2: Methods
  pooled_result <- run_pooled(stage1_results$data, CONFIG$weight_fun)
  meta_result <- run_meta_fixed(stage1_results$data, CONFIG$stage2_knots, CONFIG$weight_fun)
  stage2_results <- list(pooled = pooled_result, meta = meta_result$results)
  list(
    predictions = stage2_results,
    stage1_probs = if (sim_i == 1) stage1_results$probs else NULL,
    stage2_probs = if (sim_i == 1) meta_result$stage2_probs else NULL
  )
}

# Format and save results
save_simulation_results <- function(all_results, filename) {
  # Initialize method results storage
  method_names <- names(all_results[[1]]$predictions)
  # Age_sim <- all_results[[1]]$predictions[[1]]$Age
  mse_list <- list()   # <-- initialize mse_list!
  
  # Compute weighted MSE for each method
  Y_true <- CONFIG$data$Y
  # Save simulated predictions
  row_based_df <- data.frame()
  for (method in method_names) {
    sim_matrix <- do.call(cbind, lapply(
      all_results, function(res) res$predictions[[method]]['Y']
    ))
    
    temp_df <- data.frame(method = method, sim_matrix)
    row_based_df <- rbind(row_based_df, temp_df)

    row_mses <- rowMeans((sim_matrix - Y_true)^2)
    overall_mse <- mean(row_mses)  # average MSE for this setting

    mse_list[[paste0("MSE_", method)]] <- overall_mse
  }
  
  # Rename sim columns
  colnames(row_based_df) <- c("method", paste0("sim_", seq_len(CONFIG$n_sim)))
  # final_df <- cbind(CONFIG$data, Age_sim, row_based_df)
  final_df <- cbind(CONFIG$data, row_based_df)
  
  write_csv(final_df, paste0(CONFIG$simulated_dir, filename, ".csv"))
  cat(filename," Saved!!!\n")
  
  # Return settings info from first simulation
  list(
    stage1_probs = all_results[[1]]$stage1_probs,
    stage2_probs = all_results[[1]]$stage2_probs,
    MSE = mse_list
  )
}

# Main simulation loop
run_full_simulation <- function() {
  study_ids <- unique(CONFIG$data$Study)
  combinations <- get_valid_combinations(CONFIG$knots_list)
  settings_list <- list()
  
  for (stage2_knots in CONFIG$stage2_knots_list) {
    for (combo_idx in seq_len(nrow(combinations))) {
      stage1_knots_vec <- as.numeric(combinations[combo_idx, ])
      
      # filename
      filename <- paste0(
        "study", paste0(stage1_knots_vec, collapse = ""),
        "_second", stage2_knots
      )
      
      # Run all simulations
      all_results <- lapply(
        seq_len(CONFIG$n_sim),
        function(sim_i) run_single_simulation(
          sim_i, study_ids, stage1_knots_vec, stage2_knots
        )
      )
      
      # Save results
      results_info <- save_simulation_results(all_results, filename)
      
      # Collect settings
      prob_df <- data.frame(Name = filename)
      
      # Stage1 probabilities
      for (s in study_ids) {
        prob_df[paste0("Study", s)] <- paste(
          results_info$stage1_probs[[s]],
          collapse = ","
        )
      }
      
      prob_df$Stage2 <- paste(results_info$stage2_probs, collapse = ",")
      
      # Add MSE values
      for (col in names(results_info$MSE)) {
        prob_df[[col]] <- results_info$MSE[[col]]
      }
      
      settings_list[[filename]] <- unlist(prob_df[1, ])
    }
  }
  
  cat("Simulation completed!!!\n")
  
  # After loop — build settings dataframe
  if (length(settings_list) > 0) {
    settings_df <- do.call(rbind, settings_list)
    settings_df <- as.data.frame(settings_df, stringsAsFactors = FALSE)
    
    # Convert MSE columns to numeric
    mse_cols <- grep("^MSE_", names(settings_df), value = TRUE)
    settings_df[mse_cols] <- lapply(settings_df[mse_cols], as.numeric)
    
    # Add ranking columns (smaller MSE = better)
    for (col in mse_cols) {
      rank_col <- sub("MSE_", "Rank_", col)
      settings_df[[rank_col]] <- rank(settings_df[[col]], ties.method = "min")
    }
    
    # Save to CSV inside the correct simulated_dir
    write_csv(settings_df, paste0(CONFIG$simulated_dir, "simul_settings.csv"))
    cat("Settings saved!\n")
  }
  
  return(settings_df)
}

# Create dir if don't exist
dir.create(CONFIG$setting_dir, recursive = TRUE, showWarnings = FALSE)

for (cond in CONFIG$overlap_conditions) {
  
  CONFIG$simulated_dir <- paste0(CONFIG$setting_dir, cond, "/")
  dir.create(CONFIG$simulated_dir, recursive = TRUE, showWarnings = FALSE)
  
  input <- OVERLAP_INPUTS[[cond]]
  
  set.seed(2025)
  ipd_data <- generate_from_table(input) %>%
    rename(Age = X1) %>%
    mutate(Y = abs(Y), Study = as.factor(Study))
  write_csv(ipd_data, paste0(CONFIG$simulated_dir, "generated_ipd.csv"))
  
  CONFIG$data <- ipd_data
  
  dd <- datadist(CONFIG$data)
  options(datadist = "dd")
  
  cat("Running:", cond, "\n")
  
  # Run your simulation for this dataset
  settings_df <- run_full_simulation()
}
