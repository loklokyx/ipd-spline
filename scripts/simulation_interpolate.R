library(gtools)
library(tidyverse)
library(rms)
library(lme4)
source("data/random_data.r")
source("scripts/helpers.R")
library(future.apply)
plan(multisession, workers = parallel::detectCores() - 1)

set.seed(2025)

# Data generation
myfunc <- "10 + 5 * (1 + exp((.15 * (X1 - 60))))^(-1)"

# Configuration - keep all settings together
CONFIG <- list(
  data = NA,
  n_sim = 20,
  knots_list = c(4, 5),
  overlap_conditions = c("simulated_no", "simulated_some", "simulated_many"),
  # n_sim = 10,
  # knots_list = c(4),
  # overlap_conditions = c("simulated_many"),
  stage2_knots_list = c(5),
  noise_sd = 2,
  weight_fun = function(se) 1/se^2,  # inverse variance weighting
  setting_dir = "./data/interpolate_coverage1/"
)

# Stage 1: Individual study predictions
run_stage1 <- function(sim_data, study_ids, stage1_knots_vec) {
  stage1_probs <- list()

  pred_stage1_list <- vector("list", length(study_ids))
  
  for (s_idx in seq_along(study_ids)) {
    s <- study_ids[s_idx]
    study_data <- subset(sim_data, Study == s)
    knots_s <- stage1_knots_vec[s_idx]
    
    probs <- make_quantiles(study_data$Age, knots_s)
    stage1_probs[[as.character(s)]] <- probs
    
    pred_str <- make_rcs_formula("Age", probs)
    formula_obj <- as.formula(paste0("Y ~ ", pred_str))
    model <- ols(formula_obj, data = study_data)
    
    X <- model.matrix(delete.response(terms(model)), study_data)
    beta <- coef(model)
    V <- as.matrix(vcov(model))
    
    fit <- as.numeric(X %*% beta)
    se_fit <- sqrt(rowSums((X %*% V) * X))
    
    # confidence interval only
    # se_total <- se_fit
    se_total <- sqrt(se_fit^2 + model$stats["Sigma"]^2)
    
    pred_stage1_list[[s_idx]] <- data.frame(
      Study = s,
      study = s,
      Age   = study_data$Age,
      Y     = fit,
      SE    = se_total
    )
  }
  
  pred_stage1 <- do.call(rbind, pred_stage1_list)
  list(data = pred_stage1, probs = stage1_probs)
}

# Single simulation run
run_single_simulation <- function(sim_i, study_ids, stage1_knots_vec, stage2_knots) {
  sim_data <- CONFIG$data %>% mutate(Y = Y + rnorm(n(), 0, CONFIG$noise_sd))
  
  stage1_results <- run_stage1(sim_data, study_ids, stage1_knots_vec)
  # Stage 2: Methods
  pooled_result <- run_pooled(stage1_results$data, CONFIG$stage2_knots, CONFIG$weight_fun, "response_s")
  meta_result <- run_meta_fixed(stage1_results$data, CONFIG$stage2_knots, CONFIG$weight_fun, "random")
  stage2_results <- list(pooled = pooled_result, meta = meta_result)

  stage2_probs <- make_quantiles(stage1_results$data$Age, CONFIG$stage2_knots)
  list(
    predictions = stage2_results,
    stage1_probs = if (sim_i == 1) stage1_results$probs else NULL,
    stage2_probs = if (sim_i == 1) stage2_probs else NULL
  )
}

# Format and save results
save_simulation_results <- function(all_results, filename, is_same) {
  # Initialize method results storage
  method_names <- names(all_results[[1]]$predictions)
  mse_list <- list()   # <-- initialize mse_list!
  
  # Compute weighted MSE for each method
  Y_true <- CONFIG$data$Y
  # Save simulated predictions
  row_based_df <- data.frame()
  coverage_df <- data.frame()
  for (method in method_names) {
    sim_matrix <- do.call(cbind, lapply(
      all_results, function(res) res$predictions[[method]][['Y']]))
    sim_lower <- do.call(cbind, lapply(
      all_results, function(res) res$predictions[[method]][['lower']]))
    sim_upper <- do.call(cbind, lapply(
      all_results, function(res) res$predictions[[method]][['upper']]))
    
    temp_df <- data.frame(method = method, sim_matrix)
    row_based_df <- rbind(row_based_df, temp_df)
    
    row_mses <- rowMeans((sim_matrix - Y_true)^2)
    overall_mse <- mean(row_mses)  # average MSE for this setting
    mse_list[[paste0("MSE_", method)]] <- overall_mse
    
    # Compute coverage using lower and upper bounds
    coverage_matrix <- (sim_lower <= Y_true) & (Y_true <= sim_upper)
    coverage <- mean(rowMeans(coverage_matrix))
    mse_list[[paste0("Coverage_", method)]] <- coverage
    
    if(is_same) {
      temp_coverage_df <- data.frame(
        method = method,
        Study = rep(CONFIG$data$Study, CONFIG$n_sim),
        Age = rep(CONFIG$data$Age, CONFIG$n_sim),
        Y = as.vector(sim_matrix),
        sim = rep(seq_len(CONFIG$n_sim), each = nrow(sim_lower)),
        lower = as.vector(sim_lower),
        upper = as.vector(sim_upper)
      )
      
      coverage_df <- rbind(coverage_df, temp_coverage_df)
    }
  }
  
  # Rename sim columns
  colnames(row_based_df) <- c("method", paste0("sim_", seq_len(CONFIG$n_sim)))
  final_df <- cbind(CONFIG$data, row_based_df)
  
  if(is_same) {
    write_csv(coverage_df, 
              paste0(CONFIG$simulated_dir, "coverage_", filename, ".csv"))
  }
  
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
      is_same <- length(unique(stage1_knots_vec)) == 1
      if(!is_same) next
      
      # filename
      filename <- paste0(
        "study", paste0(stage1_knots_vec, collapse = ""),
        "_second", stage2_knots
      )
      
      # Run all simulations
      all_results <- future_lapply(seq_len(CONFIG$n_sim), function(sim_i) {
        run_single_simulation(sim_i, study_ids, stage1_knots_vec, stage2_knots)
      }, future.seed = TRUE)
      
      # Save results
      results_info <- save_simulation_results(all_results, filename, is_same)
      
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
    cat("Settings saved in ", CONFIG$simulated_dir, " !\n")
  }
  
  return(settings_df)
}

# Create dir if don't exist
dir.create(CONFIG$setting_dir, recursive = TRUE, showWarnings = FALSE)

for (cond in CONFIG$overlap_conditions) {
  
  CONFIG$simulated_dir <- paste0(CONFIG$setting_dir, cond, "/")
  dir.create(CONFIG$simulated_dir, recursive = TRUE, showWarnings = FALSE)
  
  input <- OVERLAP_INPUTS[[cond]]
  
  CONFIG$data <- generate_from_table(input) %>%
    rename(Age = X1) %>%
    mutate(
      Y = abs(Y), 
      Study = as.factor(Study))
  
  write_csv(CONFIG$data, paste0(CONFIG$simulated_dir, "generated_ipd.csv"))
  
  cat("Running:", cond, "\n")
  dd <- datadist(CONFIG$data)
  options(datadist = "dd")
  # Run your simulation for this dataset
  settings_df <- run_full_simulation()
}
