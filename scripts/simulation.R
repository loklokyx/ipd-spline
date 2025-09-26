library(gtools)
library(tidyverse)
library(rms)
library(lme4)
source("data/random_data.r")
source("scripts/helpers.R")

# Data generation
myfunc <- "10 + 5 * (1 + exp((.15 * (X1 - 60))))^(-1)"
# # some overlap
# input <- tribble(
#   ~study, ~n, ~noise, ~Y, ~X1,
#   "1", 300, 0, myfunc, "rnorm(n, 35, 1)",
#   "2", 300, 0, myfunc, "rnorm(n, 40, 2)",
#   "3", 300, 0, myfunc, "rnorm(n, 50, 3)",
#   "4", 300, 0, myfunc, "rnorm(n, 60, 4)",
#   "5", 300, 0, myfunc, "rnorm(n, 70, 5)",
#   "6", 300, 0, myfunc, "rnorm(n, 80, 5)"
# )

# no overlap 40 and 60
input <- tribble(
  ~study, ~n, ~noise, ~Y, ~X1,
  "1", 300, 0, myfunc, "rnorm(n, 30, 1)",
  "2", 300, 0, myfunc, "rnorm(n, 30, 2)",
  "3", 300, 0, myfunc, "rnorm(n, 50, 3)",
  "4", 300, 0, myfunc, "rnorm(n, 50, 4)",
  "5", 300, 0, myfunc, "rnorm(n, 80, 5)",
  "6", 300, 0, myfunc, "rnorm(n, 80, 5)"
)

# # many overlap
# input <- tribble(
#   ~study, ~n, ~noise, ~Y, ~X1,
#   "1", 100, 0, myfunc, "rnorm(n, 40, 2)",
#   "1", 100, 0, myfunc, "rnorm(n, 50, 2)",
#   "1", 100, 0, myfunc, "rnorm(n, 60, 2)",
#   "2", 100, 0, myfunc, "rnorm(n, 35, 4)",
#   "2", 100, 0, myfunc, "rnorm(n, 50, 4)",
#   "2", 100, 0, myfunc, "rnorm(n, 75, 4)",
#   "3", 100, 0, myfunc, "rnorm(n, 35, 3)",
#   "3", 100, 0, myfunc, "rnorm(n, 45, 3)",
#   "3", 100, 0, myfunc, "rnorm(n, 55, 3)",
#   "4", 100, 0, myfunc, "rnorm(n, 55, 4)",
#   "4", 100, 0, myfunc, "rnorm(n, 65, 4)",
#   "4", 100, 0, myfunc, "rnorm(n, 80, 3)",
#   "5", 100, 0, myfunc, "rnorm(n, 65, 5)",
#   "5", 100, 0, myfunc, "rnorm(n, 70, 5)",
#   "5", 100, 0, myfunc, "rnorm(n, 75, 4)",
#   "6", 100, 0, myfunc, "rnorm(n, 55, 5)",
#   "6", 100, 0, myfunc, "rnorm(n, 70, 5)",
#   "6", 100, 0, myfunc, "rnorm(n, 85, 3)"
# )

set.seed(2025)
ipd_data <- generate_from_table(input) %>%
  rename(Age = X1) %>%
  mutate(Y = abs(Y), Study = as.factor(Study))
# plot(Y~Age, col=Study,data=ipd_data)

write_csv(ipd_data, "data/generated_ipd.csv")

# Configuration - keep all settings together
CONFIG <- list(
  data = ipd_data,
  n_sim = 1000,
  knots_list = c(5, 4),
  stage2_knots_list = c(5),
  noise_sd = 2,
  include_meta_re = TRUE,
  include_meta_re = F,
  setting_dir = "./data/simulated/",
  simulated_dir = "./data/simulated/"
)

dd <- datadist(CONFIG$data)
options(datadist = "dd")

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
    model <- ols(formula_obj, data = study_data)
    
    temp <- data.frame(
      Study = s,
      Age = study_data$Age,
      Y = predict(model, newdata = study_data)
    )
    pred_stage1 <- rbind(pred_stage1, temp)
  }
  
  list(data = pred_stage1, probs = stage1_probs)
}

# Stage 2: Meta-analysis methods
run_stage2 <- function(pred_stage1, stage2_knots) {
  probs2 <- make_quantiles(pred_stage1$Age, stage2_knots)
  stage2_pred_str <- make_rcs_formula("Age", probs2)
  
  results <- list()
  
  # Method 1: Pooled
  pooled_data <- pred_stage1
  pooled_data$Study <- 1
  pooled_formula <- as.formula(paste0("Y ~ ", stage2_pred_str))
  pooled_model <- ols(pooled_formula, data = pooled_data)
  results$pooled <- predict(pooled_model, newdata = CONFIG$data)
  
  # Method 2 & 3: Meta-analysis
  meta_formula <- as.formula(paste0("Y ~ ", stage2_pred_str, " + (1 | Study)"))
  meta_model <- lmer(meta_formula, data = pred_stage1)
  results$meta <- predict(meta_model, newdata = CONFIG$data, re.form = NA)
  
  if (CONFIG$include_meta_re) {
    results$meta_re <- predict(meta_model, newdata = CONFIG$data, re.form = ~(1 | Study))
  }
  
  list(results = results, stage2_probs = probs2)
}

# Single simulation run
run_single_simulation <- function(sim_i, study_ids, stage1_knots_vec, stage2_knots) {
  set.seed(sim_i)
  sim_data <- CONFIG$data %>% mutate(Y = Y + rnorm(n(), 0, CONFIG$noise_sd))
  
  stage1_results <- run_stage1(sim_data, study_ids, stage1_knots_vec)
  stage2_results <- run_stage2(stage1_results$data, stage2_knots)
  
  list(
    predictions = stage2_results$results,
    stage1_probs = if (sim_i == 1) stage1_results$probs else NULL,
    stage2_probs = if (sim_i == 1) stage2_results$stage2_probs else NULL
  )
}

# Format and save results
save_simulation_results <- function(all_results, filename) {
  # Initialize method results storage
  method_names <- names(all_results[[1]]$predictions)
  mse_list <- list()   # <-- initialize mse_list!
  
  # Compute weighted MSE for each method
  Y_true <- CONFIG$data$Y
  # Save simulated predictions
  row_based_df <- data.frame()
  for (method in method_names) {
    sim_matrix <- do.call(cbind, lapply(all_results, function(res) res$predictions[[method]]))
    
    temp_df <- data.frame(method = method, sim_matrix)
    row_based_df <- rbind(row_based_df, temp_df)

    row_mses <- rowMeans((sim_matrix - Y_true)^2)
    overall_mse <- mean(row_mses)  # average MSE for this setting

    # # Optional: use equal weight or weighted by study size
    # weights <- table(CONFIG$data$Study)[as.character(CONFIG$data$Study)]
    # weights <- as.numeric(weights)
    # # Weighted overall MSE
    # overall_mse <- sum(weights * row_mses) / sum(weights)
    mse_list[[paste0("MSE_", method)]] <- overall_mse
  }
  
  # Rename sim columns
  colnames(row_based_df) <- c("method", paste0("sim_", seq_len(CONFIG$n_sim)))
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
study_ids <- unique(CONFIG$data$Study)
combinations <- get_valid_combinations(CONFIG$knots_list)
settings_list <- list()

for (stage2_knots in CONFIG$stage2_knots_list) {
  for (combo_idx in seq_len(nrow(combinations))) {
    stage1_knots_vec <- as.numeric(combinations[combo_idx, ])

    # # for same knots validation before run all combinations
    # if (!all(stage1_knots_vec==stage2_knots)) next

    filename <- paste0("study", paste0(stage1_knots_vec, collapse = ""), 
                       "_second", stage2_knots)
    
    # Run all simulations for this combination
    all_results <- list()
    for (sim_i in seq_len(CONFIG$n_sim)) {
      all_results[[sim_i]] <- run_single_simulation(
        sim_i, study_ids, stage1_knots_vec, stage2_knots
      )
    }

    # Save results and collect settings
    results_info <- save_simulation_results(all_results, filename)

    # Store settings with MSEs
    prob_df <- data.frame(Name = filename)
    for (s in study_ids) {
      prob_df[paste0("Study", s)] <- paste(results_info$stage1_probs[[s]], collapse = ",")
    }
    prob_df$Stage2 <- paste(results_info$stage2_probs, collapse = ",")
    
    # Add all MSEs
    for (col in names(results_info$MSE)) {
      prob_df[[col]] <- results_info$MSE[[col]]
    }
    settings_list[[filename]] <- unlist(prob_df[1, ])
  }
}
cat("Simulation completed!!!")

# Save settings with ranks
if (length(settings_list) > 0) {
  settings_df <- do.call(rbind, settings_list)
  settings_df <- as.data.frame(settings_df, stringsAsFactors = FALSE)
  
  # Convert MSE columns to numeric
  mse_cols <- grep("^MSE_", names(settings_df), value = TRUE)
  settings_df[mse_cols] <- lapply(settings_df[mse_cols], as.numeric)
  
  # Add rank for each method based on MSE (smaller is better)
  for (col in mse_cols) {
    rank_col <- sub("MSE_", "Rank_", col)
    settings_df[[rank_col]] <- rank(settings_df[[col]], ties.method = "min")
  }
  
  # Save to CSV
  write_csv(settings_df, paste0(CONFIG$setting_dir, "simul_settings.csv"))
  cat("Setting saved!!!")
}