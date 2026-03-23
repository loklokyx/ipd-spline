library(ggplot2)

# Define a custom theme
my_theme <- theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 14),           # X-axis ticks
    axis.text.y  = element_text(size = 14),           # Y-axis ticks
    axis.title.x = element_text(size = 16, face = "bold"), # X-axis label
    axis.title.y = element_text(size = 16, face = "bold"), # Y-axis label
    plot.title   = element_text(size = 18, face = "bold"), # Main title
    plot.subtitle = element_text(size = 16),                # Subtitle
    plot.caption  = element_text(size = 12),                # Caption
    legend.title  = element_text(size = 14, face = "bold"), # Legend title
    legend.text   = element_text(size = 12),                # Legend labels
    strip.text    = element_text(size = 14)
  )

# Apply it globally
theme_set(my_theme)

make_quantiles <- function(x, knot, dp = 4) {
  quantile_map <- list(
    '3' = c(0.10, 0.50, 0.90),
    '4' = c(0.05, 0.35, 0.65, 0.95),
    '5' = c(0.05, 0.275, 0.50, 0.725, 0.95),
    '6' = c(0.05, 0.23, 0.41, 0.59, 0.77, 0.95),
    '7' = c(0.025, 0.1833, 0.3417, 0.50, 0.6583, 0.8167, 0.975)
  )
  
  probs <- quantile_map[[as.character(knot)]]
  if (is.null(probs)) return(NULL)
  
  round(quantile(x, probs = probs, na.rm = TRUE), dp)
}

make_rcs_formula <- function(var = "Age", quantiles = NULL) {
  if (is.null(quantiles)) return(sprintf("rcs(%s)", var))
  sprintf("rcs(%s, c(%s))", var, paste(quantiles, collapse = ", "))
}

get_max_knots <- function(n) {
  if (n < 50) return(3)
  if (n < 200) return(4)
  return(5)
}

is_valid_knot <- function(n, knot) {
  knot <= get_max_knots(n)
}

knots_to_numeric <- function(knots_str) {
  # Handle NA or empty string
  if (is.na(knots_str) || knots_str == "") {
    return(c(0))
  }
  # Split by comma, trim whitespace, and convert to numeric
  as.numeric(trimws(unlist(strsplit(knots_str, ","))))
}

calculate_metrics <- function(df, level = 0.95) {
  # Identify simulation columns
  sim_cols <- grep("^sim_", names(df), value = TRUE)
  Y_true <- df$Y
  
  alpha <- 1 - level
  z_val <- qnorm(1 - alpha / 2)
  
  # Extract simulation matrix for easier computation
  sim_matrix <- as.matrix(df[, sim_cols])
  row_means <- rowMeans(sim_matrix)
  row_vars <- apply(sim_matrix, 1, var)
  row_mses <- rowMeans((sim_matrix - Y_true)^2)
  
  tibble(
    Age = df$Age,
    Y_true = Y_true,
    bias = row_means - Y_true,
    variance = row_vars,
    mse = row_mses,
    mean_pred = row_means,
    residual = Y_true - row_means,
    lower_ci = row_means - z_val * sqrt(row_vars),
    upper_ci = row_means + z_val * sqrt(row_vars),
    stage2_quantiles_str = df$stage2_quantiles_str[1],
    stage2_knots = df$stage2_knots,
    study_knots = df$study_knots,
    rank_method = df$rank_method,
    rank_id = df$rank_id
  )
}

# Function to check if Stage2 == all Study columns
is_same_knots <- function(row) {
  all(sapply(study_cols, function(c){
    length(knots_to_numeric(row[[c]])) == length(knots_to_numeric(row[["Stage2"]]))
  }))
}
is_same_knots <- function(row) {
  study_lengths <- sapply(study_cols, function(c) {
    length(knots_to_numeric(row[[c]]))
  })
  length(unique(study_lengths)) == 1   # TRUE if all equal, FALSE if not
}

library(stringr)
arrange_knots <- function(df, base_knot = NULL) {
  # Extract first part and remove "study"
  codes <- str_remove(str_split_fixed(df$Name, "_", 2)[,1], "study")
  # Split each code into digits
  code_digits <- str_split(codes, "", simplify = TRUE)
  # Determine unique knots used
  knots_settings <- sort(unique(as.numeric(as.vector(code_digits))))
  # Precompute number of occurrences of each knot in each code
  knot_counts <- sapply(knots_settings, function(k) rowSums(code_digits == k))
  colnames(knot_counts) <- knots_settings
  
  final_order <- c()
  # Loop over knots in ascending order
  for(k in knots_settings) {
    # Codes that contain this knot
    sub_codes <- df$Name[knot_counts[, as.character(k)] > 0]
    # Remove already added codes
    sub_codes <- setdiff(sub_codes, final_order)
    max_n <- ncol(code_digits)
    # Loop from max_n down to 1
    for(n in max_n:1) {
      codes_n <- sub_codes[knot_counts[match(sub_codes, df$Name), as.character(k)] == n]
      final_order <- c(final_order, codes_n)
      sub_codes <- setdiff(sub_codes, codes_n)
    }
  }
  
  df_ordered <- df %>% slice(match(final_order, Name))
  # optional grouping
  if(!is.null(base_knot)) {
    counts <- knot_counts[
      match(df_ordered$Name, df$Name),
      as.character(base_knot)
    ]
    df_ordered$group_base <- paste0(counts, " of k = ", base_knot) 
  }
  return(df_ordered)
}

plot_CI_coverage <- function(df, arrange_by = "name", base_knot = NULL,
                             stat_pos = NULL){
  library(ggplot2)
  library(dplyr)
  
  mean_pooled <- mean(df$Coverage_pooled, na.rm = TRUE)
  range_pooled <- range(df$Coverage_pooled, na.rm = TRUE)
  mean_meta <- mean(df$Coverage_meta, na.rm = TRUE)
  range_meta <- range(df$Coverage_meta, na.rm = TRUE)
  
  arrange_by <- tolower(arrange_by)
  if(arrange_by == "name") {
    df <- df %>% arrange(Name)
  } else if(arrange_by == "knots") {
    df <- arrange_knots(df, base_knot)
  }
  
  df$Name <- factor(df$Name, levels = df$Name)
  n <- length(df$Name)
  my_form <- function(m, r) sprintf("mean: %.3f [%.3f, %.3f]", m, r[1], r[2])
  
  p <- ggplot(df) + 
    geom_line(aes(x = Name, y = Coverage_pooled, group = 1, color = "pooled")) + 
    geom_line(aes(x = Name, y = Coverage_meta, group = 1, color = "meta")) +
    scale_color_manual(values = c("pooled" = "red", "meta" = "blue")) +
    theme(axis.text.x = element_text(angle = 90)) +
    ylab("95% CI Coverage")
  
  # add grouping only if requested
  if(arrange_by == "knots" && !is.null(base_knot)) {
    p <- p + 
      geom_point(aes(x = Name, y = Coverage_pooled, fill = factor(group_base)),
                 shape = 21, size = 2.5, color = "black") +
      geom_point(aes(x = Name, y = Coverage_meta, fill = factor(group_base)),
                 shape = 21, size = 2.5, color = "black") +
      labs(fill = "Stage 1 Knot Counts")
  }
  
  if(!is.null(stat_pos)) {
    stat_pos <- tolower(stat_pos)
    pos_index <- switch(stat_pos,
                        "start"  = 1 + 8, 
                        "middle" = ceiling(n / 2), 
                        "end"    = n - 8, 0)
    p <- p +
      annotate("text", x = pos_index, y = mean_pooled, color = "darkgreen",
               label = my_form(mean_pooled, range_pooled)) + 
      annotate("text", x = pos_index, y = mean_meta, color = "darkgreen",
               label = my_form(mean_meta, range_meta))
  }
  p
}

# Base plot template
base_plot <- function(df, metric_name, ncol, fix_y = TRUE, 
                      legend = TRUE, legend_ncol=1, legend_size=14) {
  quan <- make_quantiles(df$Age, df$stage2_knots[1])
  df_sub <- df %>% filter(metric_name == !!metric_name)
  
  p <- ggplot(df_sub , aes(x = Age, y = metric_value, color = line_id)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = quan, linetype = "dashed", color = "grey40") +
    facet_wrap(~rank_id, ncol = ncol, scales = if (fix_y) "fixed" else "free_y") +
    labs(
      x = "Age",
      y = metric_name,
      color = "Method"
    ) +
    theme(strip.background = element_rect(fill = "white"))
  
  if (legend) {
    p <- p + 
      guides(color = guide_legend(ncol = legend_ncol)) +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = legend_size)
      )
  } else {
    p <- p + theme(legend.position = "none")
  }
}

base_plot_CI <- function(df, ncol, fix_y = TRUE, title = "Predicted Values CI",
                         legend = TRUE, legend_ncol=1, legend_size=14,
                         jitter_width = 0, jitter_height = 0, show_jitter = FALSE) {
  quan <- make_quantiles(df$Age, df$stage2_knots[1])

  p <- ggplot(df, aes(x = Age, y = mean_pred, color = line_id, fill = line_id, group = line_id))
  
  if (show_jitter) {
    p <- p + geom_jitter(aes(y = Y_true), color = "black", alpha = 0.1, shape = 1,
                         width = jitter_width, height = jitter_height)
  }
  
  p <- p + geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = quan, linetype = "dashed", color = "grey40") +
    facet_wrap(~rank_id, ncol = ncol, scales = if (fix_y) "fixed" else "free_y") +
    labs(
      x = "Age",
      y = "Predicted Y",
      color = "Method",
      fill = "Method",
      title = title
    ) +
    guides(color = guide_legend(ncol = legend_ncol), fill = guide_legend(ncol = legend_ncol)) +
    theme(
      strip.background = element_rect(fill = "white"),
      panel.spacing = unit(1.5, "lines")
    )
  
  if (legend) {
    p <- p + 
      guides(color = guide_legend(ncol = legend_ncol)) +
      theme(
        legend.position = "bottom",
        legend.text   = element_text(size = legend_size)
      )
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  if (!show_jitter) {
    p <- p + geom_line(aes(y = Y_true), color = "black", linetype = "longdash")
  }
  
  return(p)
}

set_summary <- function(df_list, 
                        metrics = c("bias", "variance", "mse"), 
                        dp = 6,
                        prefix = NULL){
  full_df <- data.frame()
  for (name in names(df_list)) {
    df <- df_list[[name]]
    
    group_vars <- c("method", 
                    if ("knots" %in% names(df)) "knots" else NULL,
                    if ("rank_id" %in% names(df)) "rank_id" else NULL)
    
    df <- df %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(
        across(all_of(metrics),
               list(mean = ~round(mean(.x, na.rm = TRUE), dp),
                    sd   = ~round(sd(.x, na.rm = TRUE), dp),
                    min  = ~round(min(.x, na.rm = TRUE), dp),
                    max  = ~round(max(.x), dp)),
               .names = "{.fn}_{.col}"),
        .groups = "drop"
      )
    if (!is.null(prefix)) {
      df <- df %>% mutate(method = paste0(method, "_", prefix))
    }
    full_df <- rbind(full_df, df)
  }
  return(full_df)
}

make_rank_all <- function(simul_settings, top_n = 1, worst_n = 1){
  # Identify ranking columns
  rank_cols <- grep("^Rank_", names(simul_settings), value = TRUE)
  
  rank_all <- lapply(rank_cols, function(col) {
    
    # ----- TOP -----
    if (top_n > 0) {
      best <- simul_settings %>%
        arrange(.data[[col]]) %>%
        slice_head(n = top_n) %>%
        mutate(
          rank_method = sub("^Rank_", "", col),
          rank_id = paste0("best", seq_len(top_n))
        )
    } else {
      best <- NULL
    }
    
    # ----- WORST -----
    if (worst_n > 0) {
      worst <- simul_settings %>%
        arrange(.data[[col]]) %>%
        slice_tail(n = worst_n) %>%
        mutate(
          rank_method = sub("^Rank_", "", col),
          rank_id = paste0("worst", worst_n:1)
        )
    } else {
      worst <- NULL
    }
    
    bind_rows(best, worst)
  }) %>%
    bind_rows()
  
  # ----- Factor levels -----
  lvls <- c(
    if (top_n > 0) paste0("best", seq_len(top_n)),
    if (worst_n > 0) paste0("worst", worst_n:1)
  )
  
  if (length(lvls) > 0 && nrow(rank_all) > 0) {
    rank_all <- rank_all %>%
      mutate(rank_id = factor(rank_id, levels = lvls))
  }
  
  return(rank_all)
}

myfunc <- "10 + 5 * (1 + exp((.15 * (X1 - 60))))^(-1)"

OVERLAP_INPUTS <- list(
  
  simulated_some = tribble(
    ~study, ~n, ~noise, ~Y, ~X1,
    "1", 300, 0, myfunc, "rnorm(n, 35, 1)",
    "2", 300, 0, myfunc, "rnorm(n, 40, 2)",
    "3", 300, 0, myfunc, "rnorm(n, 50, 3)",
    "4", 300, 0, myfunc, "rnorm(n, 60, 4)",
    "5", 300, 0, myfunc, "rnorm(n, 70, 5)",
    "6", 300, 0, myfunc, "rnorm(n, 80, 5)"
  ),
  
  simulated_no = tribble(
    ~study, ~n, ~noise, ~Y, ~X1,
    "1", 300, 0, myfunc, "rnorm(n, 30, 1)",
    "2", 300, 0, myfunc, "rnorm(n, 30, 2)",
    "3", 300, 0, myfunc, "rnorm(n, 50, 3)",
    "4", 300, 0, myfunc, "rnorm(n, 50, 4)",
    "5", 300, 0, myfunc, "rnorm(n, 80, 5)",
    "6", 300, 0, myfunc, "rnorm(n, 80, 5)"
  ),
  
  simulated_many = tribble(
    ~study, ~n, ~noise, ~Y, ~X1,
    "1", 100, 0, myfunc, "rnorm(n, 40, 2)",
    "1", 100, 0, myfunc, "rnorm(n, 50, 2)",
    "1", 100, 0, myfunc, "rnorm(n, 60, 2)",
    "2", 100, 0, myfunc, "rnorm(n, 35, 4)",
    "2", 100, 0, myfunc, "rnorm(n, 50, 4)",
    "2", 100, 0, myfunc, "rnorm(n, 75, 4)",
    "3", 100, 0, myfunc, "rnorm(n, 35, 3)",
    "3", 100, 0, myfunc, "rnorm(n, 45, 3)",
    "3", 100, 0, myfunc, "rnorm(n, 55, 3)",
    "4", 100, 0, myfunc, "rnorm(n, 55, 4)",
    "4", 100, 0, myfunc, "rnorm(n, 65, 4)",
    "4", 100, 0, myfunc, "rnorm(n, 80, 3)",
    "5", 100, 0, myfunc, "rnorm(n, 65, 5)",
    "5", 100, 0, myfunc, "rnorm(n, 70, 5)",
    "5", 100, 0, myfunc, "rnorm(n, 75, 4)",
    "6", 100, 0, myfunc, "rnorm(n, 55, 5)",
    "6", 100, 0, myfunc, "rnorm(n, 70, 5)",
    "6", 100, 0, myfunc, "rnorm(n, 85, 3)"
  )
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

# Stage 2: Pooled method
run_pooled <- function(pred_data, knots = 5, weight_fun = NULL,
                       interval = "confidence", combine = FALSE) {
  pred_str <- make_rcs_formula("Age", make_quantiles(pred_data$Age, knots))
  pooled_formula <- as.formula(paste0("Y ~ ", pred_str))
  
  if (is.null(weight_fun)) {
    d <- pred_data
    pooled_model <- ols(pooled_formula, data = d)
  } else {
    d <- pred_data %>% mutate(weight = weight_fun(SE))
    pooled_model <- ols(pooled_formula, weights = weight, data = d)
  }
  
  X <- model.matrix(delete.response(terms(pooled_model)), data = d)
  beta <- coef(pooled_model)
  V <- as.matrix(vcov(pooled_model))
  
  fit <- as.numeric(X %*% beta)
  # se_fit <- sqrt(diag(X %*% V %*% t(X)))
  se_fit <- sqrt(rowSums((X %*% V) * X))
  
  if (interval == "response") {
    resid_sd <- sqrt(sum(residuals(pooled_model)^2) / pooled_model$df.residual)
    se_total <- sqrt(se_fit^2 + resid_sd^2)
  } else if(interval == "response_s") {
    resid_sd <- pooled_model$stats["Sigma"]
    se_total <- sqrt(se_fit^2 + resid_sd^2)
  } else {
    se_total <- se_fit
    # stage1_var <- tapply(pred_data$SE^2, pred_data$Age, mean)
    # se_total <- sqrt(se_fit^2 + stage1_var[match(pred_data$Age, names(stage1_var))])
  }
  
  d$Y  <- fit
  d$lower <- fit - 1.96 * se_total
  d$upper <- fit + 1.96 * se_total
  d$SE <- se_total
  if (combine) d <- combine_pointwise(d, weight_fun = weight_fun)
  d %>% slice(match(unique(d$Age), Age))
}

# Stage 2: Meta-analysis fixed effects with weights
run_meta_fixed <- function(pred_data, knots = 5, weight_fun = NULL,
                           interval = "confidence", combine = FALSE) {
  pred_str <- make_rcs_formula("Age", make_quantiles(pred_data$Age, knots))
  meta_formula <- as.formula(paste0("Y ~ ", pred_str, " + (1 | study)"))

  if (is.null(weight_fun)) {
    d <- pred_data
    meta_model <- lmer(meta_formula, data = d)
  } else {
    d <- pred_data %>% mutate(weight = weight_fun(SE))
    meta_model <- lmer(meta_formula, weights = weight, data = d)
  }
  
  X <- model.matrix(delete.response(terms(meta_model)), d)
  beta <- fixef(meta_model)
  V <- vcov(meta_model)
  
  fit <- as.numeric(X %*% beta)
  se_fixed <- sqrt(rowSums((X %*% V) * X))
  # se_fixed <- sqrt(diag(X %*% V %*% t(X)))
  
  if (interval == "response"){
    se_total <- sqrt(se_fixed^2 + sigma(meta_model)^2)
  } else if (interval == "random"){
    var_random <- as.numeric(VarCorr(meta_model)$study)
    se_total <- sqrt(se_fixed^2 + var_random)
  } else{
    se_total <- se_fixed
  }
  length(se_total)
  d$Y <- fit
  d$lower <- fit - 1.96 * se_total
  d$upper <- fit + 1.96 * se_total
  d$SE <- se_total
  if (combine) d <- combine_pointwise(d, weight_fun = weight_fun)
  d %>% slice(match(unique(d$Age), Age))
}

combine_pointwise <- function(pred_data, weight_fun = NULL) {
  d <- pred_data
  
  if (!is.null(weight_fun)) {
    d <- d %>% mutate(w = weight_fun(SE))
    
    d_summary <- d %>%
      group_by(Age) %>%
      summarise(
        Y_pred = sum(Y * w) / sum(w),
        var_weighted = 1 / sum(w),
        # var_weighted = sum(w^2 * SE^2)/sum(w)^2 + 1/sum(w),
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
    arrange(Age) %>% 
    slice(match(unique(d$Age), Age)) # reorder by original Age
}
