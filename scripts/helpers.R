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
# Base plot template
base_plot <- function(df, metric_name, ncol, fix_y = TRUE, legend_ncol=1) {
  quan <- make_quantiles(df$Age, df$stage2_knots[1])
  df_sub <- df %>% filter(metric_name == !!metric_name)
  
  ggplot(df_sub, aes(x = Age, y = metric_value, color = line_id)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = quan, linetype = "dashed", color = "grey40") +
    facet_wrap(~rank_id, ncol = ncol, scales = if (fix_y) "fixed" else "free_y") +
    labs(
      x = "Age",
      y = metric_name,
      color = "Method"
    ) +
    guides(color = guide_legend(ncol = legend_ncol)) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    )
}

base_plot_CI <- function(df, ncol, fix_y = TRUE, legend_ncol=1, 
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
      title = "Predicted Values with Confidence Intervals by Rank"
    ) +
    guides(color = guide_legend(ncol = legend_ncol)) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white"),
      panel.spacing = unit(1.5, "lines")
    )
  
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
    best <- simul_settings %>%
      arrange(.data[[col]]) %>%
      slice_head(n = top_n) %>%
      mutate(rank_method = sub("^Rank_", "", col),
             rank_id = paste0("top", seq_len(top_n)))
    
    worst <- simul_settings %>%
      arrange(.data[[col]]) %>%
      slice_tail(n = worst_n) %>%
      mutate(rank_method = sub("^Rank_", "", col),
             rank_id = paste0("worst", worst_n:1))  # ensures worst1 = absolute worst
    
    bind_rows(best, worst)
  }) %>%
    bind_rows()
  
  # Define factor levels in correct order
  lvls <- c(paste0("top", 1:top_n), paste0("worst", worst_n:1))
  
  rank_all <- rank_all %>%
    mutate(rank_id = factor(rank_id, levels = lvls))
  
  return(rank_all)
}

make_unique_by_precision <- function(x, digits = 5) {
  # Round first to the desired digits
  x <- round(x, digits)
  
  # Small offset based on precision (e.g. digits=5 → 1e-6)
  offset <- 10^(-(digits + 1))
  
  # Order to preserve original sequence
  ord <- order(x)
  x_sorted <- x[ord]
  
  # Handle duplicates
  dup_index <- ave(x_sorted, x_sorted, FUN = seq_along)
  
  # Add incremental offset to duplicates
  x_sorted <- x_sorted + (dup_index - 1) * offset
  
  # Revert to original order
  x_final <- x_sorted[order(ord)]
  
  return(x_final)
}


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
