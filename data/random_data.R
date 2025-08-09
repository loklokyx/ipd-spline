library(tidyverse)

# handle simple strings or ranges
generate_values <- function(value_str, n) {
  # if is "min:max" range
  if (grepl("^[0-9.-]+:[0-9.-]+$", value_str)) {
    parts <- str_split(value_str, ":", simplify = TRUE)
    min <- as.numeric(parts[1])
    max <- as.numeric(parts[2])
    if (any(is.na(c(min, max)))) stop("Invalid numeric range")
    
    # if contain decimal
    has_decimal <- grepl("\\.", parts[1]) || grepl("\\.", parts[2])
    
    if (has_decimal) {
      # continuous uniform
      return(runif(n, min = min, max = max))
    } else {
      # as integers
      return(sample(seq(min, max), size = n, replace = TRUE))
    }
  } else {
    # replace 'n' to actual value
    func_call <- gsub("\\bn\\b", as.character(n), value_str)
    
    # evaluate the distribution code
    tryCatch({
      result <- eval(parse(text = func_call))
      if (length(result) != n) {
        stop(paste("Function call returned", length(result), "values, expected", n))
      }
      return(result)
    }, error = function(e) {
      stop(paste("Error evaluating function call '", func_call, "':", e$message))
    })
  }
}

generate_study_data <- function(study_id, n, Y_range, noise_sd, ...) {
  dp <- 5
  # X variables (dynamic arguments)
  x_ranges <- list(...)
  x_data <- map(x_ranges, ~round(generate_values(.x, n), dp))
  x_df <- as.data.frame(x_data)

  # detect if Y is an expression
  if (grepl("X", Y_range)) {
    env <- as.list(x_df)
    Y_expr <- parse(text = Y_range)
    Y <- eval(Y_expr, envir = env)
    
    # Add noise
    if (noise_sd > 0) {
      Y <- Y + rnorm(n, mean = 0, sd = noise_sd)
    }
  } else {
    Y <- generate_values(Y_range, n)
  }
  # Create dataframe
  df_final <- tibble(
    Study = study_id,
    Y = round(Y, dp)
  ) %>% bind_cols(x_df)
  
  return(df_final)
}

generate_from_table <- function(input_table) {
  all_data <- list()
  
  for (i in 1:nrow(input_table)) {
    row <- input_table[i, ]
    # row_list <- as.list(row)

    # Get X columns dynamically
    x_cols <- names(row)[grepl("^X", names(row))]
    x_ranges <- row[x_cols]
    names(x_ranges) <- x_cols
    
    study_data <- do.call(generate_study_data, c(
      list(
        study_id = row$study,
        n = row$n,
        Y_range = row$Y,
        noise_sd = row$noise
      ),
      x_ranges
    ))
    
    all_data[[i]] <- study_data
  }
  
  bind_rows(all_data)
}

# example
# input <- tribble(
#   ~study, ~n, ~Y,                    ~X1,              ~X2,
#   "1",    8,  "rbinom(n, 10, 0.5)",  "1:2", "runif(n, 3.5, 7.5)",
#   "2",    6,  "rexp(n, 0.1)",        "2:4",        "rpois(n, 7)",
#   "3",    10, "15:25",               "1:3", "rgamma(n, 2, 1)"
# )
# 
# set.seed(2025)
# generate_from_table(input)

# example (expression)
# input <- tribble(
#   ~study, ~n,       ~Y,   ~X1,
#   "1",    8,    "X1^2", "1.0:2.0",
#   "2",    6,  "X1^2+2", "2.0:4.0",
#   "3",    10, "X1^2+4", "1.0:3.0",
# )
# 
# set.seed(2025)
# generate_from_table(input)
