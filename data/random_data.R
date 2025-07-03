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

generate_study_data <- function(study_id, n, Y_range, ...) {
  dp <- 5
  Y <- round(generate_values(Y_range, n), dp)
  # X variables (dynamic arguments)
  x_ranges <- list(...)
  x_data <- map(x_ranges, ~round(generate_values(.x, n), dp))
  
  # Create dataframe
  df <- data.frame(
    Study = study_id,
    Y = Y,
    x_data
  )
  
  return(df)
}

generate_from_table <- function(input_table) {
  all_data <- list()
  
  for (i in 1:nrow(input_table)) {
    row <- input_table[i, ]
    
    # Get X columns dynamically
    x_cols <- names(row)[grepl("^X", names(row))]
    x_ranges <- as.list(row[x_cols])
    names(x_ranges) <- x_cols
    
    study_data <- do.call(generate_study_data, c(
      list(study_id = row$study, n = row$n, Y_range = row$Y),
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
