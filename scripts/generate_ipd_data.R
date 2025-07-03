library(tidyverse)
source("data/random_data.r")

# Define study parameters
input <- tribble(
  ~study, ~n, ~Y,                    ~X1,     ~X2,
  "1",    8,  "rnorm(n, 400, 50)",   "25:35", "runif(n, 18, 25)",
  "2",    6,  "rnorm(n, 500, 70)",   "30:45", "runif(n, 20, 30)",
  "3",    10, "rnorm(n, 450, 60)",   "20:40", "runif(n, 22, 28)"
)

# input <- tribble(
#   ~study, ~n, ~Y,                    ~X1,     
#   "1",    8,  "rnorm(n, 400, 50)",   "25:35", 
#   "2",    6,  "rnorm(n, 500, 70)",   "30:45", 
#   "3",    10, "rnorm(n, 450, 60)",   "10:30", 
#   "4",    7,  "rnorm(n, 450, 60)",   "30:55", 
#   "5",    9,  "rnorm(n, 450, 60)",   "29:59"
# )


set.seed(2025)
ipd_data <- generate_from_table(input)

ipd_data <- ipd_data %>%
  rename(Age = X1, BMI = X2)

# ipd_data <- ipd_data %>%
#   mutate(
#     Y = case_when(
#       Study == "1" ~ 2.0 + 5*Age^2 + rnorm(n(), 0, 1.0),    
#       Study == "2" ~ 3.0 + 2*Age^2 + rnorm(n(), 0, 2.0),  
#       Study == "3" ~ 1.5 + 3*Age^2 + rnorm(n(), 0, 1.5),  
#       Study == "4" ~ 0.5 + 6*Age^2 + rnorm(n(), 0, 2.5),  
#       Study == "5" ~ 2.5 + 2*Age^2 + rnorm(n(), 0, 1.5),
#       TRUE ~ NA_real_
#     ),
#     Y = round(Y, 4)
#   )

write_csv(ipd_data, "data/generated_ipd.csv")
