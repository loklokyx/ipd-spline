library(tidyverse)
library(dplyr)
source("data/random_data.r")

# Define study parameters
# input <- tribble(
#   ~study, ~n, ~noise, ~Y,                    ~X1,     ~X2,
#   "1",    8,  0, "rnorm(n, 400, 50)",   "25:35", "runif(n, 18, 25)",
#   "2",    6,  0, "rnorm(n, 500, 70)",   "30:45", "runif(n, 20, 30)",
#   "3",    10, 0, "rnorm(n, 450, 60)",   "20:40", "runif(n, 22, 28)"
# )

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
  rename(Age = X1) %>%
  mutate(Y = abs(Y))

write_csv(ipd_data, "data/generated_ipd.csv")
