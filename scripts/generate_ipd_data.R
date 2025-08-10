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

myfunc  = "(1+exp(-.1*X1))^(-1)"
myfunc2 = "35*(1+exp(-X1+15))^(-1)"
myfunc3 = "30*(1+exp(.1*X1-7))^(-1)"

# input <- tribble(
#   ~study, ~n, ~noise,     ~Y,         ~X1,
#   "1", 20,    10, myfunc2, "10.0:20.0",
#   "1", 20,    10, myfunc2, "20.0:30.0",
#   
#   "2",  5,     2,  myfunc,  "0.0:10.0",
#   "2",  5,    10, myfunc2, "10.0:20.0",
#   "2",  5,    10, myfunc2, "20.0:30.0",
#   "2",  5,    10, myfunc3, "30.0:40.0",
#   "2",  5,    10, myfunc3, "40.0:50.0",
#   "2",  5,    10, myfunc3, "50.0:60.0",
#   "2",  5,    20, myfunc3, "60.0:70.0",
#   "2",  5,    20, myfunc3, "70.0:80.0",
#   
#   "3", 20,    10, myfunc2, "10.0:20.0",
#   "3", 20,    10, myfunc2, "20.0:30.0",
#   
#   "4", 10,     2,  myfunc,  "0.0:10.0",
#   "4", 10,    10, myfunc2, "10.0:20.0",
#   "4", 10,    10, myfunc2, "20.0:30.0",
#   "4", 10,    10, myfunc3, "30.0:40.0",
# 
#   "5",  5,     2,  myfunc,  "0.0:10.0",
#   "5",  5,    10, myfunc2, "10.0:20.0",
#   "5",  5,    10, myfunc2, "20.0:30.0",
#   "5",  5,    10, myfunc3, "30.0:40.0",
#   "5",  5,    10, myfunc3, "40.0:50.0",
#   "5",  5,    10, myfunc3, "50.0:60.0",
#   "5",  5,    20, myfunc3, "60.0:70.0",
#   "5",  5,    20, myfunc3, "70.0:80.0",
# )

# set.seed(2025)
# ipd_data <- generate_from_table(input)
# 
# ipd_data <- ipd_data %>%
#   rename(Age = X1) %>%
#   mutate(Y = abs(Y))
# 
# # new dataset for study 6 from study 5
# study6 <- ipd_data %>%
#   filter(Study == "5") %>%
#   filter(!(Age >= 20 & Age <= 60)) %>%  # remove Age between 20 and 60
#   mutate(Study = "6")
# ipd_data <- bind_rows(ipd_data, study6)
# 
# write_csv(ipd_data, "data/generated_ipd.csv")

input <- tribble(
  ~study, ~n, ~noise,     ~Y,         ~X1,
  "1", 20,    10, myfunc2, "10.0:20.0",
  "1", 20,    10, myfunc2, "20.0:30.0",
  
  "3",  5,     2,  myfunc,  "0.0:10.0",
  "3",  5,    10, myfunc2, "10.0:20.0",
  "3",  5,    10, myfunc2, "20.0:30.0",
  "3",  5,    10, myfunc3, "30.0:40.0",
  "3",  5,    10, myfunc3, "40.0:50.0",
  "3",  5,    10, myfunc3, "50.0:60.0",
  "3",  5,    20, myfunc3, "60.0:70.0",
  "3",  5,    20, myfunc3, "70.0:80.0",
  
  "5", 10,     2,  myfunc,  "0.0:10.0",
  "5", 10,    10, myfunc2, "10.0:20.0",
  "5", 10,    10, myfunc2, "20.0:30.0",
  "5", 10,    10, myfunc3, "30.0:40.0",
)

set.seed(2025)
ipd_data <- generate_from_table(input)

ipd_data <- ipd_data %>%
  rename(Age = X1) %>%
  mutate(Y = abs(Y))

# new dataset for study 6 from study 5
study2 <- ipd_data %>%
  filter(Study == "1") %>%
  filter(!(Age >= 15 & Age <= 25)) %>%  
  mutate(Study = "2")
study4 <- ipd_data %>%
  filter(Study == "3") %>%
  filter(!(Age >= 20 & Age <= 60)) %>%
  mutate(Study = "4")
study6 <- ipd_data %>%
  filter(Study == "5") %>%
  filter(!(Age >= 10 & Age <= 30)) %>%
  mutate(Study = "6")
ipd_data <- bind_rows(ipd_data, study2, study4, study6)

write_csv(ipd_data, "data/generated_ipd.csv")
