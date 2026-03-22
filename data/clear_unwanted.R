# true.mean = 20
# true.sd = 5
# N = 1000
# cover = rep(NA, N)
# 
# for (i in 1:N) {
#   ddddd <- rnorm(100, true.mean, true.sd)
#   uT <- mean(ddddd) + 1.96 * sd(ddddd)/sqrt(100)
#   lT <- mean(ddddd) - 1.96 * sd(ddddd)/sqrt(100)
#   cover[i] = lT < true.mean && true.mean < uT
# }
# 
# mean(cover)  # should be close to 0.95

library(dplyr)

# directories to process
my_dirs <- c(
  # "./data/simulated_many/",
  # "./data/simulated_some/",
  # "./data/simulated_no/",
  # "./data/extrapolate/simulated_many/",
  # "./data/extrapolate/simulated_some/",
  "./data/interpolate/simulated_no/",
  "./data/interpolate/simulated_many/",
  "./data/interpolate/simulated_some/"
)

rank_n <- 12

for (my_dir in my_dirs) {

  cat("\nProcessing:", my_dir, "\n")

  # read settings
  df <- read.csv(paste0(my_dir, "simul_settings.csv"))

  # select names
  selected_names <- sort(unique(c(
    df %>% arrange(Rank_pooled) %>% slice(1:rank_n) %>% pull(Name),
    df %>% arrange(desc(Rank_pooled)) %>% slice(1:rank_n) %>% pull(Name),
    df %>% arrange(Rank_meta) %>% slice(1:rank_n) %>% pull(Name),
    df %>% arrange(desc(Rank_meta)) %>% slice(1:rank_n) %>% pull(Name)
  )))

  # print(selected_names)
  cat("Selected:", length(selected_names), "\n")

  # list files
  all_files <- list.files(my_dir, pattern = "\\.csv$", full.names = TRUE)
  all_files <- all_files[basename(all_files) != "simul_settings.csv"]
  all_files <- all_files[basename(all_files) != "generated_ipd.csv"]

  # match filenames containing selected names
  pattern <- paste(selected_names, collapse="|")
  keep <- grepl(pattern, all_files)

  # files to delete
  files_to_delete <- all_files[!keep]

  # preview
  # print(basename(files_to_delete))

  # delete
  file.remove(files_to_delete)

  cat("Deleted", length(files_to_delete), "files\n")
}
