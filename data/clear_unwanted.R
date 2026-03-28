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
sub_dir <- c("simulated_many/", "simulated_some/", "simulated_no/")

my_dirs <- c(
  paste0("./data/", sub_dir),
  paste0("./data/extrapolate/", sub_dir),
  paste0("./data/extrapolate3/", sub_dir),
  paste0("./data/interpolate/", sub_dir),
  paste0("./data/interpolate1/", sub_dir),
  paste0("./data/interpolate2/", sub_dir),
  paste0("./data/interpolate3/", sub_dir),
  paste0("./data/simulate3/", sub_dir)
)

rank_n <- 3

for (my_dir in my_dirs) {

  cat("\nProcessing:", my_dir, "\n")

  # read settings
  df <- read.csv(paste0(my_dir, "simul_settings.csv"))

  # list files
  all_files <- list.files(my_dir, pattern = "^study.*\\.csv$", full.names = TRUE)

  same_knot_names <- tools::file_path_sans_ext(
    basename(all_files[
      sapply(basename(all_files), function(x) {
        digits <- sub("^study([0-9]+)_second[0-9]+\\.csv$", "\\1", x)
        length(unique(strsplit(digits, "")[[1]])) == 1
      })]
    ))
  
  selected_names <- sort(unique(c(
    df %>% arrange(Rank_pooled) %>% slice(1:rank_n) %>% pull(Name),
    df %>% arrange(desc(Rank_pooled)) %>% slice(1:rank_n) %>% pull(Name),
    df %>% arrange(Rank_meta) %>% slice(1:rank_n) %>% pull(Name),
    df %>% arrange(desc(Rank_meta)) %>% slice(1:rank_n) %>% pull(Name),
    same_knot_names
  )))
  
  # cat("Selected:", selected_names, "\n")
  cat("Selected:", length(selected_names), "\n")
    
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
