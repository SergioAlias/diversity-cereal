processAncomCsv <- function(file_path) {
  df <- read_csv(file_path) %>% 
    select(-matches("\\(Intercept\\)")) %>%
    rename_with(~ if_else(.x == "id", .x, paste0(.x,
                                                 "_",
                                                 sub("_slice\\.csv$",
                                                     "",
                                                     basename(file_path)))))
  return(df)
}

import_ancombc <- function(qza_path) {
  csv_name <- unzip(qza_path, list = TRUE, exdir = tempdir()) %>%
    filter(str_detect(.data$Name, ".csv")) %>% {.$Name}
  unzip(qza_path, files = csv_name, exdir = tempdir())
  csv_files <- setNames(lapply(file.path(tempdir(), csv_name), processAncomCsv),
                        sub("_slice\\.csv$", "", basename(csv_name)))
  merged_df <- reduce(csv_files, full_join, by = "id")
  return(merged_df)
}