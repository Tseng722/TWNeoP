library(magrittr)
library(data.table)
library(antigen.garnish)

command_line_args <- commandArgs(trailingOnly = TRUE)

if (length(command_line_args) < 1) {
  cat("Usage: Rscript myscript.R <my_local_variable>\n")
  quit(status = 1)
}

file_name <- command_line_args[1]


# file_path <- "/root/local179/myfilecopy.txt"
file_path <- paste0("/root/local179/", file_name, ".txt")

v <- readLines(file_path)

v %>% foreignness_score(db = "human") %>% print

v %>% dissimilarity_score(db = "human") %>% print


result_fs <- v %>% foreignness_score(db = "human")
output_filef <- paste0("/root/local179/f_", file_name, ".txt")
write.table(result_fs, file = output_filef, sep = "\t", quote = FALSE, row.names = FALSE) 

result_ds <- v %>% dissimilarity_score(db = "human")
# output_filef <- "/root/local179/outputf.txt"
# output_filed <- "/root/local179/outputd.txt"
output_filed <- paste0("/root/local179/d_", file_name, ".txt")
write.table(result_ds, file = output_filed, sep = "\t", quote = FALSE, row.names = FALSE) 
