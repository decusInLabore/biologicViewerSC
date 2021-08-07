## code to prepare `dfkey` dataset goes here

FNkey <- "data-raw/connect/db.txt"

if (file.exists(FNkey)){
  dfkey <- read.delim(
    FNkey, 
    sep = "\t",
    stringsAsFactors = F
  )
}

usethis::use_data(dfkey, overwrite = TRUE)
