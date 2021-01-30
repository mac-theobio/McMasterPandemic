## Get contact matrices from Mistry et al. Github repo
## automatically (currently set up to grab Canada matrices)
##
## Author: Irena Papst

library(RCurl) ## to download raw csvs from the web
library(gh) ## github requests via REST api
library(rlist) ## working with lists
library(stringr) ## working with strings

## function to get and write the csvs
get_and_write_csv <- function(filename){
  the_url <- paste0("https://raw.githubusercontent.com/mobs-lab/mixing-patterns/main/data/contact_matrices/",
                    filename)
  data_string <- getURL(the_url)
  data_df <- read.csv(text = data_string, header = FALSE)
  ## remove default column names
  colnames(data_df) <- NULL
  write.csv(data_df,
            file = file.path("..", "inst", "params",
                             "mistry-cmats", filename),
            row.names = FALSE
  )
}

## get filenames from repo corresponding to Canada

## request tree structure of entire repo
req <- gh("GET /repos/mobs-lab/mixing-patterns/git/trees/d52409fb875c08a4a34004b87175afefeddec929?recursive=1")
## get tree elementsin the data/contact_matrices subdirectory
req <- list.filter(req$tree, str_detect(path, "data/contact_matrices/"))
## get filenames
filenames <- unlist(list.map(req, path))
## remove path from filenames
filenames <- sub("data/contact_matrices/", "", filenames)
## keep only files for Canada
filenames <- filenames[str_detect(filenames, "Canada")]

## get all matrices for Canada
for (filename in filenames){
  get_and_write_csv(filename)
}
