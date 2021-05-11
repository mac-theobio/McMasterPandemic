library("McMasterPandemic")
library("stringr")
library("tidyverse")
library('tools')

care_these_files <- function(x) {
  is_R <- (tools::file_ext(x) %in% c('R', 'Rnw', 'rmd', 'Rmd'))
  good_dir <- (str_detect(dirname(x), regex('(tests)|(vignettes)|R')))
  return(is_R && good_dir)
}

files <- list.files('../..', recursive=TRUE)
interesting <- sapply(files, care_these_files)

files <- files[interesting]
df <- data.frame(character(), numeric(), numeric())


for (file in files) {
  query <- "(anytime)|(anydate)|20\\d\\d" #
  readin <- read_lines(paste('../../',file, sep = ""))
  counted <- str_count(readin, regex(query))
  instances <- sum(counted)
  line_nums <- toString(c(1:(length(counted)-1))[counted > 0])


  df <- rbind(df, c(file_name = paste('../../', file, sep=''), instances=instances, line_nums=line_nums))
}

names(df) <- c('file_name', 'instances', 'line_nums')

print(nrow(df[df$instances != 0,]))

print(sum(as.numeric(df$instances)))
write.csv(df[df$instances != 0,], 'potentially_problematic.csv')



care_these_files <- function(x) {
  is_R <- (tools::file_ext(x) %in% c('R'))
  return(is_R)
}

files <- list.files('../..', recursive=TRUE)
interesting <- sapply(files, care_these_files)

files <- files[interesting]
df <- data.frame(character(), numeric(), numeric())


for (file in files) {
  query <- "ont_all_sub" #
  readin <- read_lines(paste('../../',file, sep = ""))
  counted <- str_count(readin, regex(query))
  instances <- sum(counted)
  line_nums <- toString(c(1:(length(counted)-1))[counted > 0])


  df <- rbind(df, c(file_name = paste('../../', file, sep=''), instances=instances, line_nums=line_nums))
}

names(df) <- c('file_name', 'instances', 'line_nums')

print(nrow(df[df$instances != 0,]))

print(sum(as.numeric(df$instances)))
write.csv(df[df$instances != 0,], 'ont_all_sub.csv')
