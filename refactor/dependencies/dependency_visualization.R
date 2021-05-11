# devtools::install_github("datastorm-open/DependenciesGraphs")
library("DependenciesGraphs")
library("mvbutils")
library("McMasterPandemic")
library("stringr")
library("tidyverse")
dep <- DependenciesGraphs::envirDependencies("package:McMasterPandemic")
widget <- plot(dep, height = "1600px", width = "100%")

htmlwidgets::saveWidget(widget, "dependency_graph.html")

funs <- find.funs("package:McMasterPandemic")

files <- list.files("../../R/")

df <- data.frame(character(), character())
names(df) <- c('fun_name', 'file_name')

# Hacky way to locate every function in the package
for (fun in funs) {
  detected <- FALSE
  for (file in files) {
    # function_name [whitespace] <- [whitespace] function [whitespace] 
    query <- paste(fun, "\\s*<-\\s*function\\s*", sep = "")
    readin <- read.delim(paste("../../R/", file, sep = ""))
    
    # If we detect a function in a file,
    if (str_detect(readin, regex(query))) {
      # If it's the only instance, it's probably correct
      if (!detected){
        df <- rbind(df, c(fun_name=fun, file_name = file))
        detected <- TRUE}
      # If not, requires manual checking
      else{
        df$file_name[nrow(df)] <- paste(df$file_name[nrow(df)], file) 
      }
    }

  }
  if (!detected) {
    df <- rbind(df, c(fun_name=fun, file_name = 'Unknown'))
  }
}

df <- df[order(df$file_name, df$fun_name),]
write.csv(df, 'function_locations.csv', row.names=FALSE)