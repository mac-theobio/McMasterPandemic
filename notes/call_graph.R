library(McMasterPandemic)
library(mvbutils)
library(visNetwork)
library(colorspace)

## specify function names to include in call graph

make_callmat <- function(include_patterns = c("mle_fun|forecast|calibrate|run_|do_step"),
                         include_funs = "predict.fit_pansim",
                         add_links = list(list("predict.fit_pansim", c("forecast_sim", "forecast_ensemble")))) {
  pos <- which(search()=="package:McMasterPandemic")
  funs <- unlist(lapply(include_patterns, grep, x = ls(pos = pos), value = TRUE))
  funs <- c(funs, include_funs)
  ff <- foodweb(where=pos,
                funs=funs,
                ## rprune="forecast|calibrate|run_",
                plotting=FALSE, ancestors=FALSE,descendents=FALSE)
  ## HACK: foodweb doesn't recognize do.call() or ... ?
  M <- ff$funmat
  for (a in add_links) {
    M[a[[1]], a[[2]]] <- 1
  }
  ## HACK: calibrate effectively calls the other run_sim
  run_sim_funs <- setdiff(grep("run_sim_",rownames(M),value=TRUE), "run_sim_range")
  M["forecast_sim",run_sim_funs] <- 1
  M[run_sim_funs,"run_sim"] <- 1
  return(M)
}

vis_callmat <- function(M, files = NULL, colour = NULL) {
  e <- which(M==1, arr.ind=TRUE)
  e[] <- rownames(M)[e]
  e <- setNames(as.data.frame(e), c("from", "to"))
  dat <- data.frame(id=rownames(M), label=rownames(M))
  if (!is.null(files)) {
    if (is.null(colour)) {
      nm <- unique(files)
      colour <- colorspace::qualitative_hcl(length(nm))
      names(colour) <- nm
    }
    dat$color <- colour[files]
    dat$label <- files
  }
  v <- visNetwork(nodes = dat,
             edges=e) %>%
    visEdges(arrows = "to") %>%
    visHierarchicalLayout(sortMethod = "directed")
  if (!is.null(files)) {
    v <- v %>%
      visLegend(
          position = "right",
          addNodes = unique(dat[c("color", "label")]))
  }
  return(v)
}

##https://datastorm-open.github.io/visNetwork/layout.html


## TO DO; annotate by file location (manually)
## separate graphs by task (calibrate, forecast) ?

find_fun <- function(f) {
  f2 <- system(sprintf("grep -l '%s <- ' R/*.R", f), intern = TRUE)
  if (length(f2) != 1) stop("ambiguous")
  return(sub("R/", "", f2))
}

M <- make_callmat()

files <- vapply(rownames(M), find_fun, FUN.VALUE = character(1))

vis_callmat(M, files = files)

find_fun <- function(f) {
  f2 <- system(sprintf("grep -l '%s <- ' R/*.R", f), intern = TRUE)
  if (length(f2) != 1) stop("ambiguous")
  return(sub("R/", "", f2))
}


