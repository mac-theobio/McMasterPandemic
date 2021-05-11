# devtools::install_github("datastorm-open/DependenciesGraphs")
library("DependenciesGraphs")
library("mvbutils")
library("McMasterPandemic")
library("stringr")
library("tidyverse")
library("nos")

pos <- which(search() == "package:McMasterPandemic")
## FIXME: include predict methods?
funs <- grep(".*", ls(pos = pos), value = TRUE) # mle_fun|forecast|calibrate$|run_|do_step
funs <- c(funs, "predict.fit_pansim")
ff <- foodweb(
  where = pos,
  funs = funs,
  ## rprune="forecast|calibrate|run_",
  plotting = FALSE, ancestors = FALSE, descendents = FALSE
)
## HACK: foodweb doesn't recognize do.call() or ... ?
M <- ff$funmat
M["predict.fit_pansim", c("forecast_sim", "forecast_ensemble")] <- 1
## HACK: calibrate effectively calls the other run_sim
run_sim_funs <- setdiff(grep("run_sim_", rownames(M), value = TRUE), "run_sim_range")
M["forecast_sim", run_sim_funs] <- 1
M[run_sim_funs, "run_sim"] <- 1

edges <- nos::freqMat_2_edge(M)
edges <- data.frame(from = edges[, 1], to = edges[, 2])
Nomfun <- data.frame(id = c(1:(length(colnames(M)))), label = c(colnames(M)))
new_dep <- list(Nomfun = Nomfun, fromto = edges)
new_dep <- structure(new_dep, class = "dependenciesGraphs")
not_orphan <- Reduce("|", list((new_dep$Nomfun$id %in% new_dep$fromto$from), (new_dep$Nomfun$id %in% new_dep$fromto$to)))
new_dep$Nomfun <- new_dep$Nomfun[not_orphan, ]

widget <- plot(new_dep, height = "1600px", width = "100%")

htmlwidgets::saveWidget(widget, "dependency_graph.html")


# dep <- DependenciesGraphs::envirDependencies("package:McMasterPandemic")
# functions = c(dep$Nomfun$label)
# use_this_fn <- str_detect(functions, "(mle_fun)|(forecast)|(calibrate$)|(run_)|(do_step)|(predict.fit_pansim)")
# not_orphan <- Reduce('|', list((dep$Nomfun$id %in% dep$fromto$from),(dep$Nomfun$id %in% dep$fromto$to)))
# dep$Nomfun = dep$Nomfun[not_orphan,]
# keep_id <- dep$Nomfun$id[use_this_fn]
# keep_rows <- Reduce('&', list((dep$fromto$from %in% keep_id), (dep$fromto$to %in% keep_id)))
# dep$fromto = dep$fromto[keep_rows, ]
# dep$Nomfun = dep$Nomfun[use_this_fn,]
