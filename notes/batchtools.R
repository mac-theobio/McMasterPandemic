source("makestuff/makeRfuns.R")
commandEnvironments()

get_workers <- function(workers) {
    if (is.null(workers)) {
        if (file.exists("MAXCORES")) {
            workers <- scan("MAXCORES",what=integer(),quiet=TRUE)
        } else if (nzchar(workers <- Sys.getenv("MAXCORES"))) {
            workers <- as.numeric(workers)
        } else {
            workers <- parallel::detectCores()-1
        }
    }
    return(workers)
}
    

##' @param workers max number of workers (cores etc.) to use
##' @param walltime max time per job (in seconds)
##' @param ncpus cpus per task (== processes per node, ppn on an SGE scheduler)
batch_setup <- function(workers=NULL,walltime=6*60*60,ncpus=1) {
    require("batchtools")
    require("furrr")
    require("future")
    require("future.batchtools")
    hostname  <- system("hostname", intern=TRUE)
    if (grepl("^gra-login",hostname)) hostname <- "graham"
    workers <- get_workers(workers)
    switch(hostname,
           ## jdserv/yushan: SGE scheduler
           jdserv = {
               library(future.batchtools)
               ## get default template
               tmpl <- system.file("templates","sge-simple.tmpl",package="batchtools")	   
               fn <- "batchtools.sge.tmpl"
               if (!file.exists(fn)) {
                   ## minimal modifications to work on jdserv/yushan:
                   ##   (1) execute in current working directory, (2) use queue 'all.q'
                   r <- readLines(tmpl)
                   ## add line for processes per node
                   queue_line <- grep("-q ",r)
                   r <- c(r[1:queue_line],
                          "## allow more than one process per node",
                          "#$ -l nodes=1,ppn=<%= resources$ncpus %>,tpp=<%= resources$ncpus %>",
                          "",
                          r[(queue_line+1):length(r)]
                          )
                   ## below is OBSOLETE/wrong; single comments in template don't mean "ignore"!
                   ##  queue identity is specified in resources()
                   ## cwdline <- grep("-cwd",r)
                   ## r[cwdline] <- "$ -cwd"
                   ## qline <- grep("-q ",r)
                   ## r[qline] <- "$ -q all.q"
                   writeLines(r, con=fn)
               }
               plan(batchtools_sge,resources=list(queue="all.q"))
           },
           ## SHARCnet
           graham = {
               options(future.wait.interval = 10.0)  ## needed to avoid slurm timeouts: see
                                                     ## https://github.com/mllg/batchtools/issues/201
               tmpl <- system.file("templates","slurm-simple.tmpl",package="batchtools")	   
               fn <- "batchtools.slurm.tmpl"
               if (!file.exists(fn)) {
                   r <- readLines(tmpl)
                   ## in slurm syntax, a single hash is NOT a comment! we don't need to/shouldn't
                   ##   uncomment this line!
                   ## nwtline <- grep("walltime",r)
                   ## r[wtline] <- gsub("^#","",r[wtline])
                   writeLines(r, con=fn)
               }
               plan(batchtools_slurm,resources=list(walltime=walltime,ncpus=ncpus))
           },
           ## default: earnserv, local, or any other machine without
           ## a scheduling system
           {
               plan(strategy=multiprocess,workers=workers)
           }
           ) ## end switch
    return(NULL)
}


saveEnvironment()
