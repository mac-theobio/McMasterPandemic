
## Utilities
targetname <- function(fl = commandArgs(TRUE)){
	return(sub("\\.Rout$", "", fl[[1]]))
}

fileSelect <- function(fl = commandArgs(TRUE), exts)
{
	outl <- character(0)
	for (ext in exts){
		if(grepl("\\.", ext))
			warning("Extension", ext, "starts with . in fileSelect")
		ss <- paste0("\\.", ext, "$")
		outl <- c(outl, grep(ss, fl, value=TRUE))
	}
	return(outl)
}

### Loading and reading

## This is meant to be a default starting point for $(makeR) scripts
## wrapmake encodes the current defaults for $(run-R) scripts
commandFiles <- function(fl = commandArgs(TRUE)){
	commandEnvironments(fl)
	commandEnvirLists(fl)
	commandLists(fl)
	sourceFiles(fl, first=FALSE)
}

## Source certain files from a file list
sourceFiles <- function(fl=commandArgs(TRUE) 
	, exts=c("R", "r"), first=TRUE)
{
	fl <- fileSelect(fl, exts)
	if (!first) fl <- fl[-1]
	for (f in fl){
		source(f)
	}
}

## Read environments from a file list to a single environment
commandEnvironments <- function(fl = commandArgs(TRUE)
	, exts = c("RData", "rda"), parent=.GlobalEnv
)
{
	envl <- fileSelect(fl, exts)
	loadEnvironments(envl, parent)
	invisible(envl)
}

## Read rds lists from a file list to a single environment
commandLists <- function(fl = commandArgs(TRUE)
	, exts = c("Rds", "rds"), parent=.GlobalEnv
)
{
	varl <- fileSelect(fl, exts)
	loadVarLists(varl, parent)
	invisible(varl)
}

## Wrapper for legacy makefiles
## By default takes Rout dependencies and assumes rda environments
legacyEnvironments <- function(fl = commandArgs(TRUE)
	, dep = "Rout", ext="rda")
{
	envl <- fileSelect(fl, dep)
	if(length(envl>1)){
		ss <- paste0(dep, "$")
		envl <- sub(ss, ext, envl[-1])
	}
	loadEnvironments(envl)
	invisible(envl)
}

## Read environments from a file list to separate places
## NOT implemented
commandEnvirLists <- function(fl = commandArgs(TRUE)
	, exts = c("RData", "Rdata", "rdata", "rda")
)
{
	invisible(0)
}

## Load every environment found into GlobalEnv
## This is the simple-minded default
loadEnvironments <- function(envl, parent=.GlobalEnv)
{
	for (env in envl){
		load(env, parent)
	}
}

## Load every list found into GlobalEnv
## This is the efficient rds analogue of the simple-minded default
loadVarLists <- function(varl, parent=parent.frame())
{
	for (v in varl){
		l <- readRDS(v)
    	list2env(l, envir=parent)
	}
}

#### Graphics

makeGraphics <- function(...
	, target = commandArgs(TRUE)[[1]]
	, otype = NULL, ext = otype
)
{
	if(is.null(ext)) ext = "pdf.tmp"
	if(is.null(otype)) otype = "pdf"
	fn <- paste0(target, ".", ext)
	get(otype)(..., file=fn)
}

#### Saving

saveEnvironment <- function(target = targetname(), ext="rda"){
	save.image(file=paste(target, ext, sep="."))
}

saveVars <- function(..., target = targetname(), ext="rdata"){
	save(file=paste(target, ext, sep="."), ...)
}

## FIXME: I have the wrong environment for objects
saveList <-  function(..., target = targetname(), ext="rds"){
	l <- list(...)
	if(length(l)==0){
		names <- objects(parent.frame())
	} else {
		names <- as.character(substitute(list(...)))[-1]
	}

	outl <- list()
	for (n in names){
		outl[[n]] <- get(n)
	}
	saveRDS(outl, file=paste(target, ext, sep="."))
	return(invisible(names(outl)))
}
