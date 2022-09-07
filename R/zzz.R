.onLoad <- function(lib, pkg) {
  factory_fresh_macpan_options()
}

.onUnload <- function(libpath) {
    library.dynam.unload("McMasterPandemic", libpath)
}
