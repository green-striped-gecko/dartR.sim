#' Setting up the package
#'
#' Setting up dartR.sim
#' @rawNamespace import(adegenet, except = c(glMean, glSum))
#' @importFrom dartR.base theme_dartR gl.check.verbosity gl.check.wd utils.flag.start gl.colors gl2gi utils.plot.save utils.check.datatype gl.filter.allna gl.compliance.check gl.keep.loc gl.alf utils.reset.flags utils.history.call
#' @import dartR.data
#' @import ggplot2
#' @importFrom utils packageVersion head read.csv tail
#' @importFrom methods new
#' @importFrom data.table rbindlist
#' @importFrom methods getPackageName
#' @importFrom utils packageVersion
#' @importFrom grDevices hcl
#' @keywords internal

#needed to avoid error
zzz<-NULL

error <- crayon::red
warn <- crayon::yellow
report <- crayon::green
important <- crayon::blue
code <- crayon::cyan


# WELCOME MESSAGE
.onAttach <- function(...) {
  pn <- getPackageName()
  packageStartupMessage(important(
    paste(
      "**** Welcome to",pn,"[Version",
      packageVersion(pn),
      "] ****\n"
    )
  ))
}


# ## returns NULL if the 'fbm' slot is missing OR is NULL
# .fbm_or_null <- function(x) {
#   if (methods::.hasSlot(x, "fbm")) {
#     val <- methods::slot(x, "fbm")
#     return(if (is.null(val)) NULL else val)
#   }
#   NULL
# }
# 
# .onLoad <- function(libname, pkgname) {
#   # Only set a default if user hasn’t set it already
#   if (is.null(getOption("dartR_fbm"))) {
#     val <- Sys.getenv("dartR_fbm", "")
#     # Accept a few truthy values: 1, true, yes, on (case-insensitive)
#     if (val=="TRUE") options(dartR_fbm = TRUE) else options(dartR_fbm=FALSE)
#   }
# }


