#' Setting up the package
#'
#' Setting up dartR.sim
#' @import adegenet
#' @importFrom dartR.base theme_dartR gl.check.verbosity gl.check.wd utils.flag.start gl.He gl.colors gl2gi utils.plot.save utils.check.datatype gl.filter.allna gl.compliance.check gl.keep.loc gl.alf utils.reset.flags
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

# SET GLOBAL VARIABLES FOR SIMULATION FUNCTIONS
# for gl.sim.WF.table
utils::globalVariables(c("q_neutral","chromosome_name","chunk_number","real_loc","chunk_recombination","map_resolution","gamma_scale","gamma_shape","log_mean","log_sd","rate","exp_rate","chunk_cM","loci_mutation","mutations_factor", "chunk_neutral_loci", "deleterious_factor", "h_adv", "h_del", "h_distribution_adv", "h_distribution_del", "h_intercept_adv", "h_intercept_del","h_mean_adv", "h_mean_del","h_rate_adv","h_rate_del","h_sd_adv","h_sd_del","loci_advantageous","loci_deleterious", "loci_mut_adv", "loci_mut_del","loci_mut_neu","q_adv", "q_del", "q_distribution_adv","q_distribution_del", "q_equation_adv", "q_equation_del","s_adv","s_del", "s_distribution_adv","s_distribution_del"))
# for gl.sim.WF.run
utils::globalVariables(c("chromosome_name","phase1","same_line","number_pops_phase1","population_size_phase1","gen_number_phase1","dispersal_phase1","dispersal_type_phase1","number_transfers_phase1","transfer_each_gen_phase1","variance_offspring_phase1","number_offspring_phase1","selection_phase1","Ne_phase1","Ne_fst_phase1","number_pops_phase2","population_size_phase2","gen_number_phase2","dispersal_phase2","dispersal_type_phase2","number_transfers_phase2","transfer_each_gen_phase2","variance_offspring_phase2","number_offspring_phase2","selection_phase2","Ne_phase2","Ne_fst_phase2","real_freq","real_pop_size","real_pops","recombination","recombination_males","genetic_load","natural_selection_model","mutation","mut_rate","dispersal_rate_phase1","dispersal_rate_phase2","clinal_strength", "clinal_adap","local_adap"))

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


