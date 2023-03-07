#' Install miniconda python environment and required modules if not already installed and then load python functions
#' #' @export
getPython <- function() {
  try(reticulate::install_miniconda(), silent = T)
  if (!reticulate::py_module_available('scipy')) {
    reticulate::py_install('scipy')
  }
  if (!reticulate::py_module_available('numpy')) {
    reticulate::py_install('numpy')
  }
  reticulate::source_python(file = system.file('python/hypergeomCalculator.py', package = "burstMCMC", mustWork = T), )
  reticulate::source_python(file = system.file('python/f_i_m.py', package = "burstMCMC", mustWork = T))
  calculateHypergeom <<- calculateHypergeom
  get_f_i_m <<- get_f_i_m
  return()
}