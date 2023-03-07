#' @export
getPython <- function() {
  if (class(try(reticulate::install_miniconda(), silent = T)) != 'try-error') {
    reticulate::install_miniconda()
    reticulate::py_install('scipy')
    reticulate::py_install('numpy')
  }
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