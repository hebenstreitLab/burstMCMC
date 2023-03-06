#' @export
doPriors <- function(priors) {
  distributions <- c('gamma', 'normal', 'uniform', 'dirac', 'log-uniform')
  priorProbabilities <- list()
  priorRandomVariables <- list()
  hyperParameters <- list()
  distributionIndices <- c(which(priors %in% distributions), length(priors) + 1)
  for (index in 1:(length(distributionIndices) - 1)) { #define how to treat inputs for a few possible prior distributions, can add more in later if needed
    distribution <- unname(unlist(priors[distributionIndices[index]]))
    hyperParameters[[index]] <- as.numeric(unname(unlist(priors[(distributionIndices[index] + 1):(distributionIndices[index + 1] - 1)])))
    if (distribution == 'gamma') {
      priorProbabilities[[index]] <- function(x, hyperParameters, index) {return(dgamma(x, shape = hyperParameters[[index]][1], scale = hyperParameters[[index]][2]))}
      priorRandomVariables[[index]] <- function(hyperParameters, index) {return(rgamma(1, shape = hyperParameters[[index]][1], scale = hyperParameters[[index]][2]))}
    } else if (distribution == 'normal') {
      priorProbabilities[[index]] <- function(x, hyperParameters, index) {return(dnorm(x, mean = hyperParameters[[index]][1], sd = hyperParameters[[index]][2]))}
      priorRandomVariables[[index]] <- function(hyperParameters, index) {return(rnorm(1, mean = hyperParameters[[index]][1], sd = hyperParameters[[index]][2]))}
    } else if (distribution == 'uniform') {
      priorProbabilities[[index]] <- function(x, hyperParameters, index) {return(dunif(x, min = hyperParameters[[index]][1], max = hyperParameters[[index]][2]))}
      priorRandomVariables[[index]] <- function(hyperParameters, index) {return(runif(1, min = hyperParameters[[index]][1], max = hyperParameters[[index]][2]))}
    } else if (distribution == 'dirac') {
      priorProbabilities[[index]] <- function(x, hyperParameters, index) {if (signif(x, 4) == signif(hyperParameters[[index]], 4)) {return(1)} else {return(0)}}
      priorRandomVariables[[index]] <- function(hyperParameters, index) {return(hyperParameters[[index]])}
    } else if (distribution == 'log-uniform') {
      priorProbabilities[[index]] <- function(x, hyperParameters, index) {return(dlunif(x, min = hyperParameters[[index]][1], max = hyperParameters[[index]][2]))}
      priorRandomVariables[[index]] <- function(hyperParameters, index) {return(rlunif(1, min = hyperParameters[[index]][1], max = hyperParameters[[index]][2]))}
    }
  }
  priorOutput <- list('hyperParameters' = hyperParameters, 'priorProbabilities' = priorProbabilities, 'priorRandomVariables' = priorRandomVariables)
  return(priorOutput)
}

#homemade function for log-uniform distribution density which returns p=0 if x is outside the range
#' @export
dlunif <- function(x, min, max) {
  return(1 / (x * log(max / min)) * as.numeric(min <= x & x <= max))
}

#and for sampling from log-uniform distribution via its cdf
#' @export
rlunif <- function(n, min, max) {
  quantiles <- runif(n)
  return(exp(quantiles * log(max / min) + log(min)))
}

#' @export
dmulti <- function(x, p) { #wrote own function for multinomial density which is superior to the in-built one
  if (any(p == 0)) { #exit with likelihood == 0 if we get a count of 1 or more on an outcome that has 0 probability
    return(-Inf)
  } else
    return(lgamma(sum(x) + 1) + sum(x * log(p) - lgamma(x + 1))) #this is the form of coding the multinomial pmf that is the most efficient and is the best at dealing with massive/tiny numbers
}

#' @export
initialiseChain <- function(priorRandomVariables, hyperParameters) {
  initialisation <- sapply(1:length(priorRandomVariables), function(index, priorRandomVariables, hyperParameters) {return(priorRandomVariables[[index]](hyperParameters, index))}, priorRandomVariables, hyperParameters)
  return(initialisation)
}

#' @export
getPython <- function() {
  # reticulate::import('scipy', delay_load = TRUE)
  # reticulate::import('numpy', delay_load = TRUE)
  # reticulate::import('math', delay_load = TRUE)
  reticulate::source_python(file = system.file('python/hypergeomCalculator.py', package = "burstMCMC", mustWork = T))
  reticulate::source_python(file = system.file('python/f_i_m.py', package = "burstMCMC", mustWork = T))
  calculateHypergeom <<- calculateHypergeom
  get_f_i_m <<- get_f_i_m
  return()
}

