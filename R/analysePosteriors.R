#' @export
getModels <- function(outputs) {
  models <- sapply(genes, function(gene) {
    return(outputs[[gene]][['model']])
  })
  return(models)
}
#' @export
getPosteriors3P <- function(outputs) {
  posteriors <- sapply(genes, function(gene) {
    posterior <- outputs[[gene]][['posterior']]
    colnames(posterior) <- c('Expression level (transcripts / cell)', 'Burst rate (bursts / lifetime)', 'Transcript lifetime (minutes)')
    return(posterior)
  }, simplify = F)
  return(posteriors)
}
#' @export
getPosteriors5P <- function(outputs) {
  posteriors <- sapply(genes, function(gene) {
    return(outputs[[gene]][['posterior']])
  }, simplify = F)
  posteriors <- sapply(posteriors, function(posterior) {
    posterior <- cbind(posterior[, 1:2], 1 / posterior[, 3], posterior[, 1] / posterior[, 2], posterior[, 2] / posterior[, 3])
    colnames(posterior) <- c('Expression level (transcripts / cell)', 'Burst rate (bursts / lifetime)', 'Decay rate (1 / minutes)', 'Burst size (transcripts / burst)', 'Burst frequency (1 / minutes)')
    return(posterior)
  }, simplify = F)
  return(posteriors)
}
#' @export
plotTraces <- function(gene, posteriors) {
  posterior <- posteriors[[gene]]
  if (ncol(posterior) == 3) {
    par(mfrow = c(2, 2))
  } else if (ncol(posterior) == 5) {
    par(mfrow = c(3, 2))
  }
  for (parameter in colnames(posterior)) {
    plot(posterior[, parameter], type = 'l', xlab = 'Chain step', ylab = 'Parameter value', main = parameter)
  }
}
#' @export
plotPosteriors <- function(gene, posteriors, thinning) {
  posterior <- posteriors[[gene]]
  if (ncol(posterior) == 3) {
    par(mfrow = c(2, 2))
  } else if (ncol(posterior) == 5) {
    par(mfrow = c(3, 2))
  }
  if (nrow(posterior) == 1500) {
    stepIndices <- seq(501, 1500, thinning)
  } else if (nrow(posterior) == 5000) {
    stepIndices <- seq(2501, 5000, thinning)
  }
  for (parameter in colnames(posterior)) {
    plot(density(posterior[stepIndices, parameter]), ylab = 'Probability density', xlab = 'Parameter value', main = parameter, lwd = 2)
  }
}
#' @export
checkAcceptanceRates <- function(posteriors) {
  A <- sapply(1:3, function(i) {
    return(sapply(posteriors, function(posterior) {
      if (nrow(posterior) == 1500) {
        return(length(unique(posterior[501:1500, i])) / 1000)
      } else if (nrow(posterior) == 5000) {
        return(length(unique(posterior[2501:5000, i])) / 2500)
      }
    }))
  })
  parameters <- colnames(posteriors[[1]])
  par(mfrow = c(2, 2))
  for (i in 1:3) {
    plot(density(A[, i]), main = parameters[i], xlab = 'Acceptance rate', lwd = 2)
    abline(v = 0.574, col = 'blue', lwd = 2, lty = 'dashed')
  }
}
#' @export
getEstimates <- function(posteriors, thinning) {
  means <- t(sapply(posteriors, function(posterior) {
    if (nrow(posterior) == 1500) {
      return(colMeans(posterior[seq(501, 1500, thinning), ]))
    } else if (nrow(posterior) == 5000) {
      return(colMeans(posterior[seq(2501, 5000, thinning), ]))
    }
  }))
  SDs <- t(sapply(posteriors, function(posterior) {
    if (nrow(posterior) == 1500) {
      return(apply(posterior[seq(501, 1500, thinning), ], 2, sd))
    } else if (nrow(posterior) == 5000) {
      return(apply(posterior[seq(2501, 5000, thinning), ], 2, sd))
    }
  }))
  CVs <- SDs / means
  return(list('means' = means, 'SDs' = SDs, 'CVs' = CVs))
}






