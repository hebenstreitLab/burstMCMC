#' Simulate data in the format required for the inference algorithm based on real data and 'ground truth' parameter values to test algorithm performance on
#' @export
simulateGenes <- function(cores, Time, Data, estimates) {
  genes <- names(Data)
  simulatedData <- parallel::mclapply(genes, function(gene) {
    theta <- as.numeric(estimates[gene, ])
    mu <- theta[1]
    a <- min(c(1000, theta[2]))
    d <- 1 / theta[3]
    k <- a * d
    b <- mu / a
    p <- exp(- d * Time)
    geneData <- Data[[gene]]
    alphas <- geneData[["alphas"]]
    cells <- names(alphas)
    p_n <- geneData[["lambdaN"]]
    N <- length(cells)
    UMIcounts <- geneData[["UMIcounts"]]
    conversions <- geneData[["conversions"]]
    UsPerRead <- geneData[["UsPerRead"]]
    background <- geneData[['background']]
    
    readsPerUMI <- length(unlist(conversions)) / sum(unlist(UMIcounts))
    p_o <- background
    ztLambdas <- seq(0.001, readsPerUMI, 0.001)
    x <- ztLambdas / (1 - exp(-ztLambdas)) 
    lambda <- ztLambdas[which.min(abs(x - readsPerUMI))] #maximum likelihood estimate for lambda parameter (the zero-truncated poisson rate parameter for reads per UMI not the T>C rate parameters)
    UsPerReadCdf <- cumsum(UsPerRead)
    
    r_m <- rnbinom(N, size = a, mu = mu)
    r_o <- rbinom(N, r_m, p)
    
    bursts <- rpois(N, k * Time)
    r_n <- vector('numeric', N)
    for (index in 1:N) {
      Ts <- runif(bursts[index], 0, Time)
      sizes <- rgeom(bursts[index], 1 / (b + 1))
      r_n[index] <- sum(rbinom(bursts[index], sizes, exp(- d * Ts)))
    }
    r_ln <- rbinom(N, r_n, alphas)
    r_lo <- rbinom(N, r_o, alphas)
    
    newReadsByCell <- sapply(r_ln, function(l) {
      return(sum(actuar::rztpois(rep(1, l), lambda)))
    })
    oldReadsByCell <- sapply(r_lo, function(l) {
      return(sum(actuar::rztpois(rep(1, l), lambda)))
    })
    newUsByCell <- sapply(newReadsByCell, function(reads) {
      return(sapply(runif(reads), function(unif) {
        return(as.numeric(names(UsPerReadCdf))[which.max(unif <= UsPerReadCdf)])
      }))
    })
    oldUsByCell <- sapply(oldReadsByCell, function(reads) {
      return(sapply(runif(reads), function(unif) {
        return(as.numeric(names(UsPerReadCdf))[which.max(unif <= UsPerReadCdf)])
      }))
    })
    newCsByCell <- sapply(newUsByCell, function(U) {
      if (length(U) > 0) {
        return(rbinom(length(U), U, p_n + p_o))
      }
    })
    oldCsByCell <- sapply(oldUsByCell, function(U) {
      if (length(U) > 0) {
        return(rbinom(length(U), U, p_o))
      }
    })
    
    UMIcounts <- r_ln + r_lo
    conversions <- sapply(1:N, function(cell) {
      return(c(newCsByCell[[cell]], oldCsByCell[[cell]]))
    })
    names(UMIcounts) <- cells
    names(conversions) <- cells
    geneData <- list('UMIcounts' = unlist(UMIcounts), 'conversions' = conversions, 'UsPerRead' = UsPerRead, 'alphas' = alphas, 'background' = background, 'lambdaN' = p_n)
    return(geneData)
  }, mc.cores = cores)
  names(simulatedData) <- genes
  return(simulatedData)
}



