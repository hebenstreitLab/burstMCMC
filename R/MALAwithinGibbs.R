#' Function for sampling from posterior distributions of bursting parameters over a set of genes in parallel
#' @param model Which model to use to carry out inference (1, 2 or 3)
#' @param cores The number of cores to parallelise over
#' @param Time The 4sU pulse duration in minutes
#' @param Data The dataset, in a format consistent with the output of the 'formatData' function
#' @return The model used (tracked in case forced to switch from model 2 to 3) and the sample from the posterior distribution generated by the Markov chain
#' @export
MALAwithinGibbs <- function(model, cores, Time, Data) {
  dimension <- 3
  J <- dimension - as.numeric(model == 1)
  genes <- names(Data)
  
  outputs <- parallel::mclapply(genes, function(gene) {
    geneData <- Data[[gene]]
    alphas <- geneData[["alphas"]]
    UMIs <- geneData[["UMIcounts"]]
    noReads <- UMIs == 0
    noReadIndices <- which(noReads)
    readIndices <- which(!noReads)
    n <- length(alphas)
    ME <- mean(UMIs) / mean(alphas) #simple estimate for mean expression level
    TL <- Time
    
    #set dummy variables when only model 1 is used and T>C / U data is not needed
    P_u <- 0
    u <- 0
    I <- 0
    Y <- 0
    TCs <- 0
    Ys <- 0
    Is <- 0
    p_o <- 0
    p_n <- 0
    
    if (model != 1) {
      conversions <- geneData[["conversions"]]
      P_u <- geneData[['UsPerRead']]
      p_o <- geneData[['background']]
      p_n <- geneData[['lambdaN']]
      u <- as.numeric(names(P_u))
      I <- max(unlist(conversions))
      Y <- vector('numeric', I + 1)
      TCs <- table(unname(unlist(conversions)))
      Y[as.numeric(names(TCs)) + 1] <- TCs
      Ys <- sapply(1:length(readIndices), function(index) {
        y <- table(conversions[[readIndices[index]]])
        return(y)
      }, simplify = F)
      Is <- sapply(Ys, function(y) as.numeric(names(y)), simplify = F)
      Ys <- sapply(Ys, function(y) as.numeric(unname(y)), simplify = F)
      lambda <- sum(0:I * Y / sum(Y)) / sum(P_u * u)
      TL <- -Time / log(max(c(0.1, min(c(0.9, (1 - ((lambda - p_o) / p_n))))))) #estimate for transcript lifetime
    }
    steps <- 5000
    window <- 500
    if (ME >= 1000) { #only run for shorter time for the super high expression genes since they take longer and dont need to be run for long to compute expectations
      steps <- 1500
      window <- 100 #make high expression genes more sensitive to f_n failures to make switching quicker to avoid wasting time
    }
    #use MALA approach with 1d scaling of the time step also to adaptively make the acceptance rates converge to 0.574 for each parameter
    i <- 0
    while (i < steps) { #while the chain is not completed, retain the ability to restart a model 2 chain using model 3 when f_n failure happens
      i <- 0
      posterior <- matrix(ncol = dimension, nrow = steps)
      f_n_failures <- vector('logical', steps)
      likelihoods <- vector('numeric', steps)
      gradients <- matrix(ncol = dimension, nrow = steps)
      A <- rep(0.574, dimension)
      S <- rep(0.1, dimension) #time step
      N <- rep(0.1, dimension)
      priors <- c('dirac', ME, 'log-uniform', 1, 10, 'normal', TL, TL / 5) #expression level, burst rate, transcript lifetime
      priorOutput <- doPriors(priors)
      hyperParameters <- priorOutput[['hyperParameters']]
      priorProbabilities <- priorOutput[['priorProbabilities']]
      priorRandomVariables <- priorOutput[['priorRandomVariables']]
      likelihood <- -Inf
      gradient <- rep(0, dimension)
      while(likelihood == -Inf) {
        proposal <- initialiseChain(priorRandomVariables, hyperParameters)
        theta <- proposal
        jumpEvaluation <- estimateLikelihood(Time, proposal, theta, likelihood, model, gradient, i, S, alphas, readIndices, UMIs, P_u, u, I, Y, TCs, Ys, Is, p_o, p_n, n, hyperParameters, priorProbabilities, priorRandomVariables)
        likelihood <- jumpEvaluation[['likelihood']]
      }
      theta <- jumpEvaluation[['theta']]
      gradient <- jumpEvaluation[["gradient"]]
      posterior[1, ] <- theta
      likelihoods[1] <- likelihood
      gradients[1, ] <- gradient
      priors <- c('uniform', 0, 100000, 'uniform', 0, 100000, 'uniform', 1, 100000) #can't have TL < 1 minute ish cause it makes tau too high in the hypergeometric function calculation
      priorOutput <- doPriors(priors)
      hyperParameters <- priorOutput[['hyperParameters']]
      priorProbabilities <- priorOutput[['priorProbabilities']]
      priorRandomVariables <- priorOutput[['priorRandomVariables']]
      S <- theta / 100
      A <- rep(0, dimension)
      for (j in 1:J) {
        while (A[j] != 1) {
          proposal <- theta
          proposal[j] <- theta[j] + sqrt(2 * S[j]) * rnorm(1)
          subGradient <- vector('numeric', dimension)
          subGradient[j] <- gradient[j]
          jumpEvaluation <- estimateLikelihood(Time, proposal, theta, likelihood, model, subGradient, i, S, alphas, readIndices, UMIs, P_u, u, I, Y, TCs, Ys, Is, p_o, p_n, n, hyperParameters, priorProbabilities, priorRandomVariables, j)
          A[j] <- jumpEvaluation[['A']]
        }
        likelihood <- jumpEvaluation[['likelihood']]
        theta[j] <- jumpEvaluation[["theta"]][j]
        gradient[j] <- jumpEvaluation[["gradient"]][j]
      }
      posterior[2, ] <- theta
      likelihoods[2] <- likelihood
      gradients[2, ] <- gradient
      for (i in 3:steps) {
        # print(i)
        # print(theta)
        for (j in 1:J) {
          proposal <- theta
          proposalGeneration <- generateProposal(theta[j], A[j], S[j], N[j], gradient[j])
          proposal[j] <- proposalGeneration[["proposal"]]
          S[j] <- proposalGeneration[["S"]]
          N[j] <- proposalGeneration[["N"]]
          subGradient <- vector('numeric', dimension)
          subGradient[j] <- gradient[j]
          jumpEvaluation <- estimateLikelihood(Time, proposal, theta, likelihood, model, subGradient, i, S, alphas, readIndices, UMIs, P_u, u, I, Y, TCs, Ys, Is, p_o, p_n, n, hyperParameters, priorProbabilities, priorRandomVariables, j)
          likelihood <- jumpEvaluation[['likelihood']]
          theta[j] <- jumpEvaluation[["theta"]][j]
          gradient[j] <- jumpEvaluation[["gradient"]][j]
          A[j] <- jumpEvaluation[["A"]]
          if ('f_n_failure' %in% names(jumpEvaluation)) {
            f_n_failures[i] <- TRUE
          }
        }
        posterior[i, ] <- theta
        likelihoods[i] <- likelihood
        gradients[i, ] <- gradient
        if (i > window / 2) { #if we have any rolling X step block with 5% or more mathematical failure then use model 3
          if (length(which(f_n_failures[max(c((window / 2) + 1, (i - window + 1))):i])) >= window / 20) {
            model <- 3
            break
          }
        }
      }
    }
    output <- list('posterior' = posterior, 'model' = model)
    return(output)
  }, mc.cores = cores)
  
  names(outputs) <- genes
  return(outputs)
}



