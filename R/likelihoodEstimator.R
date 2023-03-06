#' @export
estimateLikelihood <- function(Time, proposal, theta, lastLikelihood, model, lastGradient, step, S, alphas, readIndices, UMIs, P_u, u, I, Y, TCs, Ys, Is, p_o, p_n, n, hyperParameters, priorProbabilities, priorRandomVariables, j = 1:dimension) {
  dimension <- 3
  lastPrior <- prod(sapply(1:length(theta), function(index, theta, priorProbabilities, hyperParameters) {
    return(priorProbabilities[[index]](theta[index], hyperParameters, index))
  }, theta, priorProbabilities, hyperParameters))
  
  proposalPrior <- prod(sapply(1:length(proposal), function(index, proposal, priorProbabilities, hyperParameters) {
    return(priorProbabilities[[index]](proposal[index], hyperParameters, index))
  }, proposal, priorProbabilities, hyperParameters))
  
  if (proposalPrior == 0) {
    return(list('theta' = theta, 'likelihood' = lastLikelihood, 'A' = 0, 'gradient' = rep(0, dimension)))
  }
  
  mu <- proposal[1]
  a <- proposal[2]
  d <- 1 / proposal[3]
  k <- a * d
  b <- mu / a
  p <- exp(- d * Time)
  
  M <- qnbinom(0.9999, size = a, mu = mu)
  f_m <- dnbinom(0:M, size = a, mu = mu)
  
  if (model == 1) {
    
    L1 <- sum(sapply(1:n, function(cellIndex, UMIs, alphas, M, f_m) {
      l <- UMIs[cellIndex]
      alpha <- alphas[cellIndex]
      return(log(sum(dpois(l, 0:M * alpha) * f_m)))
    }, UMIs, alphas, M, f_m))
    
    proposalLikelihood <- L1
    
    
  } else if (model == 2) {
    
    tau <- Time * d
    
    f_n  <- calculateHypergeom(a = a, b = b, tau = tau, M = M)
    
    if (length(f_n) == 1) {
      if (f_n == 0) {
        return(list('theta' = theta, 'likelihood' = lastLikelihood, 'A' = 0, 'gradient' = lastGradient, 'f_n_failure' = TRUE))
      } else { #when there is > 0.9999 probability of seeing 0 transcripts for proposal, reject automatically
        return(list('theta' = theta, 'likelihood' = lastLikelihood, 'A' = 0, 'gradient' = lastGradient))
      }
    } else {
      
      f_o <- dnbinom(0:M, size = a, mu = p * mu) #simplify the intergration of the negative binomial with the binomial distribution by simply multiplying the mu parameter by the transcript survival probability
      
      f_m_la <- sapply(readIndices, function(index, M, UMIs, alphas, f_m) {
        l <- UMIs[index]
        a <- alphas[index]
        fm_la <- f_m * dpois(l, 1:M * a) / sum(f_m * dpois(l, 1:M * a))
        return(fm_la)
      }, M, UMIs, alphas, f_m[1:M + 1])
      f_m_la[which(is.nan(f_m_la), arr.ind = T)] <- 0
      
      P_i <- sapply(0:I, function(i) {
        return(c(sum(dpois(i, u * (p_n + p_o)) * P_u), sum(dpois(i, u * p_o) * P_u)))
      })
      rownames(P_i) <- c('p_n', 'p_o')
      
      f_i_m <- get_f_i_m(M, I, f_n, f_o, P_i)
      
      f_i_la <- f_i_m %*% f_m_la
      
      L2 <- sum(sapply(1:length(readIndices), function(index, f_i_la, dmulti, Ys, Is) {
        y <- Ys[[index]]
        i <- Is[[index]]
        p_i <- f_i_la[i + 1, index]
        return(dmulti(y, p_i))
      }, f_i_la, dmulti, Ys, Is))
      
    }
    
    L1 <- sum(sapply(1:n, function(cellIndex, UMIs, alphas, M, f_m) {
      l <- UMIs[cellIndex]
      alpha <- alphas[cellIndex]
      return(log(sum(dpois(l, 0:M * alpha) * f_m)))
    }, UMIs, alphas, M, f_m))
    
    proposalLikelihood <- L1 + L2
    
  } else if (model == 3) {
    
    #just use simple bulk method of calculating probability of observed TC counts in this case
    ntr <- 1 - p
    Pi <- sapply(0:I, function(i) {
      return(sum(P_u * (ntr * dpois(i, u * (p_n + p_o)) + p * dpois(i, u * p_o))))
    })
    L2 <- dmulti(Y, Pi)
    L1 <- sum(sapply(1:n, function(cellIndex, UMIs, alphas, M, f_m) {
      l <- UMIs[cellIndex]
      alpha <- alphas[cellIndex]
      return(log(sum(dpois(l, 0:M * alpha) * f_m)))
    }, UMIs, alphas, M, f_m))
    proposalLikelihood <- L1 + L2
    
  }
  
  if (lastLikelihood > -Inf) {
    gradient <- (proposalLikelihood - lastLikelihood) / (proposal - theta)
  } else {
    gradient <- lastGradient
  }
  acceptanceProbability <- exp(min(c(0, (log(proposalPrior) + proposalLikelihood - (log(lastPrior) + lastLikelihood)))))
  
  if (!is.nan(acceptanceProbability)) {
    if (acceptanceProbability >= runif(1)) {
      return(list('theta' = proposal, 'likelihood' = proposalLikelihood, 'A' = acceptanceProbability, 'gradient' = gradient))
    } else {
      return(list('theta' = theta, 'likelihood' = lastLikelihood, 'A' = acceptanceProbability, 'gradient' = rep(0, dimension)))
    }
  } else {
    return(list('theta' = theta, 'likelihood' = lastLikelihood, 'A' = 0, 'gradient' = lastGradient))
  }
}

