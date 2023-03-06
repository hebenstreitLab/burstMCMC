#' @export
generateProposal <- function(theta, A = 0, S = 0, N = 0, gradient) {
  N <- N * 0.999
  S <- exp(log(S) + N * (A - 0.574))
  proposal <- theta + S * gradient + sqrt(2 * S) * rnorm(1)
  return(list('proposal' = proposal, 'S' = S, 'N' = N))
}

