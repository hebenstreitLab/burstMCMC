#' @export
formatData <- function(preprocessedData, alphas, cores) {
  alphas <- sapply(names(preprocessedData[["UMIcounts"]][[names(preprocessedData[["UMIcounts"]])[1]]]), function(cell, alphas) {
    return(alphas[[cell]])
  }, alphas)
  p_n <- preprocessedData[["lambdaN"]]
  preprocessedData <- preprocessedData[1:4]
  genes <- names(preprocessedData[["UMIcounts"]])
  Data <- parallel::mclapply(genes, function(gene) {
    geneData <- sapply(preprocessedData, function(dataset, gene) {
      return(dataset[[gene]])
    }, gene, simplify = F)
    UMIcounts <- geneData[["UMIcounts"]]
    conversions <- geneData[["conversions"]]
    genomicTs <- geneData[["genomicTs"]]
    background <- geneData[['backgrounds']]
    conversions <- sapply(names(UMIcounts), function(cell, conversions) {
      return(conversions[[cell]])
    }, conversions)
    genomicTs <- sapply(names(UMIcounts), function(cell, genomicTs) {
      return(genomicTs[[cell]])
    }, genomicTs)
    UsPerRead <- table(unlist(genomicTs)) / length(unlist(genomicTs))
    geneData <- list('UMIcounts' = unlist(UMIcounts), 'conversions' = conversions, 'UsPerRead' = UsPerRead, 'alphas' = alphas, 'background' = background, 'lambdaN' = p_n)
    return(geneData)
  }, mc.cores = cores)
  names(Data) <- genes
  return(Data)
}




