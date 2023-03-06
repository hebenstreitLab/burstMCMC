---
title: "Test"
output:
  html_document: default
  pdf_document: default
---
# Inferring transcriptional bursting dynamics from 4sU scRNA-seq data

## Installation and preabmle
``` {eval=F, echo=T}
install.packages("devtools")
devtools::install_github("hebenstreitLab/burstMCMC")
install.packages('actuar')
install.packages('parallel')
install.packages('reticulate')
reticulate::install_miniconda()
reticulate::py_install('scipy')
reticulate::py_install('numpy')
burstMCMC::getPython()
```

## Loading data
``` {eval=F, echo=T}
preprocessedData <- burstMCMC::preprocessedData
alphas <- burstMCMC::alphas

cores <- 10
Data <- formatData(preprocessedData, alphas, cores)
Data <- burstMCMC::sampleData
genes <- names(Data)
```

## Parameter inference
``` {eval=F, echo=T}
Time <- 240
model <- 2
burstMCMC::getPython()
outputs <- MALAwithinGibbs(model, cores, Time, Data)

means <- burstMCMC::means
CVs <- burstMCMC::CVs
estimates <- burstMCMC::estimates

outputs <- burstMCMC::outputs
```

## Analysing results
``` {eval=F, echo=T}
models <- getModels(outputs)


posteriors <- getPosteriors3P(outputs)

#only for the 3 parameter posterior checking
checkAcceptanceRates(posteriors)

gene <- genes[1]
plotTraces(gene, posteriors)
thinning <- 2
#applies thinning and burn-in removal
plotPosteriors(gene, posteriors, thinning)
estimates3P <- getEstimates(posteriors, thinning)

posteriors <- getPosteriors5P(outputs)
plotTraces(gene, posteriors)
plotPosteriors(gene, posteriors, thinning)
estimates5P <- getEstimates(posteriors, thinning)
```

## Simulation-based validation
``` {eval=F, echo=T}
#means estimates must be first three columns as expression level, burst rate and transcript lifetime, with genes as rownames
simulatedData <- simulateGenes(cores, Time, Data, means)
simulatedData <- simulateGenes(cores, Time, Data, estimates3P[["means"]])
outputs <- MALAwithinGibbs(model, cores, Time, simulatedData)
```

