# Inferring transcriptional bursting dynamics from 4sU scRNA-seq data

This is the R package associated with the paper entitled “Synergising single-cell resolution and 4sU labelling boosts inference of transcriptional bursting” to carry out Bayesian inference of transcriptional bursting dynamics parameters from 4sU scRNA-seq data using a Markov chain Monte Carlo (MCMC) algorithm. This does not include pre-processing steps for extracting the required information from the raw data since this is not a generalisable procedure and varies based on the sequencing method etc. That being said, the scripts we used for pre-processing of the 4sU scRNA-seq data obtained from the Qiu et al 2020 paper entitled "Massively parallel and time-resolved RNA sequencing in single cells with scNT-seq", which we applied our inference method to in our paper, can be found in our GitHub repository (https://github.com/hebenstreitLab/burstMCMCpreprocessing). Therefore, this package assumes we already have the information on UMI counts, T>C conversions per read, total Ts per read (in the corresponding fasta sequence) for each cell and gene, the cell-specific capture efficiencies, the gene-specific background T>C rates, and the gene-invariant 4sU-mediated T>C rate, which were obtained with our pre-processing scripts.

## Installation and preabmle

Download and installation of our package 'burstMCMC' is done via GitHub using the 'devtools' package. This will automatically install other required R packages. Then our function 'getPython' must run to install a miniconda python environment if one is not already present and install any missing python packages which are required. After this brief set-up we are in a position to use the package. 'getPython' must be run when first opening R each time, even after installations are complete in order to load the python functions which are used within the inference algorithm.
``` {R eval=F, echo=T}
install.packages("devtools")
devtools::install_github("hebenstreitLab/burstMCMC")
library(burstMCMC)
getPython()
```

## Loading data

After installation, the data must be formatted correctly to be input into the algorithm. A sample of the pre-processed data output by our pre-processing pipeline when applied to the Qiu dataset is provided for ten genes which we were able to derive high confidence parameter estimates for in the paper. We also provide the cell-specific capture efficiencies (alphas) estimated with the pre-processing scripts.
``` {R eval=F, echo=T}
preprocessedData <- burstMCMC::preprocessedData
alphas <- burstMCMC::alphas
```
These are loaded automatically with the package. The 'formatData' function can then be used to convert the pre-processed Qiu data and the alphas into the format required by the algorithm, with the option to parallelise over multiple cores. 
``` {R eval=F, echo=T}
cores <- 10
Data <- formatData(preprocessedData, alphas, cores)
```
The correctly formatted data for the ten sampled genes may also be directly loaded in through the package and its structure viewed, in which the data is partitioned by genes at the highest level. The 'genomicTs' part of the data has been pooled across cells for each gene and used to obtain an empirical probability mass function for the total number of uracils present in each read which may or may not have undergone conversion.
``` {R eval=F, echo=T}
Data <- burstMCMC::sampleData
genes <- names(Data)
View(Data)
```

## Parameter inference

``` {R eval=T, echo=F}
library(burstMCMC)
Data <- burstMCMC::sampleData
genes <- names(Data)
```

The correctly formatted data can be fed into the inference algorithm using the 'MALAwithinGibbs' function. When using model 1, only the UMIcounts and alphas need to be provided for each gene, meaning the rest can be dropped which allows for the analysis of conventional scRNA-seq data without 4sU. Otherwise, model 2 provides the the optimal results as described in the paper, although model 3 will be capable of inferring all parameters like model 2 with generally increased speed but without maximising confidence in the parameter estimates. Again, this can be parallelised over multiple cores, running one Markov chain per gene, and the 4sU pulse duration in minutes must also be provided, which in our case is 240 minutes. Output data generated using the below code, which contains the inferred posteriors and the model used in the inference for each of the ten sampled genes, can also be loaded in directly from the package for downstream analysis.
``` {R eval=F, echo=T}
Time <- 240
model <- 2
getPython()
outputs <- MALAwithinGibbs(model, cores, Time, Data)
outputs <- burstMCMC::outputs
```

We also provide the mean value and corresponding coefficient of variation (CV) for three parameters (expression level, burst rate and transcript lifetime) dervied from the posteriors for all 12276 genes we carried out inference on in the paper. In addition, we provide the mean value estimates for five parameters (expression level, burst rate and decay rate, burst size, burst frequency) just for the 584 high confidence genes from the paper.
``` {R eval=F, echo=T}
means <- burstMCMC::means
CVs <- burstMCMC::CVs
estimates <- burstMCMC::estimates
```

## Analysing results

After generating output data, there are several functions for diagnostics and analysis. Firstly, we may extract the information on which model was used in the inference with 'getModels'. This only matters when using model 2 as the algorithm may need to switch to model 3 for some genes, as described in the paper. We may also extract the Markov chain samples of the posterior distributions for each gene. This is done initially using 'getPosteriors3P', which will extract the posteriors just for the three default parameters which comprise the parameter space through which the Markov chain moves (expression level, burst rate and transcript lifetime).
``` {R eval=T, echo=T}
models <- getModels(outputs)
posteriors <- getPosteriors3P(outputs)
```

We may then generate diagnostic plots showing the distribution acceptance rates of the Markov chains for each parameter across the analysed genes, with the optimal acceptance rate indicated (0.574). They should generally be centred around this value.
``` {R eval=T, echo=T}
checkAcceptanceRates(posteriors)
```

There are also functions for examining the results for individual genes of interest to plot the Markov chain traces with 'plotTraces' and the posterior distributions we obtained using 'plotPosteriors'. Additionally we obtain mean value parameter estimates and associated CVs for all genes analysed with the 'getEstimates' function. 'plotPosteriors' and 'getEstimates' will automatically remove the pre-convergence burn-in of the Markov chain as well as apply specified thinning, as described in the methods section of the paper.
``` {R eval=T, echo=T}
gene <- genes[1]
plotTraces(gene, posteriors)
thinning <- 2
estimates3P <- getEstimates(posteriors, thinning)
plotPosteriors(gene, posteriors, thinning)
```

The same analysis may also be carried out but using posteriors calculated for the five parameters generally analysed in the paper (expression level, burst rate and decay rate, burst size, burst frequency) by executing 'getPosteriors5P' instead of 'getPosteriors3P'.
``` {R eval=T, echo=T}
posteriors <- getPosteriors5P(outputs)
plotTraces(gene, posteriors)
estimates5P <- getEstimates(posteriors, thinning)
plotPosteriors(gene, posteriors, thinning)
```

## Simulation-based validation
We may feed the mean value estimates for our three default parameters into the function 'simulateGenes', which will generate simulated data for each gene from our real dataset (Data), with the estimates as 'ground truth' values used to simulate the data. The simulated data will be in the correct format for feeding into the algorithm and can be used to test the capacity of the algorithm to recover known parameter values, as was done in the paper.
``` {R eval=F, echo=T}
simulatedData <- simulateGenes(cores, Time, Data, means)
simulatedData <- simulateGenes(cores, Time, Data, estimates3P[["means"]])
outputs <- MALAwithinGibbs(model, cores, Time, simulatedData)
```

