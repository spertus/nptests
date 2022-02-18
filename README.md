# nptests

R functions to conduct nonparametric inference. Install from the R console by running `devtools::install_github('spertus/nptests')`. 

Contains two functions:

1) `gaffke_CI()` constructs a nonparametric confidence bound on the mean of a bounded distribution using IID samples from that distribution.
2) `two_sample_gaffke()` takes IID samples from each of two bounded distributions and tests (nonparametrically) whether the mean of one distribution is larger than that of the other. 

The inference is nonparametric in the sense that, if the bounds truly contain the support of the population distribution, the resulting P-values and confidence intervals are valid: P-values are dominated by the uniform distribution and confidence intervals have coverage greater than their nominal level.   
