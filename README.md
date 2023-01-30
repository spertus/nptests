# nptests

R functions to conduct nonparametric inference. Install from the R console by running `devtools::install_github('spertus/nptests')`. 

Contains two functions:

1) `gaffke_CI()` constructs a nonparametric confidence bound on the mean of a bounded distribution using IID samples from that distribution.
2) `two_sample_gaffke()` uses IID samples from each of two bounded distributions and tests whether the mean of one distribution is larger than that of the other. 

The inference is finite-sample, nonparametric (FSNP) in the sense that, if the bounds truly contain the support of the population distribution, the resulting P-values and confidence intervals are valid: P-values are dominated by the uniform distribution and confidence intervals have coverage greater than their nominal level. 

This validity property is _conjectured_, not proven, but there is considerable evidence that it is true. For more details see:

1) Gaffke, N. 'Three test statistics for a nonparametric one-sided hypothesis on the mean of a nonnegative variable.' (2004) https://www.math.uni-magdeburg.de/institute/imst/ag_gaffke/files/pp1304.pdf
2) Learned-Miller, E and Thomas, P. 'A new confidence interval for the mean of a bounded random variable.' (2019) http://arxiv.org/abs/1905.06208
