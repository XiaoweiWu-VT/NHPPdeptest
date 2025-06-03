A statistical method to test for dependence between nonhomogeneous Poisson processes, which can be applied to detect interactions between transcription factors based on ChIP-seq data.

--nhppfun.R: contains all functions for (1) simulating NHPP sample paths, (2) estimating NHPP intensity, and (3) evaluating GoF of estimated NHPP.

--drawfig1&2.R: (simulation study) demonstrates simulation results for true/estimated NHPP intensity functions and CDFs of interarrival times of true/estimated NHPPs.

--drawfig3.R: (real data analysis) calculates a summary matrix for the average pairwise dependence among TFs across chromosome 1 using ChIP-seq data in mouse ESCs.

--maketab1.R: (simulation study) calculates MISE/PNP for evaluating intensity estimation and GoF test based on single sample path.

--maketab2.R: (simulation study) calculates FPR for testing dependence between NHPPs.

--maketab3.R: (simulation study) calculates TPR for testing dependence between NHPPs.
