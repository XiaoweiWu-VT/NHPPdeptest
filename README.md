A statistical method to test for dependence between nonhomogeneous Poisson processes, which can be applied to detect interactions between transcription factors based on ChIP-seq data. This website contains the complete source code (including both simulation and real data analysis) and datasets (including all simulated and demo data), as listed in the following table:

Content                                 |  Code	            Data
                                        |  nhppfun.R
Demo (for Figure 1 and Figure 2)        |  drawfig1&2.R      demo_arrtime.RData
                                                            demo_intensity.RData
                                                            demo_intensityestimate.RData
                                                            demo_intarrtime.RData
                                                            demo_cdfintarrtime.RData
Simulation 1 (for Table 1 and Table 2)  |  maketab1&2.R      simu1_arrtime_scenI_npath1.RData
                                                            simu1_arrtime_scenII_npath1.RData
                                                            simu1_arrtime_scenIII_npath1.RData
                                                            simu1_arrtime_scenI_npath5.RData
                                                            simu1_arrtime_scenII_npath5.RData
                                                            simu1_arrtime_scenIII_npath5.RData
Simulation 2 under H0 (for Table 3)       simu2h0sever.R    simu2h0_scenIV.RData
                                          maketab3.R        simu2h0_scenV.RData
                                                            simu2h0_scenVI.RData
Simulation 2 under Ha (for Table 4)       simu2hasever.R    simu2ha_arrtimes_scenI_rho_-1.RData
                                          maketab4.R        simu2ha_arrtimes_scenI_rho_-0.5.RData
                                                            simu2ha_arrtimes_scenI_rho_-0.3.RData
                                                            simu2ha_arrtimes_scenI_rho_0.RData
                                                            simu2ha_arrtimes_scenI_rho_0.1.RData
                                                            simu2ha_arrtimes_scenI_rho_1.RData
                                                            simu2ha_arrtimes_scenII_rho_-1.RData
                                                            simu2ha_arrtimes_scenII_rho_-0.5.RData
                                                            simu2ha_arrtimes_scenII_rho_-0.3.RData
                                                            simu2ha_arrtimes_scenII_rho_0.RData
                                                            simu2ha_arrtimes_scenII_rho_0.1.RData
                                                            simu2ha_arrtimes_scenII_rho_1.RData
                                                            simu2ha_arrtimes_scenIII_rho_-1.RData
                                                            simu2ha_arrtimes_scenIII_rho_-0.5.RData
                                                            simu2ha_arrtimes_scenIII_rho_-0.3.RData
                                                            simu2ha_arrtimes_scenIII_rho_0.RData
                                                            simu2ha_arrtimes_scenIII_rho_0.1.RData
                                                            simu2ha_arrtimes_scenIII_rho_1.RData
Real Data Application (for Figure 3)      drawfig3.R        mmc2.xls [2]
                                                            mmc3.xls [1]
                                                            mES_OCT4_SOX2_NANOG_TCF3.WIG [3]


--nhppfun.R: contains all functions for (1) simulating NHPP sample paths, (2) estimating NHPP intensity, and (3) evaluating GoF of estimated NHPP.

--drawfig1&2.R: (simulation study) demonstrates simulation results for true/estimated NHPP intensity functions and CDFs of interarrival times of true/estimated NHPPs.

--drawfig3.R: (real data analysis) calculates a summary matrix for the average pairwise dependence among TFs across chromosome 1 using ChIP-seq data in mouse ESCs.

--maketab1.R: (simulation study) calculates MISE/PNP for evaluating intensity estimation and GoF test based on single sample path.

--maketab2.R: (simulation study) calculates FPR for testing dependence between NHPPs.

--maketab3.R: (simulation study) calculates TPR for testing dependence between NHPPs.
