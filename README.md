# "Fruit odorants mediate co-specialization in a multispecies plant-animal mutualism"

Published: https://royalsocietypublishing.org/doi/10.1098/rspb.2021.0312

Required R packages:
vegan
fps
cluster
parallelsugar
geiger
phytools
mvMORPH
l1ou

0. Contact ssantana@uw.edu for 2-3, liliana.davalos@stonybrook.edu for 4-7.

1. Place all scripts in the same folder as the data.

2. VOC_multivariate_analyses.R runs the Multivariate Dimensional Scaling (MDS) analyses, plots and saves the results, and identifies the VOCs that significantly explain variation along MDS1-5.

3. VOC_phylogenetic_analyses.R time-calibrates de ITS phylogeny, plots traits on the phylogeny and conducts cophylogenetic comparisons, fits multivariate (BM, OU, EB) models on variables, fits transition rates models, and runs l1ou models and their simulations.

4. setup_VOCs.R organizes the data for multivariate models in brm/Stan. Saves setup_VOCs.RData for downstream analyses.

5. BayesVOC_v4.R runs multivariate models with +1 for log transformations. Several analyses must run more than once to make sure they are converging. Modify this script for 0.1, 1, and 10 constants and change the save file. Saves multiresponse_v3.RData for plotting.

6. plotting_models.R runs and plots leave one out posterior predictive checks and model coefficients. 

7. OPTIONAL picking_models.R if you ran multiple constants on step 5, use this script to compare leave one out (loo) statistics.
