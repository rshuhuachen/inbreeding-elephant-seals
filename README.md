# Northern elephant seal inbreeding models

This repository contains all raw data, plots and scripts for the inbreeding models implemented in Hoffman et al. (in prep) Genomic consequences of an extreme bottleneck in northern elephant seals.

To import the repository to your local device, ensure you have git installed and execute the following command in the terminal:

```
git clone https://github.com/rshuhuachen/inbreeding-elephant-seals.git
```

Please find below a description of all raw data sets and the scripts

## Raw data

In the `data/raw` sub-folder, you will find a number of files:
* NES dataset 19_04_23.xlsx: this is an excel file containing the full phenotypic data set (e.g. primary classifications
  for the cause of death, individual ID's and sex, etc.)
* 


## Scripts

The scripts for the full workflow are arranged in order of the methods described in the paper. All scripts are R scripts that have been run in R version 3.6.3. Package versions are described in the paper. The output of the models is not stored on github due to the large file sizes, but can be regenerated on your local device. 

* 1_smlh_g2_hwe.R: this script takes the raw microsatellite file and tests for Hardy_Weinberg equilibrium, filters accordingly, and subsequently calculates sMLH and g2 based on the clean genotypes
* 2_clean_phenotypes.R: this script cleans up the raw phenotypic data file and produces the clean data file, including genome-wide sMLH estimates based on the various methods used. Adjustments include renaming of variables and creating several new ones.
* 3_model_inbreeding_cat.R: this script contains the glmers implemented in the Bayesian package brms to predict the effect of sMLH on the 6 death cause categories. Diagnostic plots of the posterior distribution can be found in plots/diagnostics/
* 4_model_inbreeding_mass_blubber.R: this script contains the glmers implemented in the Bayesian package brms to predict the effect of sMLH on body mass and blubber thickness at admittance. Diagnostic plots of the posterior distribution can be found in plots/diagnostics/
* 5_model_inbreeding_perchr.R: this script calculates sMLH effects on worm and bacterial infection on a chromosome-by-chromosome level. Diagnostic plots of the posterior distribution can be found in plots/diagnostics/per_chromosome
* 6_poweranalysis.R: this script conducts a power analysis to show what sample sizes are needed to have high power in detecting inbreeding depression in cause of death
* 7_manuscript_plots.R: this script contains the code used to generate the plots (main text and supplements) shown in the manuscript