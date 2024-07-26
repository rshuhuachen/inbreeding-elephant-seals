# Northern elephant seal inbreeding models

This repository contains all raw data, plots and scripts for the inbreeding models implemented in Hoffman et al. (2024, in press) titled "Genomic and fitness consequences of a near-extinction event in the northern elephant seal"  published in Nature Ecology and Evolution. 

To import the repository to your local device, ensure you have git installed and execute the following command in the terminal:

```
git clone https://github.com/rshuhuachen/inbreeding-elephant-seals.git
```

Please find below a description of all raw and clean data sets and the scripts used to clean up data, run models, and plot figures.

## Raw data

In the `data/raw` sub-folder, you will find a number of files:
* NES dataset 19_04_23.xlsx: this is an excel file containing the full phenotypic data set (e.g. primary classifications
  for the cause of death, individual ID's and sex, etc.) as well as microsatellite genotypes (N = 219)
* NES_notImp.raw: genotype data from RAD-seq SNPs in plink raw format (N = 74)
* elephant_seal_sequence_report.tsv: sequence report containing chromosome lengths taken from NCBI: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_021288785.2/ 

## sMLH values
In the `data/smlh` sub-folder, you will find a number of files that contain sMLH (and g2) calculations. These values are calculated using the `scripts/1_smlh_g2_hwe.R` script.

* smlh_genomewide_msats.txt: contains sMLH values based on microsatellite data 
* smlh_snp.txt: contains sMLH values based on RAD-seq SNP data
* smlh_genotypelikelihoods.txt: contains sMLH values based on GL (excluded in manuscript)
* g2snp_10000.txt: contains g2 values based on RAD-seq SNP data
* g2msats_10000.txt contains g2 values based on microsatellite data

## Clean data
In the `data/clean` sub-folder, you will find a number of cleaned up files that have been generated within the R scripts from the raw data described above.
* msat_raw.stru: this is a STRUCTURE formatted file of the microsatellite genotypes, without any filtering steps implemented. This file is generated in the `scripts/1_smlh_g2_hwe.R` script.
* msat_clean.stru: this is a STRUCTURE formatted file of the microsatellite genotypes, filtered for HWE and individuals with a lot of missing data/genotypes. this file is generated in the `scripts/1_smlh_g2_hwe.R` script.
* phenotypes in .csv and .RData format: clean file of the phenotypes, with variables renamed and restructured. This file is generated in the `scripts/2_clean_phenotypes.R` script.
* phenotypes_smlh in .csv and .RData format: similar to the phenotypes file, but with the genome-wide sMLH estimates included. This file is generated in the `scripts/2_clean_phenotypes.R` script. 


## Scripts

The scripts for the full workflow are arranged in order of the methods described in the paper. All scripts are R scripts that have been run in R version 3.6.3. Package versions are described in the paper. The output of the models is not stored on github due to the large file sizes, but can be regenerated on your local device. 

* 1_smlh_g2_hwe.R: this script takes the raw microsatellite file and tests for Hardy_Weinberg equilibrium, filters accordingly, and subsequently calculates sMLH and g2 based on the clean genotypes. Next, it calculates sMLH based on RAD-seq SNP genotypes both genome-wide and on a chromosome-by-chromosome basis.
* 2_clean_phenotypes.R: this script cleans up the raw phenotypic data file and produces the clean data file, including genome-wide sMLH estimates based on the various methods used. Adjustments include renaming of variables and creating several new ones.
* 3_model_inbreeding_cat.R: this script contains the glmers implemented in the Bayesian package brms to predict the effect of sMLH on the 6 death cause categories. Diagnostic plots of the posterior distribution can be found in plots/diagnostics/
* 4_model_inbreeding_mass_blubber.R: this script contains the glmers implemented in the Bayesian package brms to predict the effect of sMLH on body mass and blubber thickness at admittance. Diagnostic plots of the posterior distribution can be found in plots/diagnostics/
* 5_model_inbreeding_perchr.R: this script calculates sMLH effects on worm and bacterial infection on a chromosome-by-chromosome level. Diagnostic plots of the posterior distribution can be found in plots/diagnostics/per_chromosome
* 6_poweranalysis.R: this script conducts a power analysis to show what sample sizes are needed to have high power in detecting inbreeding depression in cause of death categories
* 7_manuscript_plots.R: this script contains the code used to generate the plots (main text and supplements) shown in the manuscript