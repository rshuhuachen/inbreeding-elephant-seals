#### In this script we create a function to fit a Bayesian model to each 
#### death cause category using the different sMLH appraoches 

### Load packages ####
pacman::p_load(tidyverse, brms, bayesplot)

### Load data ####
load(file="data/clean/phenotype_smlh.RData")

### Modelling ####

# Here we create a function to run the Bayesian model for all 6 categories for the 3 approaches;
# calculating sMLH based on microsatellites, RAD-seq SNPs and genotype likelihoods (GL) based on whole genomes
# (GL not included in MS)

model_binomial <- function(response, smlh){
  pacman::p_load(tidyverse, brms, effectsize, bayesplot)
  formula_class <- formula(paste0(response, " ~ scale(", smlh, " ) + (1|Sex) + (1|Year) + (1|admit_month)"))
  model <- brm(formula_class,
                   data = data,
                   family = bernoulli,
                   cores = 8,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   iter = 100000, thin = 100, chains = 3,
                   set_prior("normal(0,1)", class = "b")) 
  return(model)}

# To iterate over the different models, here we make a dataframe containing the different category-method
# combinations to run the model on

models_to_run <- data.frame(response = rep(c("malnutrition_bi", "congenital_bi", "bacteria_bi",
                                  "protozoa_bi", "worms_bi", "trauma_bi"), times = 3),
                            smlh = rep(c("smlh_msat", "smlh_snp", "smlh_gl"), each = 6))
models_to_run$response <- as.character(models_to_run$response)
models_to_run$smlh <- as.character(models_to_run$smlh)

# And then we iterate over the models_to_run to actually run each model
for (i in 1:nrow(models_to_run)){
  assign(x = paste0("model_", models_to_run$response[i], "_", models_to_run$smlh[i]),
         value = model_binomial(response = models_to_run$response[i],
                 smlh = models_to_run$smlh[i]))}

# Save the output - due to the large sizes these outputs are not stored in github
save(model_malnutrition_bi_smlh_msat, file = "output/brms_binaryprimary_malnutrition_msat.RData")
save(model_malnutrition_bi_smlh_snp, file = "output/brms_binaryprimary_malnutrition_snp.RData")
save(model_malnutrition_bi_smlh_gl, file = "output/brms_binaryprimary_malnutrition_gl.RData")
save(model_congenital_bi_smlh_msat, file = "output/brms_binaryprimary_congenital_msat.RData")
save(model_congenital_bi_smlh_snp, file = "output/brms_binaryprimary_congenital_snp.RData")
save(model_congenital_bi_smlh_gl, file = "output/brms_binaryprimary_congenital_gl.RData")
save(model_bacteria_bi_smlh_msat, file = "output/brms_binaryprimary_bacteria_msat.RData")
save(model_bacteria_bi_smlh_snp, file = "output/brms_binaryprimary_bacteria_snp.RData")
save(model_bacteria_bi_smlh_gl, file = "output/brms_binaryprimary_bacteria_gl.RData")
save(model_protozoa_bi_smlh_msat, file = "output/brms_binaryprimary_protozoa_msat.RData")
save(model_protozoa_bi_smlh_snp, file = "output/brms_binaryprimary_protozoa_snp.RData")
save(model_protozoa_bi_smlh_gl, file = "output/brms_binaryprimary_protozoa_gl.RData")
save(model_trauma_bi_smlh_msat, file = "output/brms_binaryprimary_trauma_msat.RData")
save(model_trauma_bi_smlh_snp, file = "output/brms_binaryprimary_trauma_snp.RData")
save(model_trauma_bi_smlh_gl, file = "output/brms_binaryprimary_trauma_gl.RData")
save(model_worms_bi_smlh_msat, file = "output/brms_binaryprimary_worms_msat.RData")
save(model_worms_bi_smlh_snp, file = "output/brms_binaryprimary_worms_snp.RData")
save(model_worms_bi_smlh_gl, file = "output/brms_binaryprimary_worms_gl.RData")

#### Model diagnostics - only msat and snp ####
list_models <- list(model_malnutrition_bi_smlh_msat = model_malnutrition_bi_smlh_msat, 
                    model_malnutrition_bi_smlh_snp = model_malnutrition_bi_smlh_snp, 
                    model_congenital_bi_smlh_msat = model_congenital_bi_smlh_msat, 
                    model_congenital_bi_smlh_snp = model_congenital_bi_smlh_snp,
                    model_bacteria_bi_smlh_msat = model_bacteria_bi_smlh_msat, 
                    model_bacteria_bi_smlh_snp = model_bacteria_bi_smlh_snp,
                    model_protozoa_bi_smlh_msat = model_protozoa_bi_smlh_msat, 
                    model_protozoa_bi_smlh_snp = model_protozoa_bi_smlh_snp,
                    model_trauma_bi_smlh_msat = model_trauma_bi_smlh_msat, 
                    model_trauma_bi_smlh_snp = model_trauma_bi_smlh_snp,
                    model_worms_bi_smlh_msat = model_worms_bi_smlh_msat, 
                    model_worms_bi_smlh_snp = model_worms_bi_smlh_snp)

for (i in 1:length(list_models)){
  pdf(file = paste0("plots/diagnostics/diagnose_", names(list_models[i]), ".pdf"))
  # extract beta's for some plots
  betas <- subset(variables(list_models[[i]]), grepl("b_", variables(list_models[[i]])) & !grepl("_Intercept", variables(list_models[[i]])))
  #autocor
  print(mcmc_acf(list_models[[i]], lags = 10))
  #trace betas
  print(mcmc_trace(list_models[[i]], pars = c(betas)))
  #trace all
  print(mcmc_trace(list_models[[i]]))
  
  #rhat
  print(mcmc_rhat(brms::rhat(list_models[[i]])))
  #neff
  print(mcmc_neff(neff_ratio(list_models[[i]])))
  
  # areas
  print(mcmc_areas(list_models[[i]], pars = c(betas)))
  
  dev.off()
  i = i+1
}


#### Supplementary analysis: categorical model #####
#### Categorical: primary classification ####

prior <- c(set_prior("normal(-1,1)", dpar = "muBacterialinfection"),
           set_prior("normal(0,1)", dpar = "muCongenitaldefect"),
           set_prior("normal(0,1)", dpar = "muMalnutrition"),
           set_prior("normal(-1,1)", dpar = "muOtostrongylis"),
           set_prior("normal(-1,1)", dpar = "muProtozoa"))

brm_primary_msat <- brm(relevel(primary_class, ref = "Trauma") ~ scale(smlh_msat) + (1|Sex) + (1|Year) + (1|admit_month),
                        family= categorical,data = data, 
                        iter = 100000, 
                        thin = 100,
                        chains = 3,
                        prior=prior)

save(brm_primary_msat, file = "output/brms_primaryclass_smlh_msat.RData")

brm_primary_snp <- brm(relevel(primary_class, ref = "Trauma") ~ scale(smlh_snp) + (1|Sex) + (1|Year) + (1|admit_month),
                       family= categorical,data = data, 
                       iter = 100000, 
                       thin = 100,
                       chains = 3,
                       prior=prior)
save(brm_primary_snp, file = "output/brms_primaryclass_smlh_snp.RData")

for (i in c(brm_primary_msat, brm_primary_snp)){
  pdf(file = paste0("plots/diagnostics/diagnose_", names(list_models[i]), ".pdf"))
  # extract beta's for some plots
  betas <- subset(variables(list_models[[i]]), grepl("b_", variables(list_models[[i]])) & !grepl("_Intercept", variables(list_models[[i]])))
  #autocor
  print(mcmc_acf(list_models[[i]], lags = 10))
  #trace betas
  print(mcmc_trace(list_models[[i]], pars = c(betas)))
  #trace all
  print(mcmc_trace(list_models[[i]]))
  
  #rhat
  print(mcmc_rhat(brms::rhat(list_models[[i]])))
  #neff
  print(mcmc_neff(neff_ratio(list_models[[i]])))
  
  # areas
  print(mcmc_areas(list_models[[i]], pars = c(betas)))
  
  dev.off()
  i = i+1
}

