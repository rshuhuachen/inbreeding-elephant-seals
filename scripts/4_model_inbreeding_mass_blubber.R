#### Here we'll model the effects of sMLH on traits that are not disease/causes of death,
#### i.e. body mass and blubber thickness

### load packages ####

pacman::p_load(tidyverse, brms, bayesplot)

### load data
load(file="data/clean/phenotype_smlh.RData")

### Modelling ####

##### Admit weight ####
data$admit_weight <- data$`Admit weight`
brm_weight_msat <- brm(admit_weight ~ scale(smlh_msat) + (1|admit_month) + (1|Year) + (1|Sex),
                    family= gaussian, data = data, 
                    iter = 100000, 
                    thin = 100,
                    chains = 3,
                    set_prior("normal(0,1)", class = "b"))

brm_weight_snp <- brm(admit_weight ~ scale(smlh_snp) + (1|admit_month) + (1|Year) + (1|Sex),
                       family= gaussian, data = data, 
                       iter = 100000, 
                       thin = 100,
                       chains = 3,
                       set_prior("normal(0,1)", class = "b"))

##### Blubber ####
data$blubber <- data$`Blubber thickness`
brm_blubber_msat <- brm(blubber ~ scale(smlh_msat) + (1|admit_month) + (1|Year) + (1|admit_month),
                       family= gaussian, data = data, 
                       iter = 100000, 
                       thin = 100,
                       chains = 3,
                       set_prior("normal(0,1)", class = "b"))

brm_blubber_snp <- brm(blubber ~ scale(smlh_snp) + (1|admit_month) + (1|Year) + (1|Sex),
                        family= gaussian, data = data, 
                        iter = 100000, 
                        thin = 100,
                        chains = 3,
                        set_prior("normal(0,1)", class = "b"))



#save models
save(brm_weight_msat, file = "output/brms_weight_smlh_msat.RData")
save(brm_weight_snp, file = "output/brms_weight_smlh_snp.RData")
save(brm_blubber_msat, file = "output/brms_blubber_smlh_msat.RData")
save(brm_blubber_snp, file = "output/brms_blubber_smlh_snp.RData")

#### Model diagnostics ####
list_models <- list(brm_weight_msat = brm_weight_msat, 
                    brm_weight_snp = brm_weight_snp, 
                    brm_blubber_msat = brm_blubber_msat, 
                    brm_blubber_snp = brm_blubber_snp)

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
