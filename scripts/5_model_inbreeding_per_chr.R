##### Run a model for each chromosome sMLH (N = 18) ####
## No Bayesian framework
## Only look at case/control in a binomial glmer

### load packages ####

pacman::p_load(tidyverse, lme4, lmerTest, brms, data.table, effectsize)

### set theme
theme_set(theme_classic())

### load data
load("data/clean/phenotypes.RData")
smlh_snp <- fread("data/raw/74_Samples_sMLH.txt")
smlh_gl <- fread("data/raw/GL_Heterozygosity_dataset.txt")
load(file="data/clean/phenotype_smlh.RData")

#### Merge data ####
smlh_snp <- left_join(smlh_snp, class[,c("id", "Year", "Sex", "admit_month")], by = c("ID" = "id"))

## rename
smlh_snp <- smlh_snp %>% select(-c("Primary_classification"))
names(smlh_snp) <- c("id", "class", "smlh", "smlh_chr1", "smlh_chr2", "smlh_chr3", "smlh_chr4", "smlh_chr5", "smlh_chr6", "smlh_chr7", "smlh_chr8",
                     "smlh_chr9", "smlh_chr10", "smlh_chr11", "smlh_chr12","smlh_chr13", "smlh_chr14", "smlh_chr15", "smlh_chr16", "smlh_chr17", "year", "sex", "admit_month")

smlh_snp <- left_join(smlh_snp, data[,c("id", "primary_class")], by = "id")

#### Function to calculate effect size per chromosome for binary class worms ####

model_per_chr_worms_bayes <- function(data){
  ### Run model per chromosome
  for (i in 1:17){
    formula_chr <- formula(paste0("worms_bi ~ scale(smlh_chr", i, ") + (1|sex) + (1|year) + (1|admit_month)"))
    assign(x = paste0("model_chr_", i),
           value = brm(formula_chr, family = "bernoulli", data = data,
                       iter = 100000, 
                       thin = 100,
                       chains = 3,
                       set_prior("normal(0,1)", class = "b")))
  }
  
  ## Collect model output in list
  model_list_chr <- list()
  for (i in 1:17){
    model_list_chr[[i]] <- get(paste0("model_chr_", i))
  }
  
  ## Collect model summaries in list
  sum_model_list_chr <- list()
  for (i in 1:17){
    sum_model_list_chr[[i]] <- summary(get(paste0("model_chr_", i)))
  }
  
  ## Collect effect sizes in list
  eff_model_list_chr <- list()
  for (i in 1:17){
    eff_model_list_chr[[i]] <- effectsize(get(paste0("model_chr_", i)))
  }
  eff_model_chr <- plyr::rbind.fill(eff_model_list_chr)
  eff_model_chr <- subset(eff_model_chr, grepl("scale", Parameter))
  
  # add significance
  eff_model_chr <- eff_model_chr %>% mutate(Significance = case_when(
    CI_low < 0 & CI_high < 0 ~ "Doesn't cross 0",
    CI_low > 0 & CI_high > 0 ~ "Doesn't cross 0",
    TRUE ~ "Crosses 0"
  ))
  eff_model_chr$Parameter <- as.factor(1:17)
  
  output <- list(model_list_chr, sum_model_list_chr, eff_model_chr)
  return(output)
}
model_per_chr_bacteria_bayes <- function(data){
  ### Run model per chromosome
  for (i in 1:17){
    formula_chr <- formula(paste0("bacteria_bi ~ scale(smlh_chr", i, ") + (1|sex) + (1|year) + (1|admit_month)"))
    assign(x = paste0("model_chr_", i),
           value = brm(formula_chr, family = "bernoulli", data = data,
                       iter = 100000, 
                       thin = 100,
                       chains = 3,
                       set_prior("normal(0,1)", class = "b")))
  }
  
  ## Collect model output in list
  model_list_chr <- list()
  for (i in 1:17){
    model_list_chr[[i]] <- get(paste0("model_chr_", i))
  }
  
  ## Collect model summaries in list
  sum_model_list_chr <- list()
  for (i in 1:17){
    sum_model_list_chr[[i]] <- summary(get(paste0("model_chr_", i)))
  }
  
  ## Collect effect sizes in list
  eff_model_list_chr <- list()
  for (i in 1:17){
    eff_model_list_chr[[i]] <- effectsize(get(paste0("model_chr_", i)))
  }
  eff_model_chr <- plyr::rbind.fill(eff_model_list_chr)
  eff_model_chr <- subset(eff_model_chr, grepl("scale", Parameter))
  
  # add significance
  eff_model_chr <- eff_model_chr %>% mutate(Significance = case_when(
    CI_low < 0 & CI_high < 0 ~ "Doesn't cross 0",
    CI_low > 0 & CI_high > 0 ~ "Doesn't cross 0",
    TRUE ~ "Crosses 0"
  ))
  eff_model_chr$Parameter <- as.factor(1:17)
  
  output <- list(model_list_chr, sum_model_list_chr, eff_model_chr)
  return(output)
}

### Run function ####
models_snp_worms <- model_per_chr_worms_bayes(data = smlh_snp)
models_snp_bacteria <- model_per_chr_bacteria_bayes(data = smlh_snp)

save(models_snp_worms, file = "output/models_snp_worms_chromosome_bayes.RData")
save(models_snp_bacteria, file = "output/models_snp_bacteria_chromosome_bayes.RData")

#### Model diagnostics ####
list_models <- list(models_snp_worms_chr1 = models_snp_worms[[1]][[1]], 
                    brm_weight_snp = brm_weight_snp, 
                    brm_blubber_msat = brm_blubber_msat, 
                    brm_blubber_snp = brm_blubber_snp)

for (i in 1:length(list_models)){
  pdf(file = paste0("plots/diagnostics/per_chromosome/diagnose_", names(list_models[i]), ".pdf"))
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
