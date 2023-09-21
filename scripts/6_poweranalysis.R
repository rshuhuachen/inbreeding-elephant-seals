# load libraries

pacman::p_load(pwr, pwrss, prismatic)

## set colors
clr_pheno <- RColorBrewer::brewer.pal(8, "Set1")[c(1,3:5,7:8)] %>%
  color() %>% 
  set_names(nm = c("Worms", "Bacteria", "Malnutrition","Trauma", "Protozoa", "Congenital defect"))

#change brown and pink color
clr_pheno["Protozoa"] <- "#0C0A3E"
clr_pheno["Congenital defect"] <- "#8AC6D0"

# load phenotypes
load("data/clean/phenotype_smlh.RData")

# # power test for binary tests
# data %>% group_by(malnutrition_bi) %>% summarise(n = n(),
#                                                  mean_smlh = mean(smlh_snp, na.rm=T),
#                                                  sd_smlh = sd(smlh_snp, na.rm=T)) -> summary_smlh_mal
# kappa_mal <- summary_smlh_mal$n[1] / summary_smlh_mal$n[2] 
# 
# pwrss.t.2means(mu1 = summary_smlh_mal$mean_smlh[1],
#                mu2 = summary_smlh_mal$mean_smlh[2],
#                sd1 = summary_smlh_mal$sd_smlh[1],
#                sd2 = summary_smlh_mal$sd_smlh[2],
#               # n2 = summary_smlh_mal$n[2],
#                power = 0.80,
#                kappa = kappa_mal,
#                alpha = 0.05,
#                alternative = "not equal") -> mal

#pwrss.f.reg(f2=0.1, m = 1, k = 3, power = 0.8, alpha = 0.05)
#pwr.t.test(sig.level = 0.05, power = 0.8, type = "two.sample", alternative = "greater", d = 0.1)
#pwr.t2n.test(n1 = 100, n2 = 40, d = 0.1, alternative = "two.sided")
### function

power_test <- function(trait, pwr_increment, data){
  summary <- data %>% group_by({{trait}}) %>% 
    summarise(n = n(),
              mean_smlh = mean(smlh_snp, na.rm=T),
              sd_smlh = sd(smlh_snp, na.rm=T)) 
  
   kappa <- summary$n[1] / summary$n[2] 
   
  #power levels
  pwr <- seq(from = 0.1, to = 0.9, by = pwr_increment)

  pwrtest <- data.frame()
  for (i in 1:length(pwr)){
    pwrtest_single <<- pwrss.t.2means(mu1 = summary$mean_smlh[1],
                              mu2 = summary$mean_smlh[2],
                              sd1 = summary$sd_smlh[1],
                              sd2 = summary$sd_smlh[2],
                              power = pwr[i],
                              kappa = kappa,
                              alpha = 0.05,
                              alternative = "not equal",
                              verbose=F)
    pwrtest_df <- as.data.frame(t(as.data.frame(unlist(pwrtest_single))))
    pwrtest <- rbind(pwrtest, pwrtest_df)
    }
   return(pwrtest)
  }

## run power test for all 6 categories - for some reason can't loop
pwr_mal <- power_test(malnutrition_bi, pwr_increment = 0.1, data = data)
pwr_mal$trait <- "Malnutrition"
pwr_bac <- power_test(bacteria_bi, pwr_increment = 0.1, data = data)
pwr_bac$trait <- "Bacteria"
pwr_cong <- power_test(congenital_bi, pwr_increment = 0.1, data = data)
pwr_cong$trait <- "Congenital defect"
pwr_worm <- power_test(worms_bi, pwr_increment = 0.1, data = data)
pwr_worm$trait <- "Worms"
pwr_pro <- power_test(protozoa_bi, pwr_increment = 0.1, data = data)
pwr_pro$trait <- "Protozoa"
pwr_trau <- power_test(trauma_bi, pwr_increment = 0.1, data = data)
pwr_trau$trait <- "Trauma"

pwrtest_all <- rbind(pwr_mal, pwr_bac, pwr_cong, pwr_worm, pwr_pro, pwr_trau)
#reformat
pwrtest_all$n.n1 <- as.numeric(as.character(pwrtest_all$n.n1))
pwrtest_all$n.n2 <- as.numeric(as.character(pwrtest_all$n.n2))
pwrtest_all$n <- pwrtest_all$n.n1 + pwrtest_all$n.n2

save(pwrtest_all, file = "output/powertest_all.RData")

