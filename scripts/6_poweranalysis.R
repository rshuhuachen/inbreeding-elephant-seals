# load libraries

pacman::p_load(pwr, pwrss)

## set colors
clr_pheno <- RColorBrewer::brewer.pal(8, "Set1")[c(1,3:5,7:8)] %>%
  color() %>% 
  set_names(nm = c("Worms", "Bacteria", "Malnutrition","Trauma", "Protozoa", "Congenital defect"))

#change brown and pink color
clr_pheno["Protozoa"] <- "#0C0A3E"
clr_pheno["Congenital defect"] <- "#8AC6D0"

# load phenotypes
load("data/clean/phenotype_smlh.RData")
data <- data %>% 
  mutate(malnutrition_bi = case_when(primary_class == "Malnutrition" ~ 1,
                                     TRUE ~ 0),
         bacteria_bi = case_when(primary_class == "Bacterial infection" ~ 1,
                                 TRUE ~ 0),
         congenital_bi = case_when(primary_class == "Congenital defect" ~ 1,
                                   TRUE ~ 0),
         worms_bi = case_when(primary_class == "Otostrongylis" ~ 1,
                              TRUE ~ 0),
         protozoa_bi = case_when(primary_class == "Protozoa" ~ 1,
                                 TRUE ~ 0),
         trauma_bi = case_when(primary_class == "Trauma" ~ 1,
                               TRUE ~ 0))
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
## plotting

ggplot(pwrtest_all, aes(x = power, y = n, col = trait, group = trait))+
  geom_point(size=3)+ geom_line(lwd=1)+
  scale_y_log10(breaks = c(100, 1000, 10000, 100000),
                labels = c(expression("10"^2), expression("10"^3), 
                           expression("10"^4), expression("10"^5)))+
  scale_color_manual(values=clr_pheno)+
  geom_hline(aes(yintercept = min(pwrtest_all$n[which(pwrtest_all$power == 0.8)])),
             col = "red", linetype = "dashed")+
  labs(x = "Statistical power", y = "Sample size required", col = "Category")+
  theme_bw(base_family = "Arial")+
  theme(text = element_text(size = 24, color = "black"),
        axis.text.y = element_text(color = "black", hjust=0.95),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.3),
        axis.ticks = element_line(colour = "black",
                                  linewidth = 0.3),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    color = "black"),
        plot.subtitle = element_text(hjust = .5),
        strip.background = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_text(color = "black"),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),
                                    color = "black")) -> plot_power

ggsave(plot_power, file = "plots/final/Plot_power.png", width = 12, height=10)
