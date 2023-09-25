##### Plotting for MS ######

### load packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, extrafont, cowplot, performance, data.table)

### set colours categories ####
clr_pheno <- RColorBrewer::brewer.pal(8, "Set1")[c(1,3:5,7:8)] %>%
  color() %>% 
  set_names(nm = c("Worms", "Bacteria", "Malnutrition","Trauma", "Protozoa", "Congenital defect"))

#change brown and pink color
clr_pheno["Protozoa"] <- "#0C0A3E" #dark
clr_pheno["Protozoa"] <- "#DEC102" #yellow
clr_pheno["Congenital defect"] <- "#8AC6D0" #light blue
clr_pheno["Congenital defect"] <- "#E4A1A0" #pink

clr_method <- color(c("#333333", "#ACB2A9"))%>% #two shades of grey
  set_names(nm = c("22 microsatellites", "15,051 SNPs"))

clr_method <- color(c("#7570b3", "#d95f02")) %>% #original
  set_names(nm = c("22 microsatellites", "15,051 SNPs"))

clr_method <- c("#D1CE9E", "#37515F") %>% color() %>% #blue and yellow ish
  set_names(nm = c("22 microsatellites", "15,051 SNPs"))

clr_method <- color(c("#37515F", "#7198AD")) %>% #blues
  set_names(nm = c("22 microsatellites", "15,051 SNPs"))

clr_method <- color(c("#7570b3", "#DEC102")) %>% #original purple and yellow instead
  set_names(nm = c("22 microsatellites", "15,051 SNPs"))

clr_method <- color(c("grey60", "grey80")) %>% #kosmas grey
  set_names(nm = c("22 microsatellites", "15,051 SNPs"))

##### Plot 1a g2 David #####

# Import data
g2<-read.table("data/smlh/g2snp_10000.txt",h=F) # 15,051 SNPs g2
g2_m<-read.table("data/smlh/g2msats_10000.txt",h=F) # msats g2
#combine
names(g2) <- "g2"
names(g2_m) <- "g2"
g2$method <- "15,051 SNPs"
g2_m$method <- "22 microsatellites"
g2_combined <- rbind(g2, g2_m)

## plot
ggplot(g2_combined, aes(x = g2, fill = method, color = method)) +
  geom_histogram(position = "identity", bins = 100) +
  scale_fill_manual(values = alpha(clr_method, 0.4))+
  scale_color_manual(values = clr_method)+
  xlim(-0.02, 0.04)+
  ylim(0, 700)+
  geom_segment(x = 0.01108032, #point estimate bootstrap
               xend = 0.01108032,
               y = 0, yend = 650, colour = clr_method["15,051 SNPs"], linewidth = 1) +
  geom_segment(x = 0.003214729, #95% CI from inbreedR
               xend = 0.020233918,#95% CI from inbreedR
               y = 650, yend = 650, colour = clr_method["15,051 SNPs"], linewidth = 1) +
  geom_segment(x = 0.01232521,#point estimate bootstrap
               xend = 0.01232521,
               y = 0, yend = 600, colour = clr_method["22 microsatellites"], linewidth = 1) +
  geom_segment(x = -0.001410787,#95% CI from inbreedR
               xend = 0.027361777,#95% CI from inbreedR
               y = 600, yend = 600, colour = clr_method["22 microsatellites"], linewidth = 1) +
  labs(x=substitute(italic(g)[2]), y = "Frequency")+
  theme_classic(base_family = "Arial")+
  theme(legend.position="none",
        text = element_text(family = "Arial", size = 24),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.3),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    color = "black"),
        plot.margin = margin(1,1,1,1, "cm"),

        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    color = "black")) -> plot_1a_g2


ggsave(plot_1a_g2, filename = "plots/final/Plot1a_g2.png", width=8, height=8)

##### Plot 1b: other factors #####
#load models
load(file = "output/brms_weight_smlh_msat.RData")
load(file = "output/brms_weight_smlh_snp.RData")
load(file = "output/brms_blubber_smlh_msat.RData")
load(file = "output/brms_blubber_smlh_snp.RData")

#combine
brms_other_interval <- rbind(mcmc_intervals_data(brm_weight_msat, prob =0.8, prob_outer = 0.95),
                             mcmc_intervals_data(brm_weight_snp, prob =0.8, prob_outer = 0.95),
                             mcmc_intervals_data(brm_blubber_msat, prob =0.8, prob_outer = 0.95),
                             mcmc_intervals_data(brm_blubber_snp, prob =0.8, prob_outer = 0.95))

brms_other_interval$model <- c(rep("22 microsatellites", nrow(mcmc_intervals_data(brm_weight_msat))),
                               rep("15,051 SNPs", nrow(mcmc_intervals_data(brm_weight_snp))),
                               rep("22 microsatellites", nrow(mcmc_intervals_data(brm_blubber_msat))),
                               rep("15,051 SNPs", nrow(mcmc_intervals_data(brm_blubber_snp))))

brms_other_interval$response <- c(rep("Body mass", nrow(mcmc_intervals_data(brm_weight_msat))),
                                  rep("Body mass", nrow(mcmc_intervals_data(brm_weight_snp))),
                                  rep("Blubber thickness", nrow(mcmc_intervals_data(brm_blubber_msat))),
                                  rep("Blubber thickness", nrow(mcmc_intervals_data(brm_blubber_snp))))

brms_other_interval <- subset(brms_other_interval, !grepl("b_Intercept", parameter))
brms_other_interval <- subset(brms_other_interval, grepl("b_", parameter))
brms_other_interval$model <- factor(as.factor(brms_other_interval$model), levels = c("22 microsatellites", "15,051 SNPs"))

brms_other_area <- rbind(mcmc_areas_data(brm_weight_msat),
                         mcmc_areas_data(brm_weight_snp),
                         mcmc_areas_data(brm_blubber_msat),
                         mcmc_areas_data(brm_blubber_snp))

brms_other_area$model <- c(rep("22 microsatellites", nrow(mcmc_areas_data(brm_weight_msat))),
                           rep("15,051 SNPs", nrow(mcmc_areas_data(brm_weight_snp))),
                           rep("22 microsatellites", nrow(mcmc_areas_data(brm_blubber_msat))),
                           rep("15,051 SNPs", nrow(mcmc_areas_data(brm_blubber_snp))))

brms_other_area$response <- c(rep("Body mass", nrow(mcmc_areas_data(brm_weight_msat))),
                              rep("Body mass", nrow(mcmc_areas_data(brm_weight_snp))),
                              rep("Blubber thickness", nrow(mcmc_areas_data(brm_blubber_msat))),
                              rep("Blubber thickness", nrow(mcmc_areas_data(brm_blubber_snp))))

brms_other_area <- subset(brms_other_area, !grepl("b_Intercept", parameter))
brms_other_area <- subset(brms_other_area, grepl("b_", parameter))
brms_other_area$model <- factor(as.factor(brms_other_area$model), levels = c("22 microsatellites", "15,051 SNPs"))

### plot

# split by interval
data_other <- split(brms_other_area, brms_other_area$interval)

data_other$bottom <- data_other$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = subset(data_other$outer, response == "Body mass")) + 
  aes(x = .data$x, y = .data$response) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model), size = 0.8)+
  geom_segment(data = subset(brms_other_interval, response == "Body mass"), aes(x = l, xend = h, yend = response), col = "black", linewidth=3)+
  geom_segment(data = subset(brms_other_interval, response == "Body mass"), aes(x = ll, xend = hh, yend = response), col = "black")+
  geom_point(data = subset(brms_other_interval, response == "Body mass"), aes(x = m, y = response), color="black", fill = "grey60", shape=21, size = 4) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  facet_wrap(~model, ncol = 2) +
  scale_fill_manual(values =alpha(clr_method, 0.6)) +
  scale_color_manual(values =clr_method) +
  #  xlim(-1.5, 1.5)+
  labs(x = "Standardised beta coefficient")+
  theme_bw(base_family = "Arial")+
  scale_y_discrete("Trait", labels=scales::label_wrap(5))+
  theme(axis.text = element_text(size=20, color = "black"),
        text = element_text(size = 24),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.3),
        axis.ticks = element_line(colour = "black",
                                  linewidth = 0.3),
        plot.subtitle = element_text(hjust = .5),
        strip.background = element_blank(),
        legend.position = "none",
        plot.margin = margin(1,1,0,1, "cm"),      
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) -> plot_1b_weight

ggplot(data = subset(data_other$outer, response == "Blubber thickness")) + 
  aes(x = .data$x, y = .data$response) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model), size = 0.8)+
  geom_segment(data = subset(brms_other_interval, response == "Blubber thickness"), aes(x = l, xend = h, yend = response), col = "black", linewidth=3)+
  geom_segment(data = subset(brms_other_interval, response == "Blubber thickness"), aes(x = ll, xend = hh, yend = response), col = "black")+
  geom_point(data = subset(brms_other_interval, response == "Blubber thickness"), aes(x = m, y = response), color="black", fill = "grey60", shape=21, size = 4) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  facet_wrap(~model, ncol = 2) +
  scale_fill_manual(values =alpha(clr_method, 0.6)) +
  scale_color_manual(values =clr_method) +
  #  xlim(-1.5, 1.5)+
  labs(x = "Standardised beta coefficient")+
  theme_bw(base_family = "Arial")+
  scale_y_discrete("Trait", labels=scales::label_wrap(5))+
  theme(axis.text = element_text(size=20, color = "black"),
        text = element_text(size = 24),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.3),
        axis.ticks = element_line(colour = "black",
                                  linewidth = 0.3),
        plot.subtitle = element_text(hjust = .5),
        strip.background = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,1,1,1, "cm"),      
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    color = "black")) -> plot_1b_blubber


## separate axes
title <- ggdraw() + 
  draw_label(
    "Trait",
    x = 0,
    angle = 90,
    hjust = 0,
    size = 26)+ theme_bw(base_family = "Arial")+
  theme(plot.margin = margin(0, 0, 0, 1, "cm"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()
  )
    
  
plot_grid(plot_1b_weight, plot_1b_blubber, nrow = 2, align = "hv", axis = "lb", rel_heights = c(1, 0.9)) -> plot_1b_blubber_weight
plot_grid(title, plot_1b_blubber_weight, rel_widths = c(0.1, 1)) -> plot_1b_blubber_weight # add title
ggsave(plot_1b_blubber_weight, filename = "plots/final/Plot1b_other_bayes.png", width=10, height=8)

##### Plot 1c: Primary classification presence/absence bayes #####
#load data

load( file = "output/brms_binaryprimary_malnutrition_msat.RData")
load( file = "output/brms_binaryprimary_malnutrition_snp.RData")
load( file = "output/brms_binaryprimary_congenital_msat.RData")
load( file = "output/brms_binaryprimary_congenital_snp.RData")
load( file = "output/brms_binaryprimary_bacteria_msat.RData")
load( file = "output/brms_binaryprimary_bacteria_snp.RData")
load( file = "output/brms_binaryprimary_protozoa_msat.RData")
load( file = "output/brms_binaryprimary_protozoa_snp.RData")
load( file = "output/brms_binaryprimary_trauma_msat.RData")
load( file = "output/brms_binaryprimary_trauma_snp.RData")
load( file = "output/brms_binaryprimary_worms_msat.RData")
load( file = "output/brms_binaryprimary_worms_snp.RData")

### combine results in one df
brms_primary_bi_interval <- rbind(mcmc_intervals_data(model_congenital_bi_smlh_msat, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_congenital_bi_smlh_snp, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_bacteria_bi_smlh_msat, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_bacteria_bi_smlh_snp, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_protozoa_bi_smlh_msat, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_protozoa_bi_smlh_snp, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_trauma_bi_smlh_msat, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_trauma_bi_smlh_snp, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_worms_bi_smlh_msat, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_worms_bi_smlh_snp, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_malnutrition_bi_smlh_msat, prob =0.8, prob_outer = 0.95),
                                  mcmc_intervals_data(model_malnutrition_bi_smlh_snp, prob =0.8, prob_outer = 0.95))
# mcmc_intervals_data(model_malnutrition_bi_smlh_gl, prob =0.8, prob_outer = 0.9))

brms_primary_bi_interval <- subset(brms_primary_bi_interval, grepl("smlh", parameter))

brms_primary_bi_interval$model <- rep(c("22 microsatellites", "15,051 SNPs"), times = 6)

brms_primary_bi_interval$response <- rep(c("Congenital defect", "Bacteria", "Protozoa", "Trauma", "Worms", "Malnutrition"),
                                         each = 2)

#areas
brms_primary_bi_area <- rbind(mcmc_areas_data(model_congenital_bi_smlh_msat),
                              mcmc_areas_data(model_congenital_bi_smlh_snp),
                              mcmc_areas_data(model_bacteria_bi_smlh_msat),
                              mcmc_areas_data(model_bacteria_bi_smlh_snp),
                              mcmc_areas_data(model_protozoa_bi_smlh_msat),
                              mcmc_areas_data(model_protozoa_bi_smlh_snp),
                              mcmc_areas_data(model_trauma_bi_smlh_msat),
                              mcmc_areas_data(model_trauma_bi_smlh_snp),
                              mcmc_areas_data(model_worms_bi_smlh_msat),
                              mcmc_areas_data(model_worms_bi_smlh_snp),
                              mcmc_areas_data(model_malnutrition_bi_smlh_msat),
                              mcmc_areas_data(model_malnutrition_bi_smlh_snp))



brms_primary_bi_area$model <- c(rep("22 microsatellites", nrow(mcmc_areas_data(model_congenital_bi_smlh_msat))),
                                rep("15,051 SNPs", nrow(mcmc_areas_data(model_congenital_bi_smlh_snp))),
                                rep("22 microsatellites", nrow(mcmc_areas_data(model_bacteria_bi_smlh_msat))),
                                rep("15,051 SNPs", nrow(mcmc_areas_data(model_bacteria_bi_smlh_snp))),
                                rep("22 microsatellites", nrow(mcmc_areas_data(model_protozoa_bi_smlh_msat))),
                                rep("15,051 SNPs", nrow(mcmc_areas_data(model_protozoa_bi_smlh_snp))),
                                rep("22 microsatellites", nrow(mcmc_areas_data(model_trauma_bi_smlh_msat))),
                                rep("15,051 SNPs", nrow(mcmc_areas_data(model_trauma_bi_smlh_snp))),
                                rep("22 microsatellites", nrow(mcmc_areas_data(model_worms_bi_smlh_msat))),
                                rep("15,051 SNPs", nrow(mcmc_areas_data(model_worms_bi_smlh_snp))),
                                rep("22 microsatellites", nrow(mcmc_areas_data(model_malnutrition_bi_smlh_msat))),
                                rep("15,051 SNPs", nrow(mcmc_areas_data(model_malnutrition_bi_smlh_snp))))

brms_primary_bi_area$response <- c(rep("Congenital defect", nrow(mcmc_areas_data(model_congenital_bi_smlh_msat))),
                                   rep("Congenital defect", nrow(mcmc_areas_data(model_congenital_bi_smlh_snp))),
                                   rep("Bacteria", nrow(mcmc_areas_data(model_bacteria_bi_smlh_msat))),
                                   rep("Bacteria", nrow(mcmc_areas_data(model_bacteria_bi_smlh_snp))),
                                   rep("Protozoa", nrow(mcmc_areas_data(model_protozoa_bi_smlh_msat))),
                                   rep("Protozoa", nrow(mcmc_areas_data(model_protozoa_bi_smlh_snp))),
                                   rep("Trauma", nrow(mcmc_areas_data(model_trauma_bi_smlh_msat))),
                                   rep("Trauma", nrow(mcmc_areas_data(model_trauma_bi_smlh_snp))),
                                   rep("Worms", nrow(mcmc_areas_data(model_worms_bi_smlh_msat))),
                                   rep("Worms", nrow(mcmc_areas_data(model_worms_bi_smlh_snp))),
                                   rep("Malnutrition", nrow(mcmc_areas_data(model_malnutrition_bi_smlh_msat))),
                                   rep("Malnutrition", nrow(mcmc_areas_data(model_malnutrition_bi_smlh_snp))))

brms_primary_bi_area <- subset(brms_primary_bi_area, grepl("smlh", parameter)) #only select fixed effect beta estimates

brms_primary_bi_interval$model <- factor(as.factor(brms_primary_bi_interval$model), levels = c("22 microsatellites", "15,051 SNPs"))
brms_primary_bi_area$model <- factor(as.factor(brms_primary_bi_area$model), levels = c("22 microsatellites", "15,051 SNPs"))

#plot

# split by interval
brms_primary_bi <- split(brms_primary_bi_area, brms_primary_bi_area$interval)

brms_primary_bi$bottom <- brms_primary_bi$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()


## plot

ggplot(data = brms_primary_bi$outer) + aes(x = .data$x, 
                                           y = fct_relevel(.data$response, "Congenital defect","Malnutrition","Trauma", 
                                                           "Protozoa", "Bacteria", "Worms")) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = .data$response, col = .data$response), size = 0.8)+
  geom_segment(data = brms_primary_bi_interval, aes(x = l, xend = h, yend = response), col = "black", size=3)+
  geom_segment(data = brms_primary_bi_interval, aes(x = ll, xend = hh, yend = response), col = "black")+
  geom_point(data = brms_primary_bi_interval, aes(x = m, y = response), color="black", fill = "grey60", shape=21, size = 4) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  scale_fill_manual(values = alpha(clr_pheno, 0.5))+
  scale_color_manual(values = clr_pheno)+
  facet_wrap(~model, ncol = 3) +
 # xlim(-1, 1)+
  labs(x = "Standardised beta coefficient", y = "Cause of death")+
  theme_bw(base_family = "Arial")+
  theme(text = element_text(size = 24, color = "black"),
        axis.text.y = element_text(color = "black", hjust=0.95),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.3),
        axis.ticks = element_line(colour = "black",
                                  linewidth = 0.3),
        plot.subtitle = element_text(hjust = .5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    color = "black"),
        axis.text.x = element_text(color = "black"),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),
                                    color = "black")) -> plot_1c_binary_primary
plot_1c_binary_primary
ggsave(plot_1c_binary_primary, file = "plots/final/Plot1c_primary_bayes_binary.png", width=8, height = 10)

##### Plot 1d: boxplots of raw data ####

load("data/clean/phenotype_smlh.RData")

#transform to long
data_select <- data[,c("id", "smlh_msat", "smlh_snp", 
                       "malnutrition_bi", "bacteria_bi", "congenital_bi",
                       "worms_bi", "protozoa_bi", "trauma_bi")]
names(data_select) <- c("id", "smlh_msat", "smlh_snp", 
                       "Malnutrition", "Bacteria", "Congenital defect",
                       "Worms", "Protozoa", "Trauma")

data_long <- gather(data_select, condition, measurement, Malnutrition:Trauma)

## facet wrap also smlh method
data_long$smlh_msat <- scale(data_long$smlh_msat) #standardise
data_long$smlh_snp <- scale(data_long$smlh_snp) #standardise
data_long_long <- gather(data_long, method, smlh, smlh_msat:smlh_snp)
data_long_long$method <- as.factor(data_long_long$method)

data_long_long$condition <- fct_relevel(data_long_long$condition,
                                        c("Worms", "Bacteria",
                                          "Protozoa","Trauma", 
                                          "Malnutrition","Congenital defect"))

levels(data_long_long$method) <- c("22 microsatellites", "15,051 SNPs")
data_long_long$method <- factor(as.factor(data_long_long$method), levels = c("22 microsatellites", "15,051 SNPs"))

ggplot(subset(data_long_long, !is.na(smlh)), aes(x = as.factor(measurement), y = smlh, 
                      fill = as.factor(method))) + 
  geom_point(position=position_jitter(width = 0.2, height = 0),
             aes(col = as.factor(method)), size = 2) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, lwd = 0.7, aes(fill = as.factor(method))) + 
  scale_fill_manual(values=alpha(clr_method, 0.5))+
  scale_color_manual(values=alpha(clr_method, 0.6))+
  facet_grid(condition~method, scales = "free", switch = "y")+
  labs(x = "Cause of death", y = "z-transformed sMLH")+
  theme_bw(base_family = "Arial")+
  theme(text = element_text(size = 24),
        legend.position = "none",
        strip.text.y.left = element_text(angle = 0, color = "black"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                  color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),
                                    color = "black"),
      axis.line = element_line(colour = "black",
                               linewidth = 0.3),
      plot.margin = margin(1,1,1,1, "cm"),
      strip.placement = "outside",
      strip.text.y = element_text(hjust=0.95)) -> plot_1d_boxplots

plot_1d_boxplots

ggsave(plot_1d_boxplots, file = "plots/final/Plot1d_boxplot_smlh.png", height = 10, width = 8)


# ### Combine 1a - 1d into Figure 1 ####

plot_left <- plot_grid(plot_1a_g2, plot_1c_binary_primary, labels = c("a)", "c)"), ncol = 1, align = "v", axis = "l",
                      label_fontface = "plain",label_size = 28, rel_heights = c(1,2))
plot_right <- plot_grid(plot_1b_blubber_weight, plot_1d_boxplots, labels = c("b)", "d)"), 
                        ncol = 1, label_fontface = "plain",label_size = 28, rel_heights = c(1,2))

plot_grid(plot_left, plot_right, nrow = 1, align = "lv", axis = "lb") -> plot1

ggsave(plot1, file="plots/final/Plot1_v9.png", height = 18, width = 20)

#### Supplementary plot: effect size sMLH vs chromosome length ####
### Part A: worms ####
load(file = "output/models_snp_worms_chromosome_bayes.RData")

#combine
brms_chr_worms_interval <- mcmc_intervals_data(models_snp_worms[[1]][[1]], prob_outer = 0.95, prob = 0.8)
for (i in 2:17){brms_chr_worms_interval <- rbind(brms_chr_worms_interval,
                                                 mcmc_intervals_data(models_snp_worms[[1]][[i]]))}
brms_chr_worms_interval <- subset(brms_chr_worms_interval, grepl("chr", parameter))

brms_chr_worms_areas <- mcmc_areas_data(models_snp_worms[[1]][[1]])
for (i in 2:17){brms_chr_worms_areas <- rbind(brms_chr_worms_areas,
                                              mcmc_areas_data(models_snp_worms[[1]][[i]]))}
brms_chr_worms_areas <- subset(brms_chr_worms_areas, grepl("chr", parameter))

#rename
brms_chr_worms_interval$parameter <- gsub("b_scalesmlh_chr", "", brms_chr_worms_interval$parameter)
brms_chr_worms_areas$parameter <- gsub("b_scalesmlh_chr", "", brms_chr_worms_areas$parameter)

brms_chr_worms_interval$parameter <- as.factor(as.numeric(brms_chr_worms_interval$parameter ))
brms_chr_worms_areas$parameter <- as.factor(as.numeric(brms_chr_worms_areas$parameter ))


# split by interval
data_chr_worms <- split(brms_chr_worms_areas, brms_chr_worms_areas$interval)

data_chr_worms$bottom <- data_chr_worms$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()


### Part B: bacteria ####
load(file = "output/models_snp_bacteria_chromosome_bayes.RData")

#combine
brms_chr_bacteria_interval <- mcmc_intervals_data(models_snp_bacteria[[1]][[1]], prob_outer = 0.95, prob = 0.8)
for (i in 2:17){brms_chr_bacteria_interval <- rbind(brms_chr_bacteria_interval,
                                                    mcmc_intervals_data(models_snp_bacteria[[1]][[i]]))}
brms_chr_bacteria_interval <- subset(brms_chr_bacteria_interval, grepl("chr", parameter))

brms_chr_bacteria_areas <- mcmc_areas_data(models_snp_bacteria[[1]][[1]])
for (i in 2:17){brms_chr_bacteria_areas <- rbind(brms_chr_bacteria_areas,
                                                 mcmc_areas_data(models_snp_bacteria[[1]][[i]]))}
brms_chr_bacteria_areas <- subset(brms_chr_bacteria_areas, grepl("chr", parameter))

#rename
brms_chr_bacteria_interval$parameter <- gsub("b_scalesmlh_chr", "", brms_chr_bacteria_interval$parameter)
brms_chr_bacteria_areas$parameter <- gsub("b_scalesmlh_chr", "", brms_chr_bacteria_areas$parameter)

brms_chr_bacteria_interval$parameter <- as.factor(as.numeric(brms_chr_bacteria_interval$parameter ))
brms_chr_bacteria_areas$parameter <- as.factor(as.numeric(brms_chr_bacteria_areas$parameter ))


# split by interval
data_chr_bacteria <- split(brms_chr_bacteria_areas, brms_chr_bacteria_areas$interval)

data_chr_bacteria$bottom <- data_chr_bacteria$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()


### Combine in plot ####

## get chr sizes
chrsize <- fread("data/raw/elephant_seal_sequence_report.tsv")
chrsize <- chrsize %>% select(c("Chromosome name", "Seq length"))
names(chrsize) <- c("chr", "size")

## add column phenotype
brms_chr_worms_interval$trait <- "Worms"
brms_chr_bacteria_interval$trait <- "Bacteria"

## combine the two in a new df
brms_chr_interval_combined <- rbind(brms_chr_worms_interval, brms_chr_bacteria_interval)

## join with chr size
brms_chr_interval_combined <- left_join(brms_chr_interval_combined, chrsize[,c("chr", "size")], by = c("parameter" = "chr"))

## convert to Mb
brms_chr_interval_combined$size_mb <- brms_chr_interval_combined$size/1000000

## relevent to get the right order for plotting
data_long_long$condition <- fct_relevel(data_long_long$condition,
                                        c("Worms", "Bacteria",
                                          "Protozoa","Trauma", 
                                          "Malnutrition","Congenital defect"))

brms_chr_interval_combined$trait <- fct_relevel(brms_chr_interval_combined$trait ,
                                                c("Worms", "Bacteria"))
# get R2 for MS
chrsize_effect_worms <- lm(m ~ size, data = subset(brms_chr_interval_combined, trait == "Worms"))
r2(chrsize_effect_worms)

chrsize_effect_bacteria <- lm(m ~ size, data = subset(brms_chr_interval_combined, trait == "Bacteria"))
r2(chrsize_effect_bacteria)

#formula to display statistics on the plot
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
# edited by https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
# further adapted

lm_stat <- function(df){
  m <- lm(m ~ size, df);
  eq <- substitute(~~italic(r)^2~"="~r2*","~~italic("p")~~"="~p, 
                   list(r2 = format(summary(m)$r.squared, digits = 2),
                        p = format(summary(m)$coefficients[2,4], digits = 2)))
  as.character(as.expression(eq));
}

text_worms <- lm_stat(subset(brms_chr_interval_combined, trait == "Worms"))
text_bacteria <- lm_stat(subset(brms_chr_interval_combined, trait == "Bacteria"))

ann_text <- data.frame(size_mb = c(mean(brms_chr_interval_combined$size_mb), mean(brms_chr_interval_combined$size_mb)),
                       m = c(1.5, 1.5),
                       lab = c(text_worms, text_bacteria),
                       trait = factor(c("Worms", "Bacteria")))
# plot

ggplot(brms_chr_interval_combined, aes(x = size_mb, y = m))+
  geom_errorbar(aes(x = size_mb, ymin = ll, ymax = hh), col = "#777777", lwd = 0.7) + #80% CI
  geom_point(aes(col = trait, fill = trait), size = 3) + ylim(-2,2)+theme_classic() +
  
  stat_smooth(method = "lm", aes(col = trait), se = T) +
  labs(x = "Chromomome size (Mb)", y="Standardised beta coefficient")+
  theme_bw(base_family = "Arial")+
  scale_color_manual(values = c(clr_pheno["Worms"], clr_pheno["Bacteria"]))+
  scale_fill_manual(values = c(clr_pheno["Worms"], clr_pheno["Bacteria"])) +
  theme(text = element_text(size = 24, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.3),
        axis.ticks = element_line(colour = "black",
                                  linewidth = 0.3),
        plot.subtitle = element_text(hjust = .5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    color = "black"),
        axis.text.x = element_text(color = "black"),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),
                                    color = "black"))+
  facet_wrap(~trait) +
#  geom_text(data=ann_text_a, label = ann_text_a$lab, parse=T, size = 5) +
  geom_text(data=ann_text, label = ann_text$lab, parse=T, size = 5) -> plot_chrsize
plot_chrsize
ggsave(plot_chrsize, filename = "plots/final/Supp_chrsize_wormsbacteria.png", height = 10, width = 14)



##### Supplementary plot: Primary classification bayes categorical #####
### load data ###
load(file = "output/brms_primaryclass_smlh_msat.RData")
load(file = "output/brms_primaryclass_smlh_snp.RData")

### combine results in one df
brms_primary_interval <- rbind(mcmc_intervals_data(brm_primary_msat, prob =0.8, prob_outer = 0.95),
                               mcmc_intervals_data(brm_primary_snp, prob =0.8, prob_outer = 0.95))

brms_primary_interval$model <- c(rep("22 microsatellites", nrow(mcmc_intervals_data(brm_primary_msat))),
                                 rep("15,051 SNPs", nrow(mcmc_intervals_data(brm_primary_snp))))

brms_primary_interval <- subset(brms_primary_interval, grepl("smlh", parameter))

brms_primary_areas <- rbind(mcmc_areas_data(brm_primary_msat),
                            mcmc_areas_data(brm_primary_snp))

brms_primary_areas$model <- c(rep("22 microsatellites", nrow(mcmc_areas_data(brm_primary_msat))),
                              rep("15,051 SNPs", nrow(mcmc_areas_data(brm_primary_snp))))
brms_primary_areas <- subset(brms_primary_areas, grepl("smlh", parameter)) #only select fixed effect beta estimates

## rename
brms_primary_interval$parameter <- gsub("b_mu", "", brms_primary_interval$parameter)
brms_primary_interval$parameter <- gsub("_scalesmlh_msat", "", brms_primary_interval$parameter)
brms_primary_interval$parameter <- gsub("_scalesmlh_snp", "", brms_primary_interval$parameter)

brms_primary_interval$parameter <- gsub("Bacterialinfection", "Bacteria", brms_primary_interval$parameter)
brms_primary_interval$parameter <- gsub("Congenitaldefect", "Congenital defect", brms_primary_interval$parameter)
brms_primary_interval$parameter <- gsub("Otostrongylis", "Worms", brms_primary_interval$parameter)


brms_primary_areas$parameter <- gsub("b_mu", "", brms_primary_areas$parameter)
brms_primary_areas$parameter <- gsub("_scalesmlh_msat", "", brms_primary_areas$parameter)
brms_primary_areas$parameter <- gsub("_scalesmlh_snp", "", brms_primary_areas$parameter)

brms_primary_areas$parameter <- gsub("Bacterialinfection", "Bacteria", brms_primary_areas$parameter)
brms_primary_areas$parameter <- gsub("Congenitaldefect", "Congenital defect", brms_primary_areas$parameter)
brms_primary_areas$parameter <- gsub("Otostrongylis", "Worms", brms_primary_areas$parameter)

brms_primary_areas$model <- factor(as.factor(brms_primary_areas$model), levels = c("22 microsatellites", "15,051 SNPs"))
brms_primary_interval$model <- factor(as.factor(brms_primary_interval$model), levels = c("22 microsatellites", "15,051 SNPs"))

# split by interval
data_primary <- split(brms_primary_areas, brms_primary_areas$interval)

data_primary$bottom <- data_primary$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = data_primary$outer) + aes(x = .data$x, y = fct_relevel(.data$parameter, "Congenital defect","Malnutrition",
                                                                         "Protozoa", "Bacteria", "Worms")) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = .data$parameter, col = .data$parameter), size = 0.8)+
  geom_segment(data = brms_primary_interval, aes(x = l, xend = h, yend = parameter), col = "black", size=3)+
  geom_segment(data = brms_primary_interval, aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data = brms_primary_interval, aes(x = m, y = parameter), color="black", fill = "grey60", shape=21, size = 4) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  scale_fill_manual(values = alpha(clr_pheno, 0.5))+
  scale_color_manual(values = clr_pheno)+
  facet_wrap(~model, ncol = 2) +
  xlim(-4, 4)+
  labs(x = "Standardised beta coefficient", y = "Cause of death")+
  theme_bw(base_family = "Arial")+
  theme(text = element_text(size = 24, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black",
                                 linewidth = 0.3),
        axis.ticks = element_line(colour = "black",
                                  linewidth = 0.3),
        plot.subtitle = element_text(hjust = .5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    color = "black"),
        axis.text.x = element_text(color = "black"),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),
                                    color = "black")) -> plot_1b_categorical_primary

ggsave(plot_1b_categorical_primary, filename = "plots/final/Supp_primary_bayes_categorical.png", width=10, height=10)

#### Supplementary plot: Power test #####
## plotting

load(file = "output/powertest_all.RData")

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

ggsave(plot_power, file = "plots/final/Supp_power.png", width = 12, height=10)
