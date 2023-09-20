#### Cleaning phenotypic data from excel sheet ####

#### Load packages ####
if (!require("pacman")) install.packages("pacman") #pacman can be used to load and/or install required packages for each script

pacman::p_load(readxl, tidyverse, lubridate, data.table)

#### Load in data ####

class <- read_excel("data/raw/NES dataset 19_04_23.xlsx", sheet = 7, range="A3:V222")

### Clean phenotypes ####

#names
names(class)[1] <- "id"
names(class)[13] <- "primary_class_num"
names(class)[14] <- "primary_class"

#rename levels + reclassify
factors <- c(6,7,14)
class[factors] <- lapply(class[factors], as.factor)
levels(class$Sex) <- c("male", "female")
levels(class$Age) <- c("pup", "pup", "weaner", "yearling")

nums <- c(9:11)
class[nums] <- lapply(class[nums], as.numeric)

### take out hyphen from ID name
class$id <- gsub("-", "", class$id)

### Resort primary and secondary bacterial infection into one category
class$primary_class <- gsub("Primary bacterial infection", "Bacterial infection", class$primary_class)
class$primary_class <- gsub("Secondary bacterial infection", "Bacterial infection", class$primary_class)

## extract month and year of admittance 
class$admit_month <- months(class$`Admit date`)
class$Year <- year(class$`Admit date`)

### Save clean phenotypes
write.csv(class, "data/clean/phenotypes.csv", quote=F, row.names = F)
save(class, file = "data/clean/phenotypes.RData")

### Add sMLH to the phenotypic dataset ####
### load data smlh
smlh_msat <- fread("data/raw/msats_Samples_sMLH.txt")
smlh_snp <- fread("data/raw/74_Samples_sMLH.txt")
smlh_gl <- fread("data/raw/GL_Heterozygosity_dataset.txt")

## clean ID's to same format
smlh_gl$id <- gsub("-", "", smlh_gl$Sample_ID)

## combine data into a single data frame
data <- left_join(class, smlh_msat[,c("Animal_ID", "sMLH")], by = c("id" = "Animal_ID")) %>%
  left_join(smlh_snp[,c("ID", "sMLH_GenomeWide")], by = c("id" = "ID")) %>%
  left_join(smlh_gl[,c("id", "GL_Het")], by = c("id" = "id"))

#rename variables
names(data)[15:22] <- c("malnutrition", "trauma", "primary_bacterial_infect",
                        "secondary_bacterial_infect", "congenital", "protozoa",
                        "otostrongylis", "herpes")
names(data)[25:27] <- c("smlh_msat", "smlh_snp", "smlh_gl")

data$primary_class <- as.factor(data$primary_class)

#merge primary and secondary bacterial infection into one category
data <- data %>% mutate(bacterial = case_when(
  primary_bacterial_infect == 1 ~ 1,
  secondary_bacterial_infect == 1 ~ 1,
  TRUE ~ 0
))

#binary classification of each category

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

# reformat month
data$admit_month <- as.factor(data$admit_month)
data <- data %>% mutate(admit_month_no = case_when(
  admit_month == "April" ~ 4,
  admit_month == "August" ~ 8,
  admit_month == "February" ~ 2,
  admit_month == "January" ~ 1,
  admit_month == "July" ~ 7,
  admit_month == "June" ~ 6,
  admit_month == "March" ~ 3,
  admit_month == "May" ~ 5,
  admit_month == "November" ~ 11
))

#reformat sex
data <- data %>% mutate(sex_bi = case_when(
  Sex == "male" ~ 0,
  Sex == "female" ~ 1
))
save(data, file="data/clean/phenotype_smlh.RData")
write.csv(data, "data/clean/phenotypes_smlh.csv", quote=F, row.names = F)

