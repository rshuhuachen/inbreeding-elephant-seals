#### Comparing David's and Rebecca's way of calculating sMLH

if (!require("pacman")) install.packages("pacman") #pacman can be used to load and/or install required packages for each script

pacman::p_load(readxl, inbreedR, tidyverse, adegenet, pegas)

################ David's way ####

msats<-read.table("data/raw/msats_genotypes.txt",h=F,row.names=1) # N = 223 individuals
msats[msats==0]<-NA
outm<-read.table("data/raw/msats_dataset_description.txt",h=T) # N = 219 individuals
outm$Animal_ID <- gsub("-","",outm$Animal_ID)

# look at differences animal IDs msats and outm
`%!in%` = Negate(`%in%`)
subset(msats, rownames(msats) %!in% outm$Animal_ID) %>% rownames() #9 that are in msats but not outm
subset(outm, outm$Animal_ID %!in% rownames(msats)) %>% select("Animal_ID") #5 that are in outm but not msats

#reformat into inbreedR format
msats<-convert_raw(msats)

#calculate sMLH
m_smlh<-sMLH(msats)
m_smlh<-data.frame(Animal_ID=row.names(msats),sMLH=m_smlh)
out_m<-left_join(outm,m_smlh,by="Animal_ID")

#### Rebecca's way ####
msat_clean <- read_excel("data/raw/NES dataset 19_04_23.xlsx", sheet = 5, skip = 2) #datafile with HWE filtered loci
names(msat_clean)<- c("id", "hg3.6_a", "hg3.6_b", 
                      "pv9_a", "pv9_b",
                      "mang27_a", "mang27_b",
                      "g04_a", "g04_b",  "mang01_a", "mang01_b","dz441_a",  "dz441_b","hi8_a", "hi8_b" , "lw20_a",  
                      "lw20_b","mang44_a","mang44_b","f07_a" ,"f07_b", "pvc1_a" , "pvc1_b","c01_a", "c01_b" ,"mang35_a",
                      "mang35_b" ,"a12_a","a12_b", "e04_a", "e04_b" ,"mang06_a" ,"mang06_b" ,"dh4.7_a", "dh4.7_b", "dh3.6_a", 
                      "dh3.6_b" ,"dh1.8_a",  "dh1.8_b","m11a_a" ,"m11a_b",  "bg_a",  "bg_b", "mang36_a", "mang36_b" )


## Change formats to be loaded into inbreedR (binary where rownames are id's and no populations)

msat_clean.inb <- msat_clean
msat_clean.inb <- column_to_rownames(msat_clean.inb, var = "id")
msat_clean.inb[msat_clean.inb==0] <- NA

# convert to inbreedR
msat_clean.inb <- convert_raw(msat_clean.inb)

#### Calculate sMLH ####

sMLH <- sMLH(msat_clean.inb) #sMLH
## put in df
sMLH_df <- as.data.frame(sMLH)
colnames(sMLH_df)[1] <- "sMLH"
sMLH_df <- tibble::rownames_to_column(sMLH_df, "id")

#exclude some individuals due to low phenotypic uncertainty
exclude <- c("ES2937", "ES2636","ES2695","ES2932","ES3256","FaMa012912","ES3196","FaMa010712","FaMa010412")
`%!in%` = Negate(`%in%`)
sMLH_df <- subset(sMLH_df, id %!in% exclude)

write.table(sMLH_df, file = "data/smlh/smlh_genomewide_msats_2.txt", quote=F, row.names=F, sep="\t")

## import what was sent to me
sent <- fread("data/smlh/smlh_genomewide_msats.txt")

#merge to compare
merge <- left_join(sMLH_df, out_m[,c("Animal_ID", "sMLH")], by = c("id" = "Animal_ID"))
merge <- left_join(merge, sent[,c("Animal_ID", "sMLH")], by = c("id" = "Animal_ID"))
View(merge)
