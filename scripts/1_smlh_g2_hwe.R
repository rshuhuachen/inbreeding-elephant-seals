#### Cleaning microsatellite data from excel sheet ####
#### MESS- TAKE FROM DAVID ###
pacman::p_load(readxl, inbreedR, tidyverse, adegenet, pegas)

#### load in data ####

msat_raw <- read_excel("data/raw/NES dataset 19_04_23.xlsx", sheet = 2, skip = 2)

#### Transform MSAT to structure ####
# raw file is in format one row per ID, only 1/2 msats are named
# for structure format: two rows per ID
# for inbreedR: need this format (one row per ID)

## change format to .stru
names(msat_raw)<- c("id", "hg3.6_a", "hg3.6_b", 
                    "hl.10_a", "hl.10_b", "pv9_a", "pv9_b",
                    "ag9_a", "ag9_b", "mang27_a", "mang27_b",
                    "g04_a", "g04_b",  "mang01_a", "mang01_b","dz441_a",  "dz441_b","hi8_a", "hi8_b" , "lw20_a",  
                    "lw20_b","mang44_a","mang44_b","f07_a" ,"f07_b", "pvc1_a" , "pvc1_b","c01_a", "c01_b" ,"mang35_a",
                    "mang35_b" ,"a12_a","a12_b", "e04_a", "e04_b" ,"mang06_a" ,"mang06_b" ,"dh4.7_a", "dh4.7_b", "dh3.6_a", 
                    "dh3.6_b" ,"dh1.8_a",  "dh1.8_b","m11a_a" ,"m11a_b",  "bg_a",  "bg_b", "mang36_a", "mang36_b" )

msat_raw_ms <- msat_raw[,c(-1)]
msat_raw_ms <- as.data.frame(msat_raw_ms)
msat_raw_a <- cbind(msat_raw[,1],msat_raw_ms[c(T,F)])
msat_raw_b <- cbind(msat_raw[,1],msat_raw_ms[c(F,T)])

#remove suffix
names(msat_raw_a) <- sub("_a", "", names(msat_raw_a), fixed=TRUE)
names(msat_raw_b) <- sub("_b", "", names(msat_raw_a), fixed=TRUE)

#merge
msat_raw_stru <- rbind(msat_raw_a, msat_raw_b)

#add empty col for pop
msat_raw_stru$pop <- NA
msat_raw_stru <- msat_raw_stru[,c(1,ncol(msat_raw_stru), 2:(ncol(msat_raw_stru)-1))]

msat_raw_stru[msat_raw_stru == 0] <- -9

msat_raw_stru <- msat_raw_stru %>% arrange(id)
write.table(msat_raw_stru, "data/clean/msat_raw.stru", sep=" ", row.names = F, quote=F)

#### Filter for HWE ####

#### Summary stats ####
msat_stru_raw <- read.structure("data/clean/msat_raw.stru", n.ind = 223, n.loc = 24, onerowperind = F,
                                col.lab = 1, col.pop = 2, col.others = NULL,
                                row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                ask = F, quiet = FALSE)

summary_msat_raw <- summary(msat_stru_raw)
summary(summary_msat_raw$loc.n.all)

## HWE
msat_hwe <- hw.test(msat_stru_raw, B = 1000) #B = 1000 for 1000 Monte Carlo permutations 
msat_hwe <- as.data.frame(msat_hwe)
msat_hwe$locus <- names(msat_raw)[c(3:26)] #add names to hwe output
msat_hwe$qval <- p.adjust(msat_hwe$Pr.exact, method ="fdr", n = length(msat_hwe$Pr.exact))

# remove two loci
msat_stru_clean <- msat_raw_stru %>% select(-c("hl.10", "ag9"))

write.table(msat_stru_clean, "data/clean/msat_clean.stru", sep=" ", row.names = F, quote=F)


msat_clean <- fread("data/clean/msat_clean.stru")

#### Calculate sMLH ####
msat_clean <- read_excel("data/raw/NES dataset 19_04_23.xlsx", sheet = 5, skip = 2) #datafile with HWE filtered loci
names(msat_clean)<- c("id", "hg3.6_a", "hg3.6_b", 
                    "pv9_a", "pv9_b",
                    "mang27_a", "mang27_b",
                    "g04_a", "g04_b",  "mang01_a", "mang01_b","dz441_a",  "dz441_b","hi8_a", "hi8_b" , "lw20_a",  
                    "lw20_b","mang44_a","mang44_b","f07_a" ,"f07_b", "pvc1_a" , "pvc1_b","c01_a", "c01_b" ,"mang35_a",
                    "mang35_b" ,"a12_a","a12_b", "e04_a", "e04_b" ,"mang06_a" ,"mang06_b" ,"dh4.7_a", "dh4.7_b", "dh3.6_a", 
                    "dh3.6_b" ,"dh1.8_a",  "dh1.8_b","m11a_a" ,"m11a_b",  "bg_a",  "bg_b", "mang36_a", "mang36_b" )

# drop ES2636, ES3196, FN-MA_010412, FN-MA-010712
msat_clean <- subset(msat_clean, id != "ES2636" & id != "ES2695" & id != "ES2932" &
                            id != "ES2937" & id != "ES3196" &
                            id != "ES3256" & id != "FaMa010412" &
                            id != "FaMa010712" & id != "FaMa012912")

## Change formats to be loaded into inbreedR (binary where rownames are id's and no populations)

msat_clean.inb <- msat_clean
rownames(msat_clean.inb) <- msat_clean.inb$id
# convert to inbreedR
msat_clean.inb <- convert_raw(msat_clean.inb)

#### Calculate sMLH ####

sMLH <- sMLH(msat_clean.inb) #sMLH
summary(sMLH)
hist(sMLH)
het_var <- var(sMLH, na.rm=TRUE) # variance in sMLH
summary(het_var)

#### Calculate g2 ####
g2_msat <- g2_microsats(msat_clean.inb, nboot = 10000, nperm = 10000)
g2_msat

## put in df

sMLH_df <- as.data.frame(sMLH)
colnames(sMLH_df)[1] <- "sMLH"
sMLH_df <- tibble::rownames_to_column(sMLH_df, "id")

## write out table
write.csv(sMLH_df, "data/clean/sMLH.csv", row.names = F, quote=F)
