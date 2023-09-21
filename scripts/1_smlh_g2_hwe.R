#### Cleaning microsatellite data from excel sheet ####

if (!require("pacman")) install.packages("pacman") #pacman can be used to load and/or install required packages for each script

pacman::p_load(readxl, inbreedR, tidyverse, adegenet, pegas)

#### load in data ####

msat_raw <- read_excel("data/raw/NES dataset 19_04_23.xlsx", sheet = 2, skip = 2) #raw microsatellites, unfiltered

#### Testing for HWE ####
#### Transform MSAT file to structure format ####
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

#### Calculate HWE ####

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
msat_hwe$locus <- names(msat_raw_stru)[c(3:26)] #add names to hwe output
msat_hwe$qval <- p.adjust(msat_hwe$Pr.exact, method ="fdr", n = length(msat_hwe$Pr.exact))

# remove two loci out of HWE: hl.10 and ag9
msat_stru_clean <- msat_raw_stru %>% select(-c("hl.10", "ag9"))

write.table(msat_stru_clean, "data/clean/msat_clean.stru", sep=" ", row.names = F, quote=F)

#### Calculate sMLH based on microsatellites ####
msat_clean <- read_excel("data/raw/NES dataset 19_04_23.xlsx", sheet = 5, skip = 2) #datafile with HWE filtered loci in one row per ID format needed for inbreedR
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
exclude <- c("ES2937", "ES2636", "ES2932","FaMa012912","ES3196","FaMa010712","FaMa010412")
`%!in%` = Negate(`%in%`)
sMLH_df <- subset(sMLH_df, id %!in% exclude)

write.table(sMLH_df, file = "data/smlh/smlh_genomewide_msats.txt", quote=F, row.names=F, sep="\t")

### Calculate g2 ####
#calculate g2 based on 10,000 bootstraps and 100 permutations
g2_msats_10000 <- g2_microsats(msat_clean.inb, nperm = 100, nboot = 10000, CI = 0.95)

write.table(g2_msats_10000$g2_boot,"data/smlh/g2msats_10000.txt",quote=F,col.names=F,row.names=F,sep="\t")

################ SNPs

raw <- read.table("data/raw/NES_notImp.raw",h=T) # Genotype data in plink raw format
#out<-read.table("data/smlh/smlh_snp.txt",h=T) # This is the file where sMLH info will be saved (Empty 74 x 20 matrix with 
# one column containing samples IDs and one column containing cause of death).

# Rename 'raw' IDs to match those of 'out'
iid<-raw$IID
iid<-gsub('-','',iid)
iid<-gsub('_.*','',iid)
raw$IID<-iid

# Calculate genome-wide sMLH
geno<-raw[,-c(1:6)] # Remove columns with no genotypic data
geno[geno==2]<-0
smlh<-sMLH(geno)

out <- data.frame(ID= raw$IID)
# Put it onto the table
temp1<-data.frame(ID=raw$IID,sMLH=smlh)
temp2<-data.frame(ID=out$ID)
temp3<-left_join(temp2,temp1,by="ID")
out$sMLH_GenomeWide <- temp3$sMLH

### sMLH separately for each chromosome
cols<-colnames(geno)
cols<-gsub('.[0-9]*_[ACTG]','',cols) # Keep only chr name
cols[12778:15051]<-"Unplaced"
fac<-factor(cols, levels = c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5",
                             "HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12",
                             "HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","Unplaced"))
index<-table(fac) 

for (i in 1:length(index)){
  if (i==1){
    temp<-geno[,1:index[i]]
    assign(paste("chr_",i,sep=""),sMLH(temp))
  } else {
    temp<-geno[,(((index[i-1]+sum(index[1:(i-2)]))+1):(index[i]+sum(index[1:(i-1)])))]
    assign(paste("chr_",i,sep=""),sMLH(temp))
  }
}

out$sMLH_Chr1 <- chr_1
out$sMLH_Chr2 <- chr_2
out$sMLH_Chr3 <- chr_3
out$sMLH_Chr4 <- chr_4
out$sMLH_Chr5 <- chr_5
out$sMLH_Chr6 <- chr_6
out$sMLH_Chr7 <- chr_7
out$sMLH_Chr8 <- chr_8
out$sMLH_Chr9 <- chr_9
out$sMLH_Chr10 <- chr_10
out$sMLH_Chr11 <- chr_11
out$sMLH_Chr12 <- chr_12
out$sMLH_Chr13 <- chr_13
out$sMLH_Chr14 <- chr_14
out$sMLH_Chr15 <- chr_15
out$sMLH_Chr16 <- chr_16
out$sMLH_Chr17 <- chr_17

write.table(out,"data/smlh/smlh_snp.txt",quote=F,col.names=T,row.names=F,sep="\t")

