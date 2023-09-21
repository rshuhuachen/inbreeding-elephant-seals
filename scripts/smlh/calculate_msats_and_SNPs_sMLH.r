######################################################
### Calculate all sMLH on microsatellites and SNPs ###
######################################################

pacman::p_load(dplyr, inbreedR, tidyr)

################ Microsatellites

msats<-read.table("raw/msats_genotypes.txt",h=F,row.names=1)
msats[msats==0]<-NA
outm<-read.table("raw/msats_dataset_description.txt",h=T)
outm$Animal_ID <- gsub("-","",outm$Animal_ID)

#reformat into inbreedR format
msats<-convert_raw(msats)

#calculate sMLH
m_smlh<-sMLH(msats)
m_smlh<-data.frame(Animal_ID=row.names(msats),sMLH=m_smlh)
out_m<-left_join(outm,m_smlh,by="Animal_ID")

write.table(out_m,"data/smlh/smlh_genowide_msats.txt",quote=F,col.names=T,row.names=F,sep="\t")

#calculate g2 based on 10,000 bootstraps and 100 permutations
g2_msats_10000 <- g2_microsats(msats, nperm = 100, nboot = 10000, CI = 0.95)

write.table(g2_msats_10000$g2_boot,"data/smlh/g2msats_10000.txt",quote=F,col.names=F,row.names=F,sep="\t")

################ SNPs

raw<-read.table("data/raw/NES_notImp.raw",h=T) # Genotype data in plink raw format
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
