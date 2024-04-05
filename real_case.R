# http://ricediversity.org/data/sets_hydra/HDRA-G6-4-RDP1-RDP2-NIAS-2.tar.gz
# http://ricediversity.org/data/sets_hydra/HDRA-G6-4-SNP-MAP.tar.gz
# http://ricediversity.org/data/sets_hydra/phenoAvLen_G6_4_RDP12_ALL.tar.gz

#please replace all the pathway in this script to the pathway saved your data
setwd("~/Downloads")

#create a new map file for PLINK binary format to replace allele from A/B to G/C/A/T
map=read.table("HDRA-G6-4-SNP-MAP/HDRA-G6-4-final-snp-map.tsv",head=F)
bim=read.table("HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.bim",head=F)
bim[,5:6]=map[,4:5]
write.table(bim,"HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.bim",quote=F,row.names = F,col.names = F)

#convert PLINK format to VCF format
system("./plink2 --bfile HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS --export vcf --out rice")

#run beagle5 to impute missing genotype. beagle5 could be downloaded from https://faculty.washington.edu/browning/beagle/beagle.html
system("java -Xmx20000m -jar beagle.jar gt=rice.vcf out=test")

#convert imputed data to sample major PLINK text format
system("./plink2 --vcf test.vcf.gz --export A --out rice")

#create txt file
zlines <- readLines("HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.fam")
zline <- "FID IID father mother sex phenotype"
zlines <- c(zline,zlines)
writeLines(zlines,"rice.txt")

system("./blink --file test --vcf --compress --out rice")
system("./blink --file rice --hapmap --recode --out rice")
system("awk '{print $2}' rice.raw | awk -F '_' '{print $1}' > ID.txt")

#load SNP data, map, and phenotype
myGM=read.table("HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.bim",head=F)[,c(2,1,4)]
trait1=read.table("phenoAvLen_G6_4_RDP12_ALL/phenoAvLen_G6_4_RDP12_ALL.txt",head=T)

#match and convert sample ID between SNP data and phenotype data
ID=read.table("ID.txt",head=T)
trait1[,1]=paste("IRGC",trait1[,1],sep='')
trait=read.delim("HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv",head=F)
trait=merge(trait,ID,by.x="V2",by.y="IID")
trait=merge(trait,trait1,by.x="V3",by.y="FID")
myG=read.table("rice.hmp",head=T,check.names = F)
trait[,2]=paste(trait[,2],trait[,3],sep='_')
myY=trait[,c(2,12)]

#run blink via GAPIT
source("http://zzlab.net/GAPIT/gapit_functions.txt")
myGAPIT=GAPIT(Y=myY,G=myG,PCA.total=3,model="BLINK")

