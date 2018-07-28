### principal component analysis of multispecies vcf

setwd("~/Dropbox/syma_speciation")
library(vcfR);library(adegenet);library(ggplot2);
library(ape);library(strataG);library(data.table);
library(pcadapt);library("qvalue");library("OutFLANK");
library("ggplot2");library(vcfR);library(PopGenome)
library(plyr);library(ape);library(ggtree)
library(magrittr)
library(radiator)

# read in sample data
sample_data <- read.csv("data/syma_spp_master.csv")
attach(sample_data)

# 85% complete matrix
syma.med <- read.vcfR("raw_data/syma.gatk.85p.d5.maf05.recode.vcf")
dna <- vcfR2DNAbin(syma.med, unphased_as_NA = F, consensus = T, extract.haps = F)
syma <- DNAbin2genind(dna)
samples <- as.character(read.table("data/samples.txt")[[1]]) 
rownames(syma@tab) <- samples
seq.scaled <- scaleGen(syma,NA.method="mean",scale=F)
#snpID <- as.data.frame(syma.med.genind$loc.n.all)

clust.k1 <- find.clusters(syma,n.pca=95,n.clust = 1,choose.n.clust = F)
clust.k2 <- find.clusters(syma,n.pca=95,n.clust = 2,choose.n.clust = F)
clust.k3 <- find.clusters(syma,n.pca=95,n.clust = 3,choose.n.clust = F)
clust.k4 <- find.clusters(syma,n.pca=95,n.clust = 4,choose.n.clust = F)
clust.k5 <- find.clusters(syma,n.pca=95,n.clust = 5,choose.n.clust = F)
clust.k6 <- find.clusters(syma,n.pca=95,n.clust = 6,choose.n.clust = F)
clust <- cbind(sampleID=rownames(syma@tab),clust.k1=unname(clust.k1$grp),clust.k2=unname(clust.k2$grp),clust.k3=unname(clust.k3$grp),
               clust.k4=unname(clust.k4$grp),clust.k5=unname(clust.k5$grp),clust.k6=unname(clust.k6$grp)) %>% data.frame()
pca <- prcomp(seq.scaled,center=F,scale=F)
screeplot(pca)
pc <- data.frame(pca$x[,1:3])
pc$prep_ID <- rownames(pc)

# quick visualization
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust$clust.k2))+geom_text(aes(label=prepID))

# to do -- regress against missing data, batch for all PCs

# merge PCs w/ sample data, write to file
new_df <- merge(sample_data, pc, by.x = "prep_ID", by.y = "prep_ID")
write.csv(new_df, "data/syma_spp_pcs.csv")

# write file for structure w/ radiator
snps <- tidy_genomic_data("raw_data/syma.85p.1perlocus.vcf")
genomic_converter(snps, output = c("structure"))

# random seeds for structure
sample(1000000:2000000, 15) 
# [1] 1631127 1557097 1646607 1977274 1347526 
# [5] 1529654 1715085 1706929 1971822 1907208 
# [10] 1834841 1136596 1556258 1352111 1354814



