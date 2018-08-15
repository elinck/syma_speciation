### principal component analysis of multispecies vcf

setwd("~/Dropbox/syma_speciation")
library(vcfR);
library(adegenet);
library(ggplot2);
library(data.table);
library(pcadapt);
library(PopGenome);
library(magrittr)
library(radiator)
library(jaatha);library(coala)
library(ape)
library(treeio)
library(ggtree)
library(phangorn)

# read in sample data
sample_data <- read.csv("data/syma_spp_master.csv")
attach(sample_data)

# 75% complete matrix
syma.med <- read.vcfR("raw_data/syma.75.1perlocus.vcf")
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
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust$clust.k4))+geom_text(aes(label=prep_ID))

# to do -- regress against missing data, batch for all PCs

# load elevation data
png <- raster("raw_data/gis/PNG_alt.grd")
idn <- raster("raw_data/gis/IDN_alt.grd")
aus <- raster("raw_data/gis/AUS_alt.grd")
pac <- raster::merge(png,idn,aus)
pts <- as.matrix(cbind(sample_data$long, sample_data$lat))
pts <- SpatialPoints(pts)
elev <- extract(pac, pts)
sample_data$elevation <- elev

# drop samples with high MD
sample_data <- sample_data[!sample_data$prep_ID=="EL58_toro" & !sample_data$prep_ID=="EL56_toro",]

# merge PCs w/ sample data, write to file
new_df <- merge(sample_data, pc, by.x = "prep_ID", by.y = "prep_ID")
write.csv(new_df, "data/syma_spp_pcs.csv")

# write file for structure w/ radiator
snps <- tidy_genomic_data("raw_data/syma.75.1perlocus.vcf")
genomic_converter(snps, output = c("structure"))

# write file for dadi w/ radiator
syma.dadi <- read.vcfR("raw_data/syma.wgs.filtered.recode.vcf")
vcf2dadi("raw_data/syma.wgs.filtered.recode.vcf", strata = "data/id_pop_dadi.txt", imputation.method = "max")

# read dadi results
dadi <- read.csv("raw_data/dadi_results.txt", sep = '\t')

# subset to last round
dadi_optim <- dadi[grep("Round_3", dadi$Replicate),]
dadi_optim$AIC <- as.numeric(as.character(dadi_optim$AIC))
dadi_optim <- dadi_optim[order(dadi_optim$AIC),] # SC Round_3_Replicate_2 -384.94 781.88, then IMG Round_3_Replicate_3 -421.64 855.28

# create liklihood df
write.csv(dadi_optim, "data/dadi_plotting.csv")
