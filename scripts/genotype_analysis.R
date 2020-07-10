### genotypic data analysis for linck et al. 2019

# set working directory
setwd("~/Dropbox/syma_speciation")

# load libraries
library(vcfR)
library(adegenet)
library(ggplot2)
library(data.table)
library(magrittr)
library(radiator)
library(ape)
library(phangorn)
library(raster)
library(ggrepel)
library(plyr)
library(treeio)
library(ggnetworx)
library(ggtree)
library(poppr)

# read in sample data
sample_data <- read.csv("data/syma_spp_master.csv")
attach(sample_data)
syma.med <- read.vcfR("raw_data/vcf/syma.recode.vcf")
dna <- vcfR2DNAbin(syma.med, unphased_as_NA = F, consensus = T, extract.haps = F)
syma <- DNAbin2genind(dna)

# drop all hyrad samples for WGS only PCA
drop <- c(16:24,26:32)
syma.wgs <- syma[-drop]
syma.wgs.scaled <- scaleGen(syma.wgs,NA.method="mean",scale=F)
screeplot(pca.wgs)
pc.wgs <- data.frame(pca.wgs$x[,1:3])
pc.wgs$prep_ID <- rownames(syma.wgs@tab)
sample_drop <- syma[drop]@tab %>% rownames()
sample_wgs <- sample_data[!sample_data$prep_ID %in% sample_drop,]
wgs_pca <- merge(sample_wgs, pc, by.x = "prep_ID", by.y = "prep_ID")
write.csv(wgs_pca, "data/syma_spp_pcs_wgs.csv")
ggplot(data=pc.wgs, aes(x=PC1,y=PC2))+geom_text(aes(label=prep_ID))

# drop all modern samples for aDNA only PCA
drop_mod <- c(16:21)
syma.old <- syma[-drop_mod]
syma.old <- missingno(syma.old, type = "mean", quiet = FALSE, freq = FALSE)
syma.old.scaled <- scaleGen(syma.old,NA.method="zero",scale=F)
pca.old <- prcomp(syma.old.scaled, center=F,scale=F)
screeplot(pca.old)
pc.old <- data.frame(pca.old$x[,1:3])
pc.old$prep_ID <- rownames(syma.old@tab)
sample_drop <- syma[drop_mod]@tab %>% rownames()
sample_old <- sample_data[!sample_data$prep_ID %in% sample_drop,]
old_pca <- merge(sample_old, pc.old, by.x = "prep_ID", by.y = "prep_ID")
write.csv(old_pca, "data/syma_spp_pcs_old.csv")
ggplot(data=pc.old,aes(x=PC1,y=PC2))+geom_text(aes(label=prep_ID))

# 95% complete matrix for fully sampled version in main text
syma <- missingno(syma, type = "mean", quiet = FALSE, freq = FALSE)
clust.k1 <- find.clusters(syma,n.pca=95,n.clust = 1,choose.n.clust = F)
clust.k2 <- find.clusters(syma,n.pca=95,n.clust = 2,choose.n.clust = F)
clust.k3 <- find.clusters(syma,n.pca=95,n.clust = 3,choose.n.clust = F)
clust.k4 <- find.clusters(syma,n.pca=95,n.clust = 4,choose.n.clust = F)
clust.k5 <- find.clusters(syma,n.pca=95,n.clust = 5,choose.n.clust = F)
clust.k6 <- find.clusters(syma,n.pca=95,n.clust = 6,choose.n.clust = F)
clust <- cbind(sampleID=rownames(syma@tab),clust.k1=unname(clust.k1$grp),clust.k2=unname(clust.k2$grp),clust.k3=unname(clust.k3$grp),
               clust.k4=unname(clust.k4$grp),clust.k5=unname(clust.k5$grp),clust.k6=unname(clust.k6$grp)) %>% data.frame()
seq.scaled <- scaleGen(syma, NA.method="zero",scale=F)
pca <- prcomp(seq.scaled,center=F,scale=F)
screeplot(pca)
pc <- data.frame(pca$x[,1:3])
pc$prep_ID <- rownames(pc)

# quick visualization
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust$clust.k3))+geom_text(aes(label=prep_ID))

# best fit clusters
clust.best <- find.clusters(syma,n.pca=95,choose.n.clust=FALSE, criterion = "smoothNgoesup") #k3 chosen
summary(pca) #PC1; 0.1702; PC2; 0.06311

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
sample_data <- sample_data[!sample_data$prep_ID=="EL58_toro" & !sample_data$prep_ID=="EL56_toro" & !sample_data$prep_ID=="EL50_toro" & !sample_data$prep_ID=="EL5_ochr",]

# merge PCs w/ sample data, write to file
new_df <- merge(sample_data, pc, by.x = "prep_ID", by.y = "prep_ID")
write.csv(new_df, "data/syma_spp_pcs.csv")

# load likelihoods
im <- fread("raw_data/IM_realparams.txt", sep = '\t')[,6] %>% max() #-631.8897
sc <- fread("raw_data/SC_realparams.txt", sep = '\t')[,7] %>% max() #-564.062
am <- fread("raw_data/AM_realparams.txt", sep = '\t')[,7] %>% max() #-590.0233
si <- fread("raw_data/SI_realparams.txt", sep = '\t')[,4] %>% max() #-637.185
img <- fread("raw_data/IMg_realparams.txt", sep = '\t')[,6] %>% max() #-617.116
scg <- fread("raw_data/SCg_realparams.txt", sep = '\t')[,7] %>% max() #-716.0587
amg <- fread("raw_data/AMg_realparams.txt", sep = '\t')[,7] %>% max() #-653.5472
sig <- fread("raw_data/SIg_realparams.txt", sep = '\t')[,4] %>% max() #-611.0493

# test for statistical artifacts by coverage etc. 

# load data and merge with missing data info
md <- read.table("raw_data/vcf/out.imiss", sep = "\t")[-1,]
colnames(md) <- c("sample", "loci", "genotypes","filtered","f_miss")
md_df <- merge(new_df, md, by.x="prep_ID", by.y="sample")

# make missing data info numeric
md_df$f_miss <- as.numeric(as.character(md_df$f_miss))

# add ochracea as species category
md_df$sp <- as.character(md_df$sp)
md_df[12,]$sp <- "ochracea"
md_df[23,]$sp <- "ochracea"
md_df$sp <- as.factor(md_df$sp)

# assign modern and historic samples
md_df$tissue_type <- ifelse(md_df$year>2000, "modern", "historic")

# missing data and coverage mean and sd
aggregate(md_df$f_miss, by=list(md_df$sequencing_strategy), mean) #hyRAD 0.014916457, WGS 0.006966324
aggregate(md_df$f_miss, by=list(md_df$sequencing_strategy), sd) #hyRAD 0.01343492, WGS 0.01037236
aggregate(md_df$coverage, by=list(md_df$sequencing_strategy), mean) #hyRAD 4.427339, WGS 6.021934
aggregate(md_df$coverage, by=list(md_df$sequencing_strategy), sd) #hyRAD 2.794949, WGS 3.569487
aggregate(md_df$f_miss, by=list(md_df$tissue_type), mean) #historic 0.011428831, modern 0.007341647
aggregate(md_df$f_miss, by=list(md_df$tissue_type), sd) #historic 0.013013127, modern 0.009132718
aggregate(md_df$coverage, by=list(md_df$tissue_type), mean) #historic 5.452605, modern 4.426548
aggregate(md_df$coverage, by=list(md_df$tissue_type), sd) #historic 3.320965, modern 3.235294


# subset by sequencing strategy
wgs <- md_df[md_df$sequencing_strategy=="WGS",]
hyrad <- md_df[md_df$sequencing_strategy=="hyRAD",]

# pariwise wilcox test for missing data among species
pairwise.wilcox.test(md_df$f_miss, md_df$sp,
                     p.adjust.method = "bonferroni")

# pariwise wilcox test for coverage among species
pairwise.wilcox.test(md_df$coverage, md_df$sp,
                     p.adjust.method = "bonferroni")

# read in SC model params, manipulate, write to file
sc.params <- fread("raw_data/SC_realparams_boots.txt", sep = '\t')
colnames(sc.params) <- c('nTor','nMeg',
                         'T0','T1','m12','m21','ll_model','theta')

# migration rate
mig0 <- cbind.data.frame(sc.params$m12, rep("Nm_tor_meg"))
colnames(mig0) <- c("value", "parameter")
mig1 <- cbind.data.frame(sc.params$m21, rep("Nm_meg_tor"))
colnames(mig1) <- c("value", "parameter")
mig.df <- rbind.data.frame(mig0,mig1)
write.csv(mig.df, "data/migration_rate.csv")

# divergence times
T0 <- cbind.data.frame(sc.params$T0, rep("T0"))
colnames(T0) <- c("value", "parameter")
T1 <- cbind.data.frame(sc.params$T1, rep("T1"))
colnames(T1) <- c("value", "parameter")
Tsplit <- cbind.data.frame((im.params$T0+im.params$T1), rep("Tsplit"))
colnames(Tsplit) <- c("value", "parameter")
time.df <- rbind.data.frame(T0,T1,Tsplit)
write.csv(time.df, "data/divergence_times.csv")

# tototoro pop size
tor.df <- cbind.data.frame(sc.params$nTor, rep("nTor"))
colnames(tor.df) <- c("value", "parameter")
write.csv(tor.df, "data/tor_sizes.csv")

# megarhyncha pop size
meg.df <- cbind.data.frame(sc.params$nMeg, rep("nMeg"))
colnames(meg.df) <- c("value", "parameter")
write.csv(meg.df, "data/meg_sizes.csv")

# get parameter values for ms
mean(T0$value) #458469.1
sd(T0$value) #6162.569
mean(T1$value) #191231.2
sd(T1$value) #7510.339
mean(T0$value) + mean(T1$value) #649700.3

mean(mig0$value) #1.147602
sd(mig0$value) #0.02284373
mean(mig1$value) #0.1190246
sd(mig1$value) #0.006739464

mean(meg.df$value) #368390.6
sd(meg.df$value) #8528.114
mean(tor.df$value) #1909985
sd(tor.df$value) #26028.62

# coverage and mapping summary stats
cov.df <- read.table("data/coverage.txt")
mapping.df <- read.table("data/percent_mapped.txt")
mean(cov.df$V2) #5.388055
range(cov.df$V2) #1.91992 12.11580
mean(mapping.df$V2) #0.8336799
range(mapping.df$V2) #0.355944 0.923737
mapping.df <- subset(mapping.df, !(V1 %in% c('EL5_ochr_R1.clean.fq')))
mean(mapping.df$V2) #0.8459296
range(mapping.df$V2) #0.479000; 0.923737

