### genotypic data analysis for linck et al. 2019

setwd("~/Dropbox/syma_speciation")
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

# read in sample data
sample_data <- read.csv("data/syma_spp_master.csv")
attach(sample_data)

# 95% complete matrix
syma.med <- read.vcfR("raw_data/vcf/syma.recode.vcf")
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
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust$clust.k3))+geom_text(aes(label=prep_ID))

# remove outlier EL50, recalculate
syma <- syma[-26]
seq.scaled <- scaleGen(syma,NA.method="mean",scale=F)
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
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust$clust.k3))+geom_text(aes(label=prep_ID))
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

# read in d-test data; write to file
dstats1 <- as.data.frame(read.table("~/Dropbox/syma_speciation/data/d_d3_och_meg.csv")[101:200,])
dstats2 <- as.data.frame(read.table("~/Dropbox/syma_speciation/data/d_d3_och_tor.csv")[101:200,])
dstats3 <- as.data.frame(read.table("~/Dropbox/syma_speciation/data/d_d3_meg_tor.csv")[1:100,])
colnames(dstats1) <- c("D")
c1 <- cbind.data.frame(dstats1$D, rep("(ochr,mega),toro)", nrow(dstats1)))
colnames(c1) <- c("D", "topology")
colnames(dstats2) <- c("D")
c2 <- cbind.data.frame(dstats2$D, rep("(ochr,toro),mega)", nrow(dstats2)))
colnames(c2) <- c("D", "topology")
colnames(dstats3) <- c("D")
c3 <- cbind.data.frame(dstats3$D, rep("(mega,toro),ochr)", nrow(dstats3)))
colnames(c3) <- c("D", "topology")
dstat.df <- rbind(c1,c2,c3)
write.csv(dstats.df, "data/dstats.csv")

# test significance of d-test results
t.test(c1$D) #t = 877.23, df = 99, p-value < 2.2e-16
t.test(c2$D) #t = 837.36, df = 99, p-value < 2.2e-16
t.test(c3$D) #t = -19.416, df = 99, p-value < 2.2e-16
range(c3$D) #-0.023348408  0.003304606

# read in IM model params, manipulate, write to file
im.params <- fread("raw_data/IM_realparams_boots_v4.txt", sep = '\t')
colnames(im.params) <- c('nMeg_0','nTor_0','nMeg_1','nTor_1',
                         'T0','T1','mMT0','mMT1','ll_model','theta')

# migration rate
mig0 <- cbind.data.frame(im.params$mMT0, rep("m_T0"))
colnames(mig0) <- c("value", "parameter")
mig1 <- cbind.data.frame(im.params$mMT1, rep("m_T1"))
colnames(mig1) <- c("value", "parameter")
mig.df <- rbind.data.frame(mig0,mig1)
write.csv(dstats.df, "data/migration_rate.csv")

# "real" migration rate
pop0 <- (im.params$nMeg_0+im.params$nTor_0)/2
mig0_real <- cbind.data.frame(im.params$mMT0, pop0, rep("m_T0"))
colnames(mig0_real) <- c("value","pop_size","parameter")
pop1 <- (im.params$nMeg_1+im.params$nTor_1)/2
mig1_real <- cbind.data.frame(im.params$mMT1, pop1, rep("m_T1"))
colnames(mig1_real) <- c("value", "pop_size", "parameter")
mig_real.df <- rbind.data.frame(mig0_real,mig1_real)
mig_real.df$Nm <- mig_real.df$value*mig_real.df$pop_size
write.csv(dstats.df, "data/real_migration_rate.csv")

# divergence times
T0 <- cbind.data.frame(im.params$T0, rep("T0"))
colnames(T0) <- c("value", "parameter")
T1 <- cbind.data.frame(im.params$T1, rep("T1"))
colnames(T1) <- c("value", "parameter")
Tsplit <- cbind.data.frame((im.params$T0+im.params$T1), rep("Tsplit"))
colnames(Tsplit) <- c("value", "parameter")
time.df <- rbind.data.frame(T0,T1,Tsplit)
write.csv(dstats.df, "data/divergence_times.csv")

# megarhyncha pop size
meg0 <- cbind.data.frame(im.params$nMeg_0, rep("nMeg_T0"))
colnames(meg0) <- c("value", "parameter")
meg1 <- cbind.data.frame(im.params$nMeg_1, rep("nMeg_T1"))
colnames(meg1) <- c("value", "parameter")
meg.df <- rbind.data.frame(meg0,meg1)
write.csv(meg.df, "data/meg_sizes.csv")

# tototoro pop size
tor0 <- cbind.data.frame(im.params$nTor_0, rep("nTor_T0"))
colnames(tor0) <- c("value", "parameter")
tor1 <- cbind.data.frame(im.params$nTor_1, rep("nTor_T1"))
colnames(tor1) <- c("value", "parameter")
tor.df <- rbind.data.frame(tor0,tor1)
write.csv(tor.df, "data/tor_sizes.csv")

# get parameter values for ms
mean(T0$value) #22301.72
sd(T0$value) #22301.72
mean(T1$value) #22301.72
sd(T1$value) #13686.37
mean(T0$value) + mean(T1$value) #744163.3

mean(mig0_real$value*mig0_real$pop_size) #3.408352
sd(mig0_real$value*mig0_real$pop_size) #0.00106009
mean(mig1_real$value*mig1_real$pop_size) #3.694829
sd(mig1_real$value*mig1_real$pop_size) #0.01589379

mean(meg0$value) #780894.6
sd(meg0$value) #37243.88
mean(meg1$value) #346054.4
sd(meg1$value) #16517.94

mean(tor0$value) #256736.6
sd(tor0$value) #12254.61
mean(tor1$value) #596653.6
sd(tor1$value) #28464.24

# read in windowed summary stat data, overlapping 
win.df <- fread("raw_data/syma_windows_overlap.csv")
win.df <- win.df[complete.cases(win.df),]
#df <- subset(df, sites > 10)
#restrict to good scaffolds
scaf <- as.data.frame(table(win.df$scaffold))
scaf <- subset(scaf, Freq > 10)
scaf <- as.vector(scaf$Var1)
win.df <- win.df[win.df$scaffold %in% scaf,]
win.df$windowID <- seq.int(nrow(win.df))

#summarize mummer output
files <- list.files("raw_data/mummer_out/coords",full.names = T)
files <- files[grepl("\\.coords",files)]
i <- 1
for(f in files){
  if(i==1){
    dat <- fread(f,sep=" ",data.table=F)
    i=i+1
  } else {
    tmp <- fread(f,sep=" ",data.table=F)
    dat <- rbind(tmp,dat)
    i=i+1
  }
}
colnames(dat) <- c("refStart","refStop","sep1","qStart","qStop","sep2","refLength","qLength","sep3",
                   "p.identity","sep4","names")
dat$refName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[1])
dat$qName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[2])
dat <- arrange(dat,refName,refStart)
sum <- ddply(dat,
             .(qName,refName),
             summarize,
             totalMatch=sum(qLength),
             refStart=min(refStart))
sum <- arrange(sum,refName,refStart,totalMatch)
sum <- subset(sum,totalMatch>5000)                                                 
sum <- ddply(sum,.(qName),function(e){                                             #get contigs with 1 hit > 1000bp
  a <- subset(e,refName==e$refName[e$totalMatch==max(totalMatch)])
  if(nrow(a)==1){
    a
  }
})

#chromosome order for pretty plots
chr_order <- c("1","1A","1B","2","3","4","4A",as.character(5:28),"Z","M","NA")

#merge mummer info with angsd windowed Fst's
win.df <- merge(win.df,sum,by.x="scaffold",by.y="qName",all.x=T,all.y=F)
win.df$chr <- gsub("chr","",win.df$refName)
win.df$chr[!win.df$chr %in% chr_order] <- "NA"
win.df$chr <- factor(win.df$chr,levels=chr_order)
win.df <- arrange(win.df,chr,refStart)
win.df$row <- 1:nrow(win.df)

# second dataset w/ chr labels
chr_labels <- ddply(win.df,.(chr),summarize,mid=median(row),start=min(row),stop=max(row))
chr_labels$chr <- as.character(chr_labels$chr)
chr_labels$chr[chr_labels$chr %in% as.character(21:27)] <- "21-27"
chr_labels$mid[chr_labels$chr=="21-27"] <- median(chr_labels$mid[chr_labels$chr=="21-27"],na.rm=T)
chr_labels$start[chr_labels$chr=="21-27"] <- min(chr_labels$start[chr_labels$chr=="21-27"],na.rm=T)
chr_labels$stop[chr_labels$chr=="21-27"] <- max(chr_labels$stop[chr_labels$chr=="21-27"],na.rm=T)
chr_labels$start[chr_labels$chr=="1B"] <- 0
chr_labels$stop[chr_labels$chr=="1B"] <- 0
chr_labels$win.df <- min(win.df$Fst_meg_toro)-.25*min(win.df$Fst_meg_toro)
chr_labels <- subset(chr_labels,!is.na(chr) & !duplicated(chr))

# ready for faceting
fst <- win.df
win.df$windowID <- seq.int(nrow(win.df))
colnames(win.df) <- c("scaffold","start","end","mid","sites","pi_meg","pi_toro","dxy_meg_toro","fst","windowID", "refName","totalMatch","refStart","chr","row")
fst.df <- cbind.data.frame(win.df$windowID, win.df$start, win.df$end, win.df$chr, win.df$scaffold, win.df$fst,win.df$row)
dxy.df <- cbind.data.frame(win.df$windowID, win.df$start, win.df$end, win.df$chr, win.df$scaffold, win.df$dxy_meg_toro,win.df$row)
pi.tor.df <- cbind.data.frame(win.df$windowID, win.df$start, win.df$end, win.df$chr, win.df$scaffold, win.df$pi_toro,win.df$row)
pi.meg.df <- cbind.data.frame(win.df$windowID, win.df$start, win.df$end, win.df$chr, win.df$scaffold, win.df$pi_meg,win.df$row)
colnames(fst.df) <- c("windowID", "start", "end", "chr", "scaffold", "value","row")
fst.df$stat <- rep("fst", nrow(fst.df))
fst.df$rollmean <- zoo::rollmean(fst.df$value,150,na.pad = TRUE)
colnames(dxy.df) <- c("windowID", "start", "end", "chr", "scaffold", "value","row")
dxy.df$stat <- rep("dxy", nrow(dxy.df))
dxy.df$rollmean <- zoo::rollmean(dxy.df$value,150,na.pad = TRUE)
colnames(pi.tor.df) <- c("windowID", "start", "end", "chr", "scaffold", "value","row")
pi.tor.df$stat <- rep("pi_tor", nrow(pi.tor.df))
pi.tor.df$rollmean <- zoo::rollmean(pi.tor.df$value,150,na.pad = TRUE)
colnames(pi.meg.df) <- c("windowID", "start", "end", "chr", "scaffold", "value","row")
pi.meg.df$stat <- rep("pi_meg", nrow(pi.meg.df))
pi.meg.df$rollmean <- zoo::rollmean(pi.meg.df$value,150,na.pad = TRUE)
win.df <- rbind.data.frame(fst.df, dxy.df, pi.meg.df, pi.tor.df)
win.df$col = ifelse(as.numeric(win.df$chr) %% 2 == 0, 0, 1)
fst_df <- win.df[win.df$stat=="fst",]
dxy_df <- win.df[win.df$stat=="dxy",]
pi_tor_df <- win.df[win.df$stat=="pi_tor",]
pi_meg_df <- win.df[win.df$stat=="pi_meg",]
fst_df$outlier <- fst_df$value>=quantile(fst_df$value,0.995, na.rm = TRUE)
quantile(fst_df$value,0.995, na.rm = TRUE) #0.4084545
dxy_df$outlier <- dxy_df$value>=quantile(dxy_df$value,0.995, na.rm = TRUE)
quantile(dxy_df$value,0.995, na.rm = TRUE) #0.17  
pi_tor_df$outlier <- pi_tor_df$value>=quantile(pi_tor_df$value,0.995, na.rm = TRUE)
quantile(pi_tor_df$value,0.995, na.rm = TRUE) #0.168
pi_meg_df$outlier <- pi_meg_df$value>=quantile(pi_meg_df$value,0.995, na.rm = TRUE)
quantile(pi_meg_df$value,0.995, na.rm = TRUE) #0.1514545
plot.df <- rbind(fst_df,dxy_df,pi_tor_df,pi_meg_df)

# parsing for special characters and labels
plot.df$stat <- factor(plot.df$stat, levels = c("dxy","fst","pi_meg","pi_tor"),
                  ordered = TRUE, labels=c(expression('D[XY]'), expression('F[ST]'),paste(expression(pi),expression('[megarhyncha]')),paste(expression(pi),expression('[torotoro]'))))
#plot.df.seg <- cbind(plot.df, chr_labels)

#subset for revision
plot.df1 <- subset(plot.df, stat=="F[ST]")
plot.df2 <- subset(plot.df, stat=="D[XY]")
plot.df <- rbind.data.frame(plot.df1,plot.df2)
plot.df$chr <- factor(plot.df$chr, levels=chr_order)
plot.df$stat <- factor(plot.df$stat, levels=c("F[ST]","D[XY]"))
write.csv(plot.df,"data/window_stats_chr.csv")

# write chr data file
write.csv(chr_labels,"data/chr_labels.csv")

#subset for chromosome specific plots, FST
win.5 <- subset(plot.df, chr==5 | chr==11 | chr==9 | chr==23)
win.5$chr_order <- factor(win.5$chr, levels=c('5','9','11','23'))
win.5 <- subset(win.5, stat=="F[ST]")
win.5$fst_rollmean <- zoo::rollmean(win.5$value,100,na.pad = TRUE)
ylab.text = expression('F'['ST'])
write.csv(win.5,"data/chromsomes_fst.csv")

#subset for chromosome specific plots, dxy
win.6 <- subset(ploit.df, chr==5 | chr==11 | chr==9 | chr==23)
win.6$chr_order <- factor(win.6$chr, levels=c('5','9','11','23'))
win.6 <- subset(win.6, stat=="D[XY]")
ylab.text = expression('D'['XY'])
win.6$dxy_rollmean <- zoo::rollmean(win.6$value,100,na.pad = TRUE)
write.csv(win.6,"data/chromsomes_dxy.csv")

#subset for correlation plot
df2 <- read.csv("data/window_correlations.csv")[,-1]
df3 <- df2[which(df2$windowID %in% win.5$windowID),]
df3$fst_outlier <- df3$Fst_meg_toro>=quantile(df3$Fst_meg_toro,0.995, na.rm = TRUE)
df3$dxy_outlier <- df3$dxy_meg_toro>=quantile(df3$dxy_meg_toro,0.995, na.rm = TRUE)
df3$fst_outlier_05 <- df3$Fst_meg_toro>=quantile(df3$Fst_meg_toro,0.95, na.rm = TRUE)
df3$dxy_outlier_05 <- df3$dxy_meg_toro>=quantile(df3$dxy_meg_toro,0.95, na.rm = TRUE)
df3$both_outlier <- df3$dxy_outlier_05==TRUE & df3$fst_outlier_05==TRUE
outlier.df <- df3
write.csv(outlier.df, "data/outlier_correlations.csv")

# correlation in outlier windows
genome.mod <- lm(Fst_meg_toro ~ dxy_meg_toro, df2)
summary(genome.mod)
outlier.mod <- lm(Fst_meg_toro ~ dxy_meg_toro, outlier.df)
summary(outlier.mod)

# windowed summary stats
fst.df <- win.df[win.df$stat=="F[ST]",]
dxy.df <- win.df[win.df$stat=="D[XY]",]
mean.fst <- mean(fst.df$value)
mean.dxy <- mean(dxy.df$value)

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