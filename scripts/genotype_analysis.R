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

# write file for structure w/ radiator
snps_str <- tidy_genomic_data("raw_data/vcf/syma.str.recode.vcf")
genomic_converter(snps_str, output = c("plink"))

# write file for dadi (moments) w/ radiator
vcf2dadi("raw_data/syma.moments.revised.unlinked.recode.vcf", strata = "raw_data/id_pop_dadi.txt", imputation.method = "max")

# read dadi formatted file, edit for moments bootsrapping, rewrite
dadi <- fread("raw_data/syma.moments.revised.txt", sep = '\t')
dadi$MARKERS <- gsub('scaffold_', 'scaffold', dadi$MARKERS)
dadi$MARKERS <- gsub('__.*__', '_', dadi$MARKERS)
dadi <- as.data.frame(dadi)
write.table(dadi, "raw_data/syma.moments.revised.txt", sep = " ", quote=FALSE, row.names = F)

# read results, export ll_model for plotting, SI
si.moments <- fread("raw_data/moments_results/SI_modelparams.txt", sep = '\t')
colnames(si.moments) <- c('nu1_0','nu2_0','nu1','nu2','T0','T1','theta','Fst','ll_model','L')
si.moments$model <- rep("SI",nrow(si.moments))
si.ll <- cbind.data.frame(si.moments$ll_model, si.moments$model)
colnames(si.ll) <- c("ll_model","model")

# read results, export ll_model for plotting, SC
sc.moments <- fread("raw_data/moments_results/SC_modelparams.txt", sep = '\t')
colnames(sc.moments) <- c('nu1_0','nu2_0','nu1','nu2','T0','T1','theta','m12','m21','Fst','ll_model','L')
sc.moments$model <- rep("SC",nrow(sc.moments))
sc.ll <- cbind.data.frame(sc.moments$ll_model,sc.moments$model)
colnames(sc.ll) <- c("ll_model","model")

# read results, export ll_model for plotting, IM
im.moments <- fread("raw_data/moments_results/IM_modelparams.txt", sep = '\t')
colnames(im.moments) <- c('nu1_0','nu2_0','nu1','nu2','T0','T1','theta','m12','m21','Fst','ll_model','L')
im.moments$model <- rep("IM",nrow(im.moments))
im.ll <- cbind.data.frame(im.moments$ll_model,im.moments$model)
colnames(im.ll) <- c("ll_model","model")

# read results, export ll_model for plotting, IIM
iim.moments <- fread("raw_data/moments_results/IIM_modelparams.txt", sep = '\t')
colnames(iim.moments) <- c('nu1_0','nu2_0','nu1','nu2','T0','T1','theta','m12','m21','Fst','ll_model','L')
iim.moments$model <- rep("IIM",nrow(iim.moments))
iim.ll <- cbind.data.frame(iim.moments$ll_model,iim.moments$model)
colnames(iim.ll) <- c("ll_model","model")

ll.models <- rbind.data.frame(si.ll,sc.ll,im.ll,iim.ll)
write.csv(ll.models, "data/moments_plotting.csv")

# make neighbor joining tree 
consensus.dist <- dist(syma)
consensus <- nj(consensus.dist)
plot(consensus)
nodelabels() # reroot appropriately
consensus <- root(consensus, node = 58)
consensus.smoothed <- chronopl(consensus, lambda = 1)

# make 1000 subsampled SNP matricies
matrix <- list()
for(i in 1:500){
  tmp <- sample(31577, 10000, replace = TRUE)
  matrix[[i]] <- syma[,tmp]
}

dist <- list()
for(i in 1:500){
  dist[[i]] <- dist(matrix[[i]])
}

trees <- list()
for(i in 1:500){
  trees[[i]] <- nj(dist[[i]])
}

# evaluate biparition support
clad <- prop.clades(consensus, trees)
layout(1)
par(mar = rep(2, 4))
plot(consensus, main = "Bipartition vs. Clade Support Values")
drawSupportOnEdges(clad)
legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
       pt.bg = c("green", "lightblue"), pt.cex = 2.5)

test <- trees[[1]]
plot(test)
nodelabels() # reroot appropriately
for(i in 1:500){
  trees[[i]] <- root(trees[[i]], node = 62)
}

# save as multiphylo and export
class(trees) <- "multiPhylo"
writeNexus(trees, "raw_data/densitree.input.nex")

# make and plot network from distance matrix
rownames(syma@tab)
syma_new <- syma[-c(26:31),]
consensus.dist <- dist(syma_new)
write.table(consensus.dist, "data/consensus.dist.tsv")
nnet <- neighborNet(consensus.dist)
# confusing order:
# EL13_toro,EL1_mega,EL11_toro,EL19_mega,EL48_ochr,EL39_toro,EL46_toro,EL4_mega,
# EL10_toro,EL21_toro,EL32_toro,EL9_toro,EL27_mega,EL43_toro,EL44_toro,EL13_toro,
# EL6_mega,EL24_mega,EL47_toro,EL40_mega,EL42_toro,EL19_mega,EL41_toro,EL18_mega,
# EL49_toro,EL23_mega,EL29_ochr,EL20_mega
cols <- c("#443B84FF","#20938CFF","#443B84FF","#443B84FF","#20938CFF","#443B84FF","#443B84FF","#73D056FF",
          "#443B84FF","#443B84FF","#443B84FF","#443B84FF","#73D056FF","#443B84FF","#443B84FF","#443B84FF",
          "#73D056FF","#73D056FF","#443B84FF","#73D056FF","#443B84FF","#73D056FF","#443B84FF","#73D056FF",
          "#443B84FF","#73D056FF","#73D056FF","#73D056FF")
  
pdf("figures/test_phylonet.pdf",width = 7, height = 5)
ggnetworx::ggnetworx(nnet,branch.length = "none",yscale = 0.1,right = TRUE) +
  geom_tiplab(size=3,color=cols)
dev.off()

# read in smcpp plots, save single file for facetting
smc_toro <- read.csv("raw_data/toro.plot.csv",stringsAsFactors = F)
smc_toro$species <- rep("torotoro", nrow(smc_toro))
smc_mega <- read.csv("raw_data/mega.plot.csv",stringsAsFactors = F)
smc_mega$species <- rep("megarhyncha", nrow(smc_mega))
smcpp.df <- rbind.data.frame(smc_toro,smc_mega)
write.csv(smcpp.df, "data/smcpp_df.csv")

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
write.csv(plot.df,"data/window_stats_chr.csv")

# write chr data file
write.csv(chr_labels,"data/chr_labels.csv")

# read in windowed summary stat data, nonoverlapping (for outlier window analysis)
win.df <- fread("raw_data/syma_windows.csv")
win.df <- win.df[complete.cases(win.df),]
#df <- subset(df, sites > 10)

#restrict to good scaffolds
scaf <- as.data.frame(table(win.df$scaffold))
scaf <- subset(scaf, Freq > 10)
scaf <- as.vector(scaf$Var1)
win.df <- win.df[win.df$scaffold %in% scaf,]
win.df$windowID <- seq.int(nrow(win.df))

#merge mummer info with angsd windowed Fst's
win.df <- merge(win.df,sum,by.x="scaffold",by.y="qName",all.x=T,all.y=F)
win.df$chr <- gsub("chr","",win.df$refName)
win.df$chr[!win.df$chr %in% chr_order] <- "NA"
win.df$chr <- factor(win.df$chr,levels=chr_order)
win.df <- arrange(win.df,chr,refStart)
win.df$row <- 1:nrow(win.df)

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
fst.df$rollmean <- zoo::rollmean(fst.df$value,50,na.pad = TRUE)
colnames(dxy.df) <- c("windowID", "start", "end", "chr", "scaffold", "value","row")
dxy.df$stat <- rep("dxy", nrow(dxy.df))
dxy.df$rollmean <- zoo::rollmean(dxy.df$value,50,na.pad = TRUE)
colnames(pi.tor.df) <- c("windowID", "start", "end", "chr", "scaffold", "value","row")
pi.tor.df$stat <- rep("pi_tor", nrow(pi.tor.df))
pi.tor.df$rollmean <- zoo::rollmean(pi.tor.df$value,50,na.pad = TRUE)
colnames(pi.meg.df) <- c("windowID", "start", "end", "chr", "scaffold", "value","row")
pi.meg.df$stat <- rep("pi_meg", nrow(pi.meg.df))
pi.meg.df$rollmean <- zoo::rollmean(pi.meg.df$value,50,na.pad = TRUE)
win.df <- rbind.data.frame(fst.df, dxy.df, pi.meg.df, pi.tor.df)
win.df$col = ifelse(as.numeric(win.df$chr) %% 2 == 0, 0, 1)
fst_df <- win.df[win.df$stat=="fst",]
dxy_df <- win.df[win.df$stat=="dxy",]
pi_tor_df <- win.df[win.df$stat=="pi_tor",]
pi_meg_df <- win.df[win.df$stat=="pi_meg",]
fst_df$outlier <- fst_df$value>=quantile(fst_df$value,0.995, na.rm = TRUE)
quantile(fst_df$value,0.995, na.rm = TRUE) #0.409031 
dxy_df$outlier <- dxy_df$value>=quantile(dxy_df$value,0.995, na.rm = TRUE)
quantile(dxy_df$value,0.995, na.rm = TRUE) #0.1727 
pi_tor_df$outlier <- pi_tor_df$value>=quantile(pi_tor_df$value,0.995, na.rm = TRUE)
quantile(pi_tor_df$value,0.995, na.rm = TRUE) #0.1701
pi_meg_df$outlier <- pi_meg_df$value>=quantile(pi_meg_df$value,0.995, na.rm = TRUE)
quantile(pi_meg_df$value,0.995, na.rm = TRUE) #0.1527655
plot.df <- rbind(fst_df,dxy_df,pi_tor_df,pi_meg_df)

# write csv, reload
write.csv(plot.df,"data/window_stats_outliers.csv")
plot.df <- read.csv("data/window_stats_outliers.csv")

# get Fst outlier windows
fst.out <- plot.df[plot.df$stat=="F[ST]",]
fst.out <- fst.out[fst.out$outlier=="TRUE",]
outliers <- cbind.data.frame(fst.out$start, fst.out$end, fst.out$scaffold, fst.out$chr)
colnames(outliers) <- c("start","end","scaffold","chr")
write.table(outliers, "~/Dropbox/syma_speciation/data/outlier_coordinates.csv",sep=",",row.names = F)
fst.out$scaffold <- droplevels(fst.out$scaffold)
outlier.scaf <- levels(fst.out$scaffold)
for(i in outlier.scaf){ #print list of outliers
  print(paste0("--chr ",i))
}

# see which scaffolds are on chr 5
fst.5 <- subset(fst.out,chr==5)
fst.5$chr <- droplevels(fst.5$chr)
levels(fst.5$scaffold)
fst.5$scaffold 

# get random sample of nonoutliers 
fst.low <- plot.df[plot.df$stat=="",]
fst.low <- fst.low[fst.low$outlier=="FALSE",]
fst.low <- fst.low[sample(nrow(fst.low), 33), ]
fst.low$scaffold <- droplevels(fst.low$scaffold)
low.scaf <- levels(fst.low$scaffold)
for(i in low.scaf){ #print list of outliers
  print(paste0("--chr ",i))
}

#prepare selective sweep df
meg.out <- fread("raw_data/RAiSD_Report.mega.outliers")
meg.out$stat <- rep("mega_outliers", nrow(meg.out))
meg.out$species <- rep("mega", nrow(meg.out))
colnames(meg.out) <- c("location","start","end","var","sfs","ld","mu","stat","species")
tor.out <- fread("raw_data/RAiSD_Report.toro.outliers")
tor.out$stat <- rep("toro_outliers", nrow(tor.out))
tor.out$species <- rep("toro", nrow(tor.out))
colnames(tor.out) <- c("location","start","end","var","sfs","ld","mu","stat","species")
meg.norm <- fread("raw_data/RAiSD_Report.mega.normal")
meg.norm$stat <- rep("mega_normal", nrow(meg.norm))
meg.norm$species <- rep("mega", nrow(meg.norm))
colnames(meg.norm) <- c("location","start","end","var","sfs","ld","mu","stat","species")
tor.norm <- fread("raw_data/RAiSD_Report.toro.normal")
tor.norm$stat <- rep("toro_normal", nrow(tor.norm))
tor.norm$species <- rep("toro", nrow(tor.norm))
colnames(tor.norm) <- c("location","start","end","var","sfs","ld","mu","stat","species")
sweep.df <- rbind.data.frame(meg.out,tor.out,meg.norm,tor.norm)
write.csv(sweep.df,"data/sweeps.csv",row.names = FALSE)

# subset top, compare distributions
mega.out <- subset(sweeps.df, stat=="mega_outliers", select="mu") 
mega.out <- as.vector(mega.out$mu)
mega.norm <- subset(sweeps.df, stat=="mega_normal", select="mu") 
mega.norm <- as.vector(mega.norm$mu)

# wilcox test between outliers and normal values
wilcox.test(mega.out, mega.norm) #W = 439380000, p-value < 2.2e-16

# load windowed pi values
sim.pi <- read.table("raw_data/sim.windowed.pi")[-1,]
meg.pi <- read.table("raw_data/meg.windowed.pi")[-1,]
compare <- nrow(sim.pi)
meg.pi <- meg.pi[sample(nrow(meg.pi), compare), ]
colnames(sim.pi) <- c("CHROM","BIN_START","BIN_END","N_VARIANTS","PI")
colnames(meg.pi) <- c("CHROM","BIN_START","BIN_END","N_VARIANTS","PI")
sim.pi$dataset <- rep("sim",nrow(sim.pi))
meg.pi$dataset <- rep("meg",nrow(meg.pi))
pi.df <- rbind.data.frame(sim.pi,meg.pi)
pi.df$PI <- as.numeric(as.character(pi.df$PI))
write.csv(pi.df,"data/pi.csv",row.names = FALSE)

# test difference
sim.pi.val <- as.numeric(as.character(sim.pi$PI))
meg.pi.val <- as.numeric(as.character(meg.pi$PI))
wilcox.test(sim.pi.val, meg.pi.val) # W = 771640, p-value < 2.2e-16

# redo windowed stats df for correlations
win.df <- fread("raw_data/syma_windows_overlap.csv")
win.df <- win.df[complete.cases(win.df),]
#df <- subset(df, sites > 10)
#restrict to good scaffolds
scaf <- as.data.frame(table(win.df$scaffold))
scaf <- subset(scaf, Freq > 10)
scaf <- as.vector(scaf$Var1)
win.df <- win.df[win.df$scaffold %in% scaf,]
win.df$windowID <- seq.int(nrow(win.df))
write.csv(win.df,"data/window_correlations.csv")

# get means and correlations (update ms)
lm1 <- lm(Fst_meg_toro ~ dxy_meg_toro, win.df)
summary(lm1) #p-value: < 2.2e-16, Adjusted R-squared:  0.1271
lm2 <- lm(Fst_meg_toro ~ pi_meg, win.df)
summary(lm2) #p-value: < 2.2e-16, 0.3785 
lm3 <- lm(Fst_meg_toro ~ pi_toro, win.df)
summary(lm3) #p-value: < 2.2e-16, 0.2427
lm4 <- lm(dxy_meg_toro ~ pi_toro, win.df)
summary(lm4) #p-value: < 2.2e-16, 0.9214  
lm5 <- lm(dxy_meg_toro ~ pi_meg, win.df)
summary(lm5) #p-value: < 2.2e-16, 0.8495 

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

#ID candidate sweep regions
subset(sum, refName == "chr11")
subset(sum, refName == "chr5")

