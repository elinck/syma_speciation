setwd("~/Dropbox/syma_wgs")
library(vcfR);library(adegenet);library(ggplot2);
library(ape);library(strataG);library(data.table);
library(pcadapt);library("qvalue");library("OutFLANK");
library("ggplot2");library(vcfR);library(PopGenome)
library(plyr);library(ape);library(ggtree)
library(magrittr)

# 85% complete matrix
syma.med <- read.vcfR("syma.gatk.85p.d5.maf05.recode.vcf")
dna <- vcfR2DNAbin(syma.med, unphased_as_NA = F, consensus = T, extract.haps = F)
syma <- DNAbin2genind(dna)
seq.scaled <- scaleGen(syma,NA.method="mean",scale=F)
#snpID <- as.data.frame(syma.med.genind$loc.n.all)

samples <- as.character(read.table("samples.txt")[[1]]) 
rownames(syma@tab) <- samples

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
pc$sampleID <- samples
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust$clust.k2))+geom_text(aes(label=sampleID))

# sliding window analyses
#VCF_split_into_scaffolds("syma.gatk.85p.d5.maf05.recode.vcf","syma_split_vcf") 
syma.popgenome <- readData("syma_split_vcf",format="VCF") # read in folder of scaffolds


pop1 <- get.individuals(syma.popgenome)[[1]] %>% grep("mega",.,value=T)
pop2 <- get.individuals(syma.popgenome)[[1]] %>% grep("toro",.,value=T)
pop3 <- get.individuals(syma.popgenome)[[1]] %>% grep("ochr",.,value=T)
syma.popgenome <- set.populations(syma.popgenome,list(pop1,pop2,pop3),diploid=TRUE) 
syma.popgenome <- F_ST.stats(syma.popgenome)
get.F_ST(syma.popgenome)
new <- sliding.window.transform(syma.popgenome,500,500,type=2,whole.data=FALSE)
new <- set.populations(syma.popgenome,list(pop1,pop2,pop3),diploid=TRUE) 
genome.pos <- sapply(sliding@region.names, function(x){
  split <- strsplit(x," ")[[1]][c(1,3)]
  val <- mean(as.numeric(split))
  return(val)
})       
length(new@region.names)
sliding <- diversity.stats(sliding)

win_snp <- F_ST.stats(sliding)

win_fst <- win_snp@nucleotide.F_ST[,1]
bb_div  <- win_snp@nuc.diversity.within[,1] # diversity among B (bb = "big B")
lb_div  <- win_snp@nuc.diversity.within[,2] # diversity among B (lb = "little b")


plot(1:length(win_fst), win_fst)

par(mfrow=c(2,1))
win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(bb_div), bb_div)

win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(lb_div), lb_div)

