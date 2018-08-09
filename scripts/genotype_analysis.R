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

# random seeds for structure
sample(1000000:2000000,20) 
#  [1] 1717248 1012603 1652698 1531375 1211848 1417797 1126032 1086360
#  [9] 1051791 1892453 1973568 1954315 1491487 1437422 1052922 1761692
# [17] 1997388 1430525 1250546 1296037

### demographic inference

# conversion from popgenome format -- not working
# VCF_split_into_scaffolds("raw_data/syma.85p.1perlocus.vcf","vcf_genome")
# syma.genome <- readData("raw_data/vcf_genome", format = "VCF")
# syma.segsites <- as.segsites(syma.genome@)
# ss <- as.segsites.list(syma.genome)
# is_segsites(ss)

# each locus -> matrix entry in list, paired  with position data
# number of loci in model must match list
# question -- can you save likelihood values from each run?


# pcadapt format work around
pcadapt::vcf2pcadapt("raw_data/syma.85p.1perlocus.vcf", "syma85.snps")
snps <- t(as.matrix(read.table("syma85.snps")))
pos <- read.table("positions.txt")
pos <- as.vector(pos$V1)

n <- 1
col <- ncol(snps)
list.df <- list()
list.df <- split(snps, rep(1:ceiling(col/n), each=n, length.out = col))
list.pos <- split(pos, rep(1:ceiling(col/n), each=n, length.out = col))

segsites <- list()
for(i in 1:300){
  a <- matrix(list.df[[i]],nrow=40,ncol=1)
  b <- as.vector(list.pos[[i]])
  segsites[[i]] <- create_segsites(a, b, check = TRUE)
}

#models 
im <- coal_model(c(28, 12), loci_number = 300, loci_length = 1, ploidy = 2) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_migration(par_range("m", 0, 3), symmetric = TRUE) +
  feat_pop_merge(par_range("t_split", 0.1, 2), 2, 1) + 
  sumstat_jsfs()

hybrid <- coal_model(c(10, 10), loci_number = 300, loci_length = 1, ploidy = 2) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_migration(par_range("m", 0, 3), symmetric = TRUE, time = 2) +
  feat_pop_merge(par_range("t_split", 0.1, 2), 2, 1) + 
  sumstat_jsfs()

i <- coal_model(c(10, 10), loci_number = 300, loci_length = 1, ploidy = 2) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_pop_merge(par_range("t_split", 0.1, 2), 2, 1) + 
  sumstat_jsfs()

## sumstats 
sumstats_im <- calc_sumstats_from_data(im, segsites)
sumstats_hybrid <- calc_sumstats_from_data(hybrid, segsites)
sumstats_i <- calc_sumstats_from_data(i, segsites)

#hybrid model test
hybrid_sim <- create_jaatha_model(hybrid)
hybrid_empirical <- create_jaatha_data(sumstats_hybrid, hybrid_sim)
estimates_hybrid <- jaatha(hybrid_sim, hybrid_empirical, 
                           sim = 100, repetitions = 2, verbose = FALSE)
#im model test
im_sim <- create_jaatha_model(im)
im_empirical <- create_jaatha_data(sumstats_im, jaatha_im)
estimates_im <- jaatha(im_sim, im_empirical, 
                       sim = 100, repetitions = 2, verbose = FALSE)

#isolation model test
i_sim <- create_jaatha_model(i)
i_empirical <- create_jaatha_data(sumstats_i, i_sim)
estimates_i <- jaatha(i_sim, i_empirical, 
                      sim = 100, repetitions = 2, verbose = FALSE)

estimates_hybrid$loglikelihood #[1] -3000.998
estimates_im$loglikelihood #[1] -2136.74
estimates_i$loglikelihood #[1] -3008.813

#plot phylogenies
mtdna <- read.tree("syma_mtDNA_14k.tre")
nuc <- read.tree("syma_species_tree.tre")
nuc <- root(nuc, outgroup = "ochracea", resolve.root = TRUE)

ggtree(nuc) +
  geom_tiplab(size=5, color="purple")

pdf("mtdna.tree.pdf")
plot(mtdna, use.edge.length = FALSE)
dev.off()

pdf("species.tree.pdf")
plot(nuc)
dev.off()

bs <- read.nexus("data/syma_lin.svd.tre", force.multi = TRUE, tree.names = NULL)
tree <- bs$B_2.1
clad <- prop.clades(tree, bs, rooted = FALSE)
boot <- prop.clades(tree, bs)
boot <- (boot/365)*100
options(digits = 2)
layout(1)
par(mar = rep(2, 4))
plot(tree, main = "Bipartition vs. Clade Support Values")
drawSupportOnEdges(boot)
legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
       pt.bg = c("green", "lightblue"), pt.cex = 2.5)

plot(tree)

p <- ggtree(bs, layout = "rectangular", color = "lightblue")

tree <- plotBS(tree,bs$trees,type="phylogram")

