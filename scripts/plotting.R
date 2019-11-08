# plots for linck et al. 2019

library(ggtree)
library(treeio)
library(ape)
library(cowplot)
library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(raster)
library(rgdal)
library(ggrepel)
library(grid)
library(dplyr)
library(plyr)
library(png)
library(data.table)
library(gghighlight)

setwd("~/Dropbox/syma_speciation/")

###  Figure 1

# prepare PC and map data
sample_data <- read.csv("data/syma_spp_pcs.csv")
attach(sample_data)
sample_data$sp <- as.character(sample_data$sp)
sample_data$sp[which(sample_data$ssp == "ochracea")] <- "ochracea"
sample_data$sp <- as.factor(sample_data$sp)
sp_means <- aggregate(PC1 ~ sp, data = sample_data, mean)
meg_val <- sp_means$PC1[1]
och_val <- sp_means$PC1[2]
tor_val <- sp_means$PC1[3]
sample_data$col <- ifelse(sample_data$sp == "megarhyncha", 
                          meg_val, ifelse(sample_data$sp == "torotoro",
                                          tor_val,och_val))
sample_data$temp <- as.factor(paste(sample_data$lat, sample_data$long)) # messy sample size assignment
table(sample_data$temp) #see which points have multiple samples
sample_data$size <- 1
sample_data$size[sample_data$temp == "-5.95 134.57"] <- 3
sample_data$size[sample_data$temp == "-2.72 134.5"] <- 2
sample_data$size[sample_data$temp == "-3.683 143.633"] <- 2
sample_data$size[sample_data$temp == "-6.4583 147.4333"] <- 2
sample_data$size[sample_data$temp == "-8.25 142.98"] <- 2
sample_data$size[sample_data$temp == "9.5 150.666667"] <- 2

map <- map_data("world")
toro_range <- shapefile("data/Syma_torotoro_22683566.shp")
mega_range <- shapefile("data/Syma_megarhyncha_22683569.shp")
toro_range <- fortify(toro_range)
mega_range <- fortify(mega_range)

# plot map
a <- ggplot()+coord_map()+theme_classic()+ylim(-14.5,2)+xlim(129,154)+
  scale_fill_viridis()+
  scale_shape_discrete(solid = F)+
  #scale_fill_manual(name="species",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=toro_range, aes(x=long,y=lat,group=group),fill="gray75")+
  geom_polygon(data=mega_range, aes(x=long,y=lat,group=group),fill="gray40")+
  #geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=sample_data,aes(x=long,y=lat,fill=PC1,size=sample_data$size), pch=21, alpha=0.85,stroke=1) +
  scale_size_continuous(range = c(3,8), breaks = c(1,2,3),name = "Samples") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) #+
#theme(legend.justification=c(1,1), legend.position=c(1,1))

# convex hull for PCA
hull_tor <- sample_data[sample_data$sp=="torotoro",] %>%
  slice(chull(PC1, PC2))
hull_meg <- sample_data[sample_data$sp=="megarhyncha",] %>%
  slice(chull(PC1, PC2))
hull_ochr <- sample_data[sample_data$sp=="ochracea",] %>%
  slice(chull(PC1, PC2))

# pca
b <- ggplot(data=sample_data,aes(x=PC1,y=PC2,fill=sp,color=sp)) + 
  geom_point(pch=1) +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"))+
  scale_color_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE)+
  theme_classic() +
  xlab("Genotype PC1") +
  ylab("Genotype PC2") +
  geom_polygon(data = hull_tor, alpha = 0.5) +
  geom_polygon(data = hull_meg, alpha = 0.5) +
  geom_polygon(data = hull_ochr, alpha = 0.5) +
  labs(fill="Species") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))

vdata <- read.csv("data/syma_spp_calls.csv")
c <- ggplot(vdata, aes(x=PC1, fill = two_species)) + 
  theme_classic() +
  scale_fill_manual(values = c("#73D056FF","#443B84FF"))+
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  geom_histogram(bins=20, alpha=0.9,color="black")+
  theme(legend.position = "none") +
  xlab("Acoustic PC1") +
  ylab("Density") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))

# plot morphology
mdata <- read.csv("data/syma_spp_morphology.csv")
d <- ggplot(mdata, aes(x=PC1, fill = three_species_model)) + 
  theme_classic() +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"))+
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  geom_histogram(bins=20, alpha=0.9,color="black")+
  theme(legend.position = "none") +
  xlab("Morphological PC1") +
  ylab("Density") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))

top <- plot_grid(a,b,c,d,labels="AUTO",nrow=1,rel_widths=c(1.5,0.75,0.5,0.5))

# read mcc ND2 tree from beast
beast <- read.beast("~/Dropbox/syma_speciation/data/syma_ND2.tree")
beast@data[isTip(beast, beast@data$node),]$length_0.95_HPD <- NA

# fortify data, change tip labels
beast <- fortify(beast)
beast$label <- c("EL10_toro","EL11_toro","EL13_toro","EL18_mega","EL19_mega","EL1_mega","EL20_mega",  
                 "EL21_toro","EL23_mega","EL24_mega","EL27_mega","EL29_ochr","EL32_toro","EL39_toro",  
                 "EL40_mega","EL41_toro","EL42_toro","EL43_toro","EL44_toro","EL45_toro","EL46_toro",  
                 "EL47_toro","EL48_ochr","EL4_mega" ,"EL50_toro","EL51_toro","EL52_toro",  
                 "EL53_toro","EL54_toro","EL55_toro","EL57_toro","EL59_toro",  
                 "EL5_ochr" ,"EL60_toro","EL6_mega" ,"EL8_toro" ,"EL9_toro" ,"t_sanctus", rep("NA",37))
# set up base tree, plot nodes
tmp1 <- ggtree(beast, mrsd="2000-12-12") + 
  theme_tree2() +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

# drop node labels
tmp1 <- ggtree(beast, mrsd="2000-12-12") + 
  theme_tree2()

# add colors
clr <- beast$label[1:38]
clr[grepl("toro",clr)]<-"#443B84FF" 
clr[grepl("mega",clr)]<-"#73D056FF"
clr[grepl("ochr",clr)]<-"#20938CFF"
clr[38] <- "black"

# plot final tree
e <- revts(tmp1) +
  #geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, color=clr) + 
  geom_range("height_0.95_HPD", color='black', size=2, alpha=.5, branch.length="height")  + 
  xlim(-2500000,400000) +
  theme(axis.text.x = element_text(size=15)) +
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)>.99, 
                 x=branch), hjust=1.25,vjust=-0.5, size=5) +
  geom_strip('2', '34', barsize=20, color="#443B84FF", 
             label="",alpha=0.5) +
  geom_strip('10', '5', barsize=20, color="#73D056FF", 
             label="",alpha=0.5) +
  geom_strip('23', '33', barsize=20, color="#20938CFF", 
             label="",alpha=0.5) 

#CI for initial divergence date
e$data$height_0.95_HPD[[40]]

#estimate of initial divergence date
(607421.8+972618.8)/2 #790020.3

#CI for meg / ochr split
e$data$height_0.95_HPD[[64]]

#estimate of meg / ochr divergence date
(243261.9+512732.0)/2 #377997

#CI for intra torotoro divergence
e$data$height_0.95_HPD[[41]]

# estimate of intra torotoro divergence
(304804.8+552355.0)/2 #428579.9

# SVDquartets trees
svdq <- ape::read.nexus("~/Dropbox/syma_speciation/data/syma.svd.tre")
consens <- ape::consensus(svdq, p=0.5)
consens <- root(consens, 1)
bs <- prop.clades(consens, svdq, rooted = TRUE)
consens <- fortify(consens)
consens$bs <- c(rep(NA, 20), bs)
consens$label[15] <- "EL40_mega"

clr2 <- consens$label[1:20]
clr2[grepl("toro",clr2)]<-"#443B84FF" 
clr2[grepl("mega",clr2)]<-"#73D056FF"
clr2[grepl("ochr",clr2)]<-"#20938CFF"
clr2[1] <- "black"

# root tree
f <- ggtree(consens) + 
  theme_tree() +
  #geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, color=clr2) + 
  xlim(-1,10) +
  geom_nodelab(aes(label=bs,size=5),vjust=-0.4,hjust=1.1)+
  geom_tippoint(color=clr2,alpha=0.5,size=5,shape=15) #+

# plot svdq species trees
sptree <- ape::read.nexus("~/Dropbox/syma_speciation/data/syma_sptree_concordant.tre")
consens <- ape::consensus(sptree, p=0.5)
consens <- root(consens, 1)
bs <- prop.clades(consens, sptree, rooted = TRUE)
consens <- fortify(consens)
consens$bs <- c(rep(NA, 4), bs)
clr3 <- consens$label[1:4]
clr3[grepl("toro",clr3)]<-"#443B84FF" 
clr3[grepl("mega",clr3)]<-"#73D056FF"
clr3[grepl("ochr",clr3)]<-"#20938CFF"
clr3[1] <- "black"

g <- ggtree(consens) +
  theme_tree() +
  xlim(-1,5) +
  geom_tippoint(color=clr3,alpha=0.5,size=15,shape=15) +
  geom_nodelab(aes(label = bs),size=5,vjust=-0.4,hjust=1.1)

h <- readPNG("~/Dropbox/syma_speciation/figures/syma_2x1.png")
h <- rasterGrob(h)

top <- plot_grid(a,h,nrow=1,labels=c("A","B"),rel_widths = c(1.25,1))
middle <- plot_grid(b,c,d, nrow=1,labels=c("C","D","E"), rel_widths = c(1.75,1,1))
bottom <- plot_grid(e,f,g,nrow=1,labels=c("F","G","H"),rel_widths = c(1.5,1,1))
fig1 <- plot_grid(top,middle,bottom,labels=NULL,nrow=3,rel_heights = c(1,0.8,1))

# write to file
pdf("figures/syma_1.pdf", width=14,height=12)
fig1
dev.off()

# figure 2
dstat.df <- read.csv("data/dstats.csv")
ititle <- expression(paste(italic("D"), "-tests"))
i <-  ggplot(dstat.df, aes(x=D,fill=topology)) + 
  theme_classic() +
  guides(fill=guide_legend(title="Topology"))+
  scale_fill_manual(values=c("grey90","grey60","grey30")) +
  xlim(-0.1,0.5) +
  #scale_x_continuous(breaks=c(-0.1,0.0,0.1,0.2,0.3,0.4,0.5)) +
  geom_vline(xintercept=0.0,linetype="dashed",color="red",size=1) +
  geom_histogram(bins=90, alpha=0.6,color="black")+
  labs(title=ititle) +
  ylab("Count") +
  xlab(expression(italic("D")))+
  #facet_wrap(~topology, scales="free") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15))

mig_real.df <- read.csv("data/real_migration_rate.csv")
tmp2 <- ggplot(mig_real.df, aes(x=value, fill = parameter)) + 
  theme_classic() +
  guides(fill=guide_legend(title="Parameter"))+
  scale_fill_manual(values=c("#32a852","#8c32a8")) +
  geom_histogram(bins=50, alpha=0.6,color="black")+
  ylab("Count") +
  xlab("Migration rate") +
  labs(title="Migration rate by time period") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15))

time.df <- read.csv("data/divergence_times.csv")
tmp3 <- ggplot(time.df, aes(x=value, fill = parameter)) + 
  theme_classic() +
  guides(fill=guide_legend(title="Parameter"))+
  scale_fill_manual(values=c("#32a852","#8c32a8","#ffbf00")) +
  geom_histogram(bins=50, alpha=0.6,color="black")+
  ylab("Count") +
  xlab("Years") +
  labs(title="Divergence time") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15))

meg.df <- read.csv("data/meg_sizes.csv")
tmp4title <- expression(paste(italic("S. megarhyncha"), " population size"))
tmp4 <- ggplot(meg.df, aes(x=value, fill = parameter)) + 
  theme_classic() +
  guides(fill=guide_legend(title="Parameter"))+
  scale_fill_manual(values=c("#32a852","#8c32a8","#ffbf00")) +
  geom_histogram(bins=50, alpha=0.6,color="black")+
  ylab("Count") +
  xlab("Population") +
  labs(title = tmp4title) +
  scale_x_continuous(breaks=c(3e+05,5e+05,7e+05,9e+05)) +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15))

meg.df <- read.csv("data/tor_sizes.csv")
tmp5title <- expression(paste(italic("S. torotoro"), " population size"))
tmp5 <- ggplot(tor.df, aes(x=value, fill = parameter)) + 
  theme_classic() +
  guides(fill=guide_legend(title="Parameter"))+
  scale_fill_manual(values=c("#32a852","#8c32a8","#ffbf00")) +
  geom_histogram(bins=50, alpha=0.6,color="black")+
  ylab("Count") +
  xlab("Population") +
  labs(title = tmp5title) +
  scale_x_continuous(breaks=c(3e+05,5e+05,7e+05,9e+05)) +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15))

j  <- plot_grid(tmp2,tmp3,tmp4,tmp5, nrow=2)
fig2 <- plot_grid(i,j,labels=c("A","B"),nrow=2,rel_heights = c(0.3,0.7))
pdf("figures/syma_2.pdf", width=9,height=10)
fig2
dev.off()

### figure 3
win.df <- read.csv("data/window_stats_chr.csv")[,-1]
n <- ggplot(data=win.df,aes(x=row,y=value,col=chr))+
  theme_bw()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=15),
        strip.text = element_text(size=15),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank())+  
  scale_y_continuous(breaks=seq(0,0.6,0.2),minor_breaks = NULL)+
  scale_color_manual(values=rep(c("grey45","grey70"),length(unique(win.df$chr))/2+1))+
  geom_point(data=subset(win.df,value>0.0),size=0.75,shape=21,alpha=0.3)+
  geom_point(data=subset(win.df,chr==5),size=0.75,shape=21,alpha=0.7)+
  geom_point(data=subset(win.df, outlier=="TRUE"),size=1,shape=21,col="red")+
  geom_hline(data=subset(win.df, stat=="F[ST]"),aes(yintercept=0.4084545),linetype="dashed",lwd=0.25)+
  geom_hline(data=subset(win.df, stat=="D[XY]"),aes(yintercept=0.17),linetype="dashed",lwd=0.25)+
  geom_line(aes(y=rollmean),lwd=0.15,col="black")+
  geom_segment(data=chr_labels,aes(x=start+250,xend=stop-250,y=-0.025,yend=-0.025,col=NA),col="black")+
  #geom_vline(xintercept=c(103017,114851)) +
  geom_text(data=chr_labels,aes(label=chr,x=mid,y=-0.1,col=NA),
            col="black",size=5,angle=0,direction="y",box.padding = 0.15,
            segment.size=0.2) +
  theme(strip.background = element_rect(colour="black", fill="grey100")) +
  facet_grid(stat ~ ., labeller=label_parsed,scales = "free_y") 
annotation_custom(ggplotGrob(q),
                  xmin=10, xmax=10000, ymin=0, ymax=1)


win.chr.fst.df <- read.csv("data/chromsomes_fst.csv")[,-1]
r <- ggplot(data=win.chr.fst.df,aes(x=row,y=value))+theme_bw()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1))+  
  ylab(ylab.text) +
  xlab("Position") +
  scale_y_continuous(breaks=seq(0,1.0,0.2),minor_breaks = NULL)+
  geom_line(aes(y=fst_rollmean),lwd=0.5,col="black")+
  #geom_hline(aes(yintercept=0.409031),linetype="dashed",lwd=0.25)+
  facet_grid(. ~ chr_order, scales="free_x") +
  theme(strip.background = element_rect(colour="black", fill="grey100")) +
  theme(legend.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.title=element_text(size=15))

win.chr.dxy.df <- read.csv("data/chromsomes_dxy.csv")[,-1]
x <- ggplot(data=win.chr.dxy.df,aes(x=row,y=value))+theme_bw()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1))+  
  ylab(ylab.text) +
  xlab("Position") +
  scale_y_continuous(breaks=c(0.1),minor_breaks = NULL)+
  geom_line(aes(y=dxy_rollmean),lwd=0.5,col="black")+
  #geom_hline(aes(yintercept=quantile(df2$dxy_meg_toro,0.995, na.rm = TRUE)),linetype="dashed",lwd=0.25)+
  facet_grid(. ~ chr_order, scales="free_x") +
  theme(strip.background = element_rect(colour="black", fill="grey100")) +
  theme(legend.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.title=element_text(size=15))

bottom <- plot_grid(r, x, nrow=2)

outlier.df <- read.csv("data/outlier_correlations.csv")[,-1]
s <- ggplot(data=outlier.df,aes(x=dxy_meg_toro,y=Fst_meg_toro)) +
  theme_classic()+
  #scale_fill_distiller(palette = "BrBG",direction = -1,name="log(n. windows)",guide=FALSE)+
  geom_point(data=subset(outlier.df,fst_outlier=="FALSE" & dxy_outlier=="FALSE" & both_outlier=="FALSE"),shape=21,aes(col="grey80"))+
  geom_point(data=subset(outlier.df,fst_outlier_05=="TRUE"),shape=21,aes(col="grey60"))+
  geom_point(data=subset(outlier.df,dxy_outlier_05=="TRUE"),shape=21,aes(col="grey60"))+
  geom_point(data=subset(outlier.df,fst_outlier=="TRUE"),shape=21,aes(col="grey40"))+
  geom_point(data=subset(outlier.df,dxy_outlier=="TRUE"),shape=21,aes(col="grey40"))+
  #geom_point(data=subset(df2,both_outlier=="TRUE"),shape=21,aes(col="red")) +
  geom_hline(aes(yintercept=quantile(outlier.df$Fst_meg_toro,0.995, na.rm = TRUE),linetype="dashed"),lwd=0.25)+
  geom_hline(aes(yintercept=quantile(outlier.df$Fst_meg_toro,0.95, na.rm = TRUE),linetype="solid"),lwd=0.25)+
  geom_vline(aes(xintercept=quantile(outlier.df$dxy_meg_toro,0.995, na.rm = TRUE),linetype="dashed"),lwd=0.25)+
  geom_vline(aes(xintercept=quantile(outlier.df$dxy_meg_toro,0.95, na.rm = TRUE),linetype="solid"),lwd=0.25)+
  xlab(expression(D[XY])) +
  ylab(expression(F[ST])) +
  #geom_smooth(method="lm",fill=NA,col="black",lwd=0.65) +
  scale_color_manual(name = element_blank(), 
                     values =c('grey40'='grey40','grey60'='grey60','grey80'='grey80','red'='red'), 
                     labels = c('outlier at 0.005','outlier at 0.05','nonoutlier','joint outlier at 0.05')) +
  scale_linetype_manual(name = element_blank(), 
                        values=c("dashed","solid"),
                        labels=c("outlier at 0.005","outlier at 0.05")) +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))


top <- plot_grid(n, labels=c("A"))
bottom2 <- plot_grid(bottom, s, labels=c("B","C"), ncol=2, rel_widths = c(0.75,0.5))
pdf("figures/syma_3.pdf",width=18,height=12)
plot_grid(top, bottom2,nrow=2, rel_heights=c(0.75,0.5))
dev.off()


