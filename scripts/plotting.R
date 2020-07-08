# plots for linck et al. 2019

library(ggtree)
library(treeio)
library(tidyr)
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
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"))+
  scale_shape_manual(values=c(22, 24, 21)) +
  #scale_fill_manual(name="species",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=toro_range, aes(x=long,y=lat,group=group),fill="gray75")+
  geom_polygon(data=mega_range, aes(x=long,y=lat,group=group),fill="gray40")+
  #geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=sample_data,aes(x=long,y=lat,fill=sp,shape=sp,size=sample_data$size), alpha=0.85,stroke=1) +
  scale_size_continuous(range = c(3,8), breaks = c(1,2,3),name = "Samples") +
  labs(fill="Species",shape="Species") +
  guides(fill = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15))
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
  geom_point(aes(shape=sp)) +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"),guide=FALSE)+
  scale_color_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE)+
  scale_shape_manual(values=c(22, 24, 21),guide=FALSE) +
  theme_classic() +
  xlab("Genotype PC1") +
  ylab("Genotype PC2") +
  geom_polygon(data = hull_tor, alpha = 0.5) +
  geom_polygon(data = hull_meg, alpha = 0.5) +
  geom_polygon(data = hull_ochr, alpha = 0.5) +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))

# subset to show clustering results unaffected by hyrad data
sample_data_wgs <- read.csv("data/syma_spp_pcs_wgs.csv")
sample_data_wgs$sp <- as.character(sample_data_wgs$sp)
sample_data_wgs$sp[which(sample_data_wgs$ssp == "ochracea")] <- "ochracea"
sample_data_wgs$sp <- as.factor(sample_data_wgs$sp)

# convex hull for PCA, subset data
hull_tor <- sample_data_wgs[sample_data_wgs$sp=="torotoro",] %>%
  slice(chull(PC1, PC2))
hull_meg <- sample_data_wgs[sample_data_wgs$sp=="megarhyncha",] %>%
  slice(chull(PC1, PC2))
hull_ochr <- sample_data_wgs[sample_data_wgs$sp=="ochracea",] %>%
  slice(chull(PC1, PC2))

# version for supplement, wgs data only: 
b.2 <- ggplot(data=sample_data_wgs,
              aes(x=PC1,y=PC2,fill=sp,color=sp)) + 
  geom_point(aes(shape=sp)) +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"),guide=FALSE)+
  scale_color_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE)+
  scale_shape_manual(values=c(22, 24, 21), guide=FALSE) +
  theme_classic() +
  xlab("Genotype PC1") +
  ylab("Genotype PC2") +
  geom_polygon(data = hull_tor, alpha = 0.5) +
  geom_polygon(data = hull_meg, alpha = 0.5) +
  geom_polygon(data = hull_ochr, alpha = 0.5) +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))

# subset to show clustering results unaffected by excluding modern data
sample_data_old <- read.csv("data/syma_spp_pcs_old.csv")
sample_data_old$sp <- as.character(sample_data_old$sp)
sample_data_old$sp[which(sample_data_old$ssp == "ochracea")] <- "ochracea"
sample_data_old$sp <- as.factor(sample_data_old$sp)

# convex hull for PCA, subset data
hull_tor <- sample_data_old[sample_data_old$sp=="torotoro",] %>%
  slice(chull(PC1, PC2))
hull_meg <- sample_data_old[sample_data_old$sp=="megarhyncha",] %>%
  slice(chull(PC1, PC2))
hull_ochr <- sample_data_old[sample_data_old$sp=="ochracea",] %>%
  slice(chull(PC1, PC2))

# version for supplement, wgs data only: 
b.3 <- ggplot(data=sample_data_old,
              aes(x=PC1,y=PC2,fill=sp,color=sp)) + 
  geom_point(aes(shape=sp)) +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"),guide=FALSE)+
  scale_color_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE)+
  scale_shape_manual(values=c(22, 24, 21), guide=FALSE) +
  theme_classic() +
  xlab("Genotype PC1") +
  ylab("Genotype PC2") +
  geom_polygon(data = hull_tor, alpha = 0.5) +
  geom_polygon(data = hull_meg, alpha = 0.5) +
  geom_polygon(data = hull_ochr, alpha = 0.5) +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))

# plot vocalization data
vdata <- read.csv("data/syma_spp_calls.csv")

# PC1 only
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

hull_tor <- vdata[vdata$two_species=="torotoro",] %>%
  slice(chull(PC1, PC2))
hull_meg <- vdata[vdata$two_species=="megarhyncha",] %>%
  slice(chull(PC1, PC2))
hull_ochr <- vdata[vdata$two_species=="ochracea",] %>%
  slice(chull(PC1, PC2))

c.2 <- ggplot(data=vdata,aes(x=PC1,y=PC2,fill=two_species,color=two_species)) + 
  geom_point(aes(shape=two_species)) +
  scale_fill_manual(values = c("#73D056FF","#443B84FF"),guide=FALSE)+
  scale_color_manual(values = c("#73D056FF","#443B84FF"), guide=FALSE)+
  scale_shape_manual(values=c(22, 24, 21),guide=FALSE) +
  theme_classic() +
  xlab("Acoustic PC1") +
  ylab("Acoustic PC2") +
  geom_polygon(data = hull_tor, alpha = 0.5) +
  geom_polygon(data = hull_meg, alpha = 0.5) +
  geom_polygon(data = hull_ochr, alpha = 0.5) +
  labs(fill="Species") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))  

# plot morphology
mdata <- read.csv("data/syma_spp_morphology_log.csv")

# PC1 only
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

# convex hull for full PCA plot
hull_tor <- mdata[mdata$three_species_model=="torotoro",] %>%
  slice(chull(PC1, PC2))
hull_meg <- mdata[mdata$three_species_model=="megarhyncha",] %>%
  slice(chull(PC1, PC2))
hull_ochr <- mdata[mdata$three_species_model=="ochracea",] %>%
  slice(chull(PC1, PC2))

d.2 <- ggplot(data=mdata,aes(x=PC1,y=PC2,fill=three_species_model,color=three_species_model)) + 
  geom_point(aes(shape=three_species_model)) +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"),guide=FALSE)+
  scale_color_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE)+
  scale_shape_manual(values=c(22, 24, 21),guide=FALSE) +
  theme_classic() +
  xlab("Morphological PC1") +
  ylab("Morphological PC2") +
  geom_polygon(data = hull_tor, alpha = 0.5) +
  geom_polygon(data = hull_meg, alpha = 0.5) +
  geom_polygon(data = hull_ochr, alpha = 0.5) +
  labs(fill="Species") +
  theme(legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.title=element_text(size=15))  

# faceted traits for supplement
mdata_sub <- cbind.data.frame(mdata$three_species_model,mdata$sex,mdata$bill_from_nostril,mdata$bill_width,mdata$bill_depth,
                              mdata$tarsus,mdata$wing_chord,mdata$tail_length)
colnames(mdata_sub) <- c("three_species_model", "sex", "Bill from nostril", "Bill width", "Bill depth", "Tarsus",
                     "Wing chord", "Tail length")
mdata_tidy <- gather(mdata_sub, `Bill from nostril`, `Bill width`, `Bill depth`, `Tarsus`, `Wing chord`, `Tail length`, key="trait",value="Millimeters")

d.3 <- ggplot(mdata_tidy) +
  theme_bw() +
  #geom_density(aes(x=Millimeters,fill=three_species_model), alpha=0.5) +
  geom_histogram(aes(x=Millimeters,fill=three_species_model), bins=22) +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF")) +
  facet_wrap(~trait, scales="free") +
  labs(fill="Species")  +
  ylab("Density") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

d.4 <- ggplot(mdata_tidy[mdata_tidy$three_species_model=="torotoro",]) +
  theme_bw() +
  #geom_density(aes(x=Millimeters,fill=sex), alpha=0.5) +
  geom_histogram(aes(x=Millimeters,fill=sex), bins=22) +
  #scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE) +
  facet_wrap(~trait, scales="free") +
  labs(fill="Sex")  +
  ylab("Density") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

d.5 <- ggplot(mdata_tidy[mdata_tidy$three_species_model=="megarhyncha",]) +
  theme_bw() +
  #geom_density(aes(x=Millimeters,fill=sex), alpha=0.5) +
  geom_histogram(aes(x=Millimeters,fill=sex), bins=22) +
  #scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE) +
  facet_wrap(~trait, scales="free") +
  labs(fill="Sex")  +
  ylab("Density") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

d.6 <- ggplot(mdata_tidy[mdata_tidy$three_species_model=="ochracea",]) +
  theme_bw() +
  #geom_density(aes(x=Millimeters,fill=sex), alpha=0.5) +
  geom_histogram(aes(x=Millimeters,fill=sex), bins=22) +
  #scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"), guide=FALSE) +
  facet_wrap(~trait, scales="free") +
  labs(fill="Sex")  +
  ylab("Density") +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

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

# load kingfisher picture
h <- readPNG("~/Dropbox/syma_speciation/ignore/figures/syma_2x1.png")
h <- rasterGrob(h)

# generate panels for figure 1
top <- plot_grid(a,e,h,nrow=1,labels=c("A","B","C"),rel_widths = c(1.9,1,1))
bottom <- plot_grid(b,c.2,d.2,nrow=1,labels=c("D","E","F"),rel_widths = c(1,1,1))
fig1 <- plot_grid(top,bottom,labels=NULL,nrow=2,rel_heights = c(1.25,1))

# write to file
png("manuscript/syma_1_rev.png",res=300,width=14.5,height=7.5,units="in")
fig1
dev.off()

# write supplement pca plots to file
png("manuscript/syma_S1.png",res=300,width=8,height=4,units="in")
plot_grid(b.2,b.3,labels="AUTO",nrow=1)
dev.off()

# write supplement trait info to file
png("manuscript/syma_S2.png",res=300,width=12,height=9,units="in")
plot_grid(d.3,d.4,d.5,d.6,labels="AUTO",nrow=2)
dev.off()
