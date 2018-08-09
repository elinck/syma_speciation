### plots for linck et al. 2018

library(ggplot2);library(viridis);library(maps);library(mapdata);
library(raster);library(rgdal);library(ggtree);library(cowplot);library(ape)
library(png);library(gtable)
library(grid)

sample_data <- read.csv("data/syma_spp_pcs.csv")
#attach(sample_data)

# messy sample size assignment
sample_data$temp <- as.factor(paste(sample_data$lat, sample_data$long))
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

a <- ggplot()+coord_map()+theme_bw()+ylim(-14.5,2)+xlim(129,154)+
  scale_fill_viridis()+
  scale_shape_discrete(solid = F)+
  #scale_fill_manual(name="species",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=toro_range, aes(x=long,y=lat,group=group),fill="gray85")+
  geom_polygon(data=mega_range, aes(x=long,y=lat,group=group),fill="gray50")+
  #geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=sample_data,aes(x=long,y=lat,fill=PC1),pch=21,size=sample_data$size*2.5, alpha=0.85,stroke=1) +
  #scale_size(guide = "none") +
  annotate("segment", x = 147, xend = 149, y = -5.6, yend = -1) +
  annotate("text", x = 149, y = -0.5, label = "sellamontis") +
  annotate("segment", x = 150.5, xend = 153, y = -9.3, yend = -8) +
  annotate("text", x = 153, y = -7.5, label = "ochracea") +
  annotate("segment", x = 131.5, xend = 140, y = 0, yend = 1) +
  annotate("segment", x = 144, xend = 140, y = -3.3, yend = 1) +
  annotate("text", x = 140, y = 1.6, label = "torotoro") +
  annotate("segment", x = 147, xend = 148.5, y = -9.5, yend = -12) +
  annotate("text", x = 149, y = -12.3, label = "meeki") +
  annotate("segment", x = 142, xend = 140, y = -9.4, yend = -10) +
  annotate("segment", x = 143, xend = 140, y = -12.6, yend = -10) +
  annotate("text", x = 137.8, y = -10, label = "flavirostris") +
  annotate("segment", x = 144.5, xend = 144.5, y = -7.75, yend = -9.4) +
  annotate("text", x = 144.5, y = -9.75, label = "pseuestes") +
  annotate("segment", x = 134.2, xend = 134.2, y = -7, yend = -8.75) +
  annotate("text", x = 134.2, y = -9, label = "tentelare") +
  annotate("segment", x = 134.5, xend = 132.5, y = -4.1, yend = -4.5) +
  annotate("text", x = 131, y = -4.5, label = "wellsi") +
  annotate("segment", x = 136.5, xend = 136.3, y = -4.8, yend = -6.75) + 
  annotate("text", x = 136.5, y = -7.0, label = "pseuestes") +
  annotate("segment", x = 145, xend = 145, y = -5.5, yend = -0) +
  annotate("text", x = 145, y = 0.5, label = "megarhyncha")
  
# plot phylogenies
spp_tree <- read.nexus("data/syma_spp.svd.tre")
best_tree <- read.nexus("data/syma_best.tre")
best_tree <- ape::root(best_tree, node = 12, resolve.root = TRUE) # force root for spp. relationships

b <- ggtree(best_tree) + 
  geom_tiplab(color='black') +
  xlim(NA, 10) 

c <- readPNG("~/Dropbox/syma_speciation/data/syma.png")
c <- rasterGrob(c)


d <- ggplot(sample_data, aes(elevation, PC1)) +
  scale_fill_viridis()+
  geom_point(aes(fill=PC1),pch=21,size=4,alpha=0.85,stroke=1) +
  geom_smooth(method = "lm",se=F,linetype="dashed",col="black")

e <- ggplot(sample_data, aes(lat, PC2)) +
  scale_fill_viridis()+
  geom_point(aes(fill=PC1),pch=21,size=4,alpha=0.85,stroke=1) + 
  geom_smooth(method = "lm",se=F,linetype="dashed",col="black")


# gridded w/ cowplot
plot_grid(a, b, c, d, e, labels=c("A","B","C","D","E"), nrow = 2, ncol = 3, rel_widths = c(1,1,1,1,1))

# violin plot - reminder for later
m1 <- ggplot(data, aes(three_species_model, PC1, fill = three_species_model, guide=F)) +
  theme_light() +
  geom_point(size=.2,col="grey",position="jitter")+
  geom_violin(draw_quantiles = c(0.5)) +
  guides(fill=guide_legend(title="Species"))+
  scale_fill_manual(labels=c("megarhyncha","ochracea","torotoro"), values=symashade) + 
  theme(strip.background = element_rect(colour="black", fill="grey100")) +
  theme(axis.title.x=element_blank())

#histogram
m2 <- ggplot(data, aes(x=PC1, fill = three_species_model)) + 
  geom_histogram(aes(y=..density..), alpha=0.2, 
                 position="dodge",color="black", bins = 30)+
  theme_light() +
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  scale_fill_manual(labels=c("megarhyncha","ochracea","torotoro"), values=symashade) + 
  geom_density(alpha=0.9)+
  xlab("PC1") +
  ylab("Density")

pdf("~/Dropbox/remark_presentations/figures/morphology.pdf", width = 5, height = 7)
grid.arrange(m1,m2,ncol=1)
dev.off()

