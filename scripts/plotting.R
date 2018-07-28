### plots for linck et al. 2018

library(ggplot2);library(viridis);library(maps);library(mapdata);
library(raster);library(rgdal);library(ggtree);library(cowplot)
sample_data <- read.csv("data/syma_spp_pcs.csv")
attach(sample_data)
map <- map_data("world2Hires")
toro_range <- shapefile("data/Syma_torotoro_22683566.shp")
mega_range <- shapefile("data/Syma_megarhyncha_22683569.shp")
toro_range <- fortify(toro_range)
mega_range <- fortify(mega_range)

a <- ggplot()+coord_map()+theme_bw()+ylim(-15,2)+xlim(130,154)+
  scale_color_viridis()+
  scale_shape_discrete(solid = F)+
  #scale_fill_manual(name="species",values = c("grey45","grey85"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=toro_range, aes(x=long,y=lat,group=group),fill="gray85")+
  geom_polygon(data=mega_range, aes(x=long,y=lat,group=group),fill="gray50")+
  #geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=sample_data,aes(x=long,y=lat,size=4,fill=PC1),pch=21,alpha=0.85,stroke=1)

# plot phylogenies

spp_tree <- read.nexus("data/syma_spp.svd.tre")
best_tree <- read.nexus("data/syma_best.tre")

#quartz()
b <- ggtree(best_tree, layout="daylight", branch.length = 'none',  inherit.aes = FALSE) + 
  geom_tiplab(color='black') 
  #geom_text2(aes(label=node), hjust=-.3) +
  #geom_cladelabel2(node = 1, label = "megarhyncha")
#dev.off()

# visualize uncertainty from bootstrapped trees
ggtree(spp_tree, layout="circular", color="lightblue", alpha=.3)
df <- fortify(best_tree, branch.length='none')
p+geom_tree(data=df, color='firebrick', layout="circular")

plot_grid(a, b, labels=c("A","B"))
