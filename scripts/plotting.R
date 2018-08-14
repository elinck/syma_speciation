### plots for linck et al. 2018

library(ggplot2);library(viridis);library(maps);library(mapdata);
library(raster);library(rgdal);library(ggtree);library(cowplot);library(ape)
library(png);library(gtable)
library(grid);library(scales)

sample_data <- read.csv("data/syma_spp_pcs.csv")
#attach(sample_data)

sample_data$sp <- as.character(sample_data$sp)
sample_data$sp[which(sample_data$ssp == "ochracea")] <- "ochracea"
sample_data$sp <- as.factor(sample_data$sp)

#get mean PC values
sp_means <- aggregate(PC1 ~ sp, data = sample_data, mean)
meg_val <- sp_means$PC1[1]
och_val <- sp_means$PC1[2]
tor_val <- sp_means$PC1[3]
sample_data$col <- ifelse(sample_data$sp == "megarhyncha", 
                          meg_val, ifelse(sample_data$sp == "torotoro",
                                                         tor_val,och_val))

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
  geom_point(aes(fill=PC1),pch=21,size=3,alpha=0.85,stroke=1) +
  geom_smooth(method = "lm",se=F,linetype="dashed",col="black") +
  xlab("Elevation") +
  theme(legend.position = "none")

e <- ggplot(sample_data, aes(lat, PC2)) +
  scale_fill_viridis()+
  geom_point(aes(fill=PC1),pch=21,size=3,alpha=0.85,stroke=1) + 
  geom_smooth(method = "lm",se=F,linetype="dashed",col="black") +
  xlab("Latitude") +
  theme(legend.position = "none")

# plot AIC
aic_df <- read.csv("data/aic_plotting.csv")
f <- ggplot(aic_df, aes(stat, aic)) + 
  theme_bw() +
  xlab("Model") +
  ylab("AIC") +
  geom_jitter(alpha=0.3) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot morphology
mdata <- read.csv("data/syma_spp_morphology.csv")

# scale palette by PC
mdata$col <- ifelse(mdata$three_species_model == "megarhyncha", 
                          meg_val, ifelse(mdata$three_species_model == "torotoro",
                                          tor_val,och_val))

# function to map colors to viridis palette for discrete histogram
num_vec <- sample_data$PC1
map_viridis <- function(vec, num) {
  
  vector_expanded <-round(vec, 1) * 10 # expand to allow for decimal precision
  vector_exp_range <- max(vector_expanded) - min(vector_expanded)
  
  colour_vector <- viridis(vector_exp_range + 1) # get vector of colour values for all possible decimals between min and max value
  value_to_colour <- colour_vector[num * 10 - min(vector_expanded) + 1] # retrieve colour value for number
  
  return(value_to_colour)
  
}

map_viridis(num_vec, meg_val) # returns "#73D056FF"
map_viridis(num_vec, tor_val) # returns #443B84FF"
map_viridis(num_vec, och_val) # returns "#20938CFF"

# plot song
vdata <- read.csv("data/syma_spp_calls.csv")

g <- ggplot(vdata, aes(x=PC1, fill = two_species)) + 
  theme_bw() +
  scale_fill_manual(values = c("#73D056FF","#443B84FF"))+
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  geom_density(alpha=0.9)+
  theme(legend.position = "none") +
  xlab("PC1") +
  ylab("Density")

# plot morphology
h <- ggplot(mdata, aes(x=PC1, fill = three_species_model)) + 
  theme_bw() +
  scale_fill_manual(values = c("#73D056FF","#20938CFF","#443B84FF"))+
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  geom_density(alpha=0.9)+
  theme(legend.position = "none") +
  xlab("PC1") +
  ylab("Density")

# gridded w/ cowplot
top <- plot_grid(a, d, e, labels=c("A","B","C"), nrow = 1, ncol = 3, rel_widths = c(2.5,1,1), rel_heights = c(1,0.8,0.8))
#left <- plot_grid(d, e, labels=c("D","E"), nrow = 1, ncol = 2, rel_widths = c(1,1))
right <- plot_grid(g, h, labels=c("G","H"), nrow = 1, ncol = 2, rel_widths = c(1,1))
#legend_lm <- get_legend(e + theme(legend.position="bottom"))
#lm <- plot_grid(left, legend_lm, ncol = 1, rel_heights = c(1, .2))
legend_pheno <- get_legend(h + theme(legend.position="bottom"))
phylo <- plot_grid(b, c, labels=c("B","C"), ncol = 2)
pheno <- plot_grid(right, legend_pheno, ncol = 1, rel_heights = c(1, .2))
bottom <- plot_grid(phylo, pheno, ncol = 2, rel_widths = c(0.9,1.1))
joint <- plot_grid(top, bottom, nrow = 2, rel_heights = c(1,0.8))
