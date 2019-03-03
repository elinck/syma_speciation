### plots for linck et al. 2019

library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(raster)
library(rgdal)
library(ggtree)
library(cowplot)
library(ape)
library(png)
library(gtable)
library(grid)
library(scales)
library(phytools)
library(cowplot)
library(ape)
library(phangorn)
library(ggridges)
library(data.table)

# set working directory
setwd("~/Dropbox/syma_speciation/")

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
  geom_point(data=sample_data,aes(x=long,y=lat,fill=PC1),pch=21,size=sample_data$size*2.5, alpha=0.85,stroke=1) #+
  #scale_size(guide = "none") +
  #annotate("segment", x = 147, xend = 149, y = -5.6, yend = -1) +
  #annotate("text", x = 149, y = -0.5, label = "sellamontis") +
  #annotate("segment", x = 150.5, xend = 153, y = -9.3, yend = -8) +
  #annotate("text", x = 153, y = -7.5, label = "ochracea") +
  #annotate("segment", x = 131.5, xend = 140, y = 0, yend = 1) +
  #annotate("segment", x = 144, xend = 140, y = -3.3, yend = 1) +
  #annotate("text", x = 140, y = 1.6, label = "torotoro") +
  #annotate("segment", x = 147, xend = 148.5, y = -9.5, yend = -12) +
  #annotate("text", x = 149, y = -12.3, label = "meeki") +
  #annotate("segment", x = 142, xend = 140, y = -9.4, yend = -10) +
  #annotate("segment", x = 143, xend = 140, y = -12.6, yend = -10) +
  #annotate("text", x = 137.8, y = -10, label = "flavirostris") +
  #annotate("segment", x = 144.5, xend = 144.5, y = -7.75, yend = -9.4) +
  #annotate("text", x = 144.5, y = -9.75, label = "pseuestes") +
  #annotate("segment", x = 134.2, xend = 134.2, y = -7, yend = -8.75) +
  #annotate("text", x = 134.2, y = -9, label = "tentelare") +
  #annotate("segment", x = 134.5, xend = 132.5, y = -4.1, yend = -4.5) +
  #annotate("text", x = 131, y = -4.5, label = "wellsi") +
  #annotate("segment", x = 136.5, xend = 136.3, y = -4.8, yend = -6.75) + 
  #annotate("text", x = 136.5, y = -7.0, label = "pseuestes") +
  #annotate("segment", x = 145, xend = 145, y = -5.5, yend = -0) +
  #annotate("text", x = 145, y = 0.5, label = "megarhyncha")

# pca
b <- ggplot(data=sample_data,aes(x=PC1,y=PC2,col=sp)) + 
  geom_text(aes(label=prep_ID, color=sp)) +
  scale_color_manual(values = c("#73D056FF","#20938CFF","#443B84FF"),guide=FALSE)+
  theme_bw() +
  ylim(-29,10)

  
c <- readPNG("~/Dropbox/syma_speciation/figures/syma.png")
c <- rasterGrob(c)

# plot admixture
d <- readPNG("~/Dropbox/syma_speciation/figures/structure.png")
d <- rasterGrob(d)
pdf("figures/final/S2_admix.pdf",width=6,height=3)
plot_grid(d)
dev.off()

# plot morphology
mdata <- read.csv("data/syma_spp_morphology.csv")

# scale palette by PC
mdata$col <- ifelse(mdata$three_species_model == "megarhyncha", 
                          meg_val, ifelse(mdata$three_species_model == "torotoro",
                                          tor_val,och_val))
#drop ochracea
mdata <- mdata[!mdata$three_species_model=="ochracea",]

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
e <- ggplot(vdata, aes(x=PC1, fill = two_species)) + 
  theme_bw() +
  scale_fill_manual(values = c("#73D056FF","#443B84FF"))+
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  geom_density(alpha=0.9)+
  theme(legend.position = "none") +
  xlab("PC1") +
  ylab("Density") 

# plot morphology
f <- ggplot(mdata, aes(x=PC1, fill = three_species_model)) + 
  theme_bw() +
  scale_fill_manual(values = c("#73D056FF","#443B84FF"))+
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  geom_density(alpha=0.9)+
  theme(legend.position = "none") +
  xlab("PC1") +
  ylab("Density")

# mtDNA tree
mtdna.tree <- read.tree("data/syma_RAxML_bipartitions.tre")
g <- ggtree(mtdna.tree) + 
  #geom_tiplab(color='black') +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  xlim(0, 0.2) +
  geom_hilight(57, fill="#20938CFF") +
  geom_hilight(50, fill="#73D056FF") +
  geom_hilight(32, fill="#443B84FF")
  #geom_cladelabel(57, "ochracea", offset = 0.0,barsize=0) +
  #geom_cladelabel(50, "megarhyncha", offset = 0.0,barsize=0) +
  #geom_cladelabel(32, "torotoro", offset = 0.0,barsize=0)

#plot densitree
trees <- read.nexus("raw_data/densitree.input.nex")
trees[[1]]$tip.label # show order of tip labels
tip_colors <- c(rep("#443B84FF",24),rep("#20938CFF",2),rep("#73D056FF",10))
densiTree(trees, scale.bar = FALSE, font = 2, width = 1.75, tip.color = tip_colors, cex = 0.75)

h <- readPNG("~/Dropbox/syma_speciation/figures/densitree.png")
h <- rasterGrob(h)

top <- plot_grid(a, b, c, labels = "AUTO", align = 'h', label_size = 12, nrow = 1, rel_widths = c(2.5,1,1))
bottom <- plot_grid(e, f, g, h, labels = c("D","E","F","G"), align = 'h', label_size = 12, nrow = 1)
pdf("figures/final/Fig_2_species_limits.pdf",width=11,height=6)
plot_grid(top, bottom, align = 'h', rel_heights = c(1.25,1),nrow=2)
dev.off()

#plot abbababa
dstat <- read.csv("~/Dropbox/syma_speciation/data/d_stat_results.csv")
dtable <- read.csv("~/Dropbox/syma_speciation/data/d_summary.csv")
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
mag <- map2color(dtable$prop_sig,magma(8))
i <- ggplot(data=dstat,aes(x=abs(Z),y=tree,fill=tree))+
  theme_bw() +
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=8))+
  scale_fill_manual(values=mag,guide=FALSE)+
  ylab("Test")+xlab("abs(Z)")+
  xlim(-5,100)+
  #geom_vline(aes(xintercept=4.3),col="red",lwd=0.5)+
  geom_vline(aes(xintercept=3.59),col="black",linetype="dashed",lwd=0.5)+
  geom_density_ridges(scale=1.25,lwd=0.35)+
  geom_point(data=dtable,shape=8,size=.75,aes(x=100,y=as.integer(factor(tree))+.3,col=median_Z>3.59))+
  scale_color_manual(values=rep(c("white","black"),nrow(dtable)),guide=F)

#plot moments
ll.models <- read.csv("data/moments_plotting.csv")
j <- ggplot(ll.models, aes(x=model, y=ll_model)) + 
  theme_bw()+
  theme(axis.title=element_text(size=12)) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  #geom_jitter(shape=16, position=position_jitter(0.2), color="#8A008A") +
  geom_boxplot(alpha=0.8) +
  ylab("Log-likelihood") +
  xlab("Model") +
 ylim(-800,-500)

#plot treemix
source("~/Desktop/treemix/src/plotting_funcs.R")
#pdf(file = "treemix.pdf", width = 5, height = 5)
#pdf("figures/treemix.pdf",width=5,height=5)
#plot_tree("raw_data/m3") # mL migration surface
#dev.off()
k <- readPNG("~/Dropbox/syma_speciation/figures/treemix.png")
k <- rasterGrob(k)
pdf("figures/final/S3_treemix.pdf",width=5,height=5)
plot_grid(k)
dev.off()

#plot smcpp
smc <- read.csv("data/smcpp_df.csv",stringsAsFactors = F)
#pdf("figures/fig3.pdf",width=3.25,height=4)
l <- ggplot(data=smc,aes(x=x,y=y))+facet_grid(species~.)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        strip.text=element_text(size=12),
        axis.text=element_text(size=8),
        panel.grid.major=element_line(color="grey",size=0.25),
        strip.background = element_blank())+
  xlab("Generations")+ylab("Ne")+
  ylim(0,4.5e4)+
  xlim(2e3,5.5e4)+
  #scale_color_manual(values=c("#443B84FF","#73D056FF"),guide=F)+
  geom_rect(aes(xmin=9000,xmax=13250,ymin=0,ymax=4.5e4),fill="grey90",alpha=0.75)+
  geom_path(lwd=.35,alpha=0.8, aes(color=species))+
  scale_color_manual(values=c("#73D056FF","#443B84FF"),guide=F)
#dev.off()

# plot network
m <- readPNG("~/Dropbox/syma_speciation/figures/phylonet.png")
m <- rasterGrob(m)

pdf("figures/final/Fig_3_introgression.pdf",width=8,height=7)
plot_grid(i,m,j,l, labels = "AUTO", align = 'h', label_size = 12, nrow = 2)
dev.off()

# manhattan plots
win.df <- read.csv("data/window_stats_chr.csv")[,-1]
chr_labels <- read.csv("data/chr_labels.csv")
#quartz()
#png(width=6,height=1.5,units="in",res=600,file="~/Desktop/both.png")
n <- ggplot(data=win.df,aes(x=row,y=value,col=chr))+theme_bw()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank())+  
  scale_y_continuous(breaks=seq(0,1.0,0.2),minor_breaks = NULL)+
  geom_point(data=subset(win.df,value>0.0),size=1,shape=16,alpha=0.7,aes(color=factor(col)))+
  geom_point(data=subset(win.df, outlier=="TRUE"),size=1,shape=16,col="red")+
  scale_color_manual(values=rep(c("grey70","grey85"),length(levels(factor(win.df$chr)))/2+1))+
  geom_hline(data=subset(win.df, stat=="F[ST]"),aes(yintercept=0.409031),linetype="dashed",lwd=0.25)+
  geom_hline(data=subset(win.df, stat=="D[XY]"),aes(yintercept=0.1727),linetype="dashed",lwd=0.25)+
  geom_hline(data=subset(win.df, stat=="pi [megarhyncha]"),aes(yintercept=0.1701),linetype="dashed",lwd=0.25)+
  geom_hline(data=subset(win.df, stat=="pi [torotoro]"),aes(yintercept=0.1527655),linetype="dashed",lwd=0.25)+
  #geom_point(data=subset(win.df,value>=quantile(win.df$value,0.999)),shape=16,stroke=0.9,size=0.5,col="red")+
  geom_line(aes(y=rollmean),lwd=0.15,col="black")+
  geom_segment(data=chr_labels,aes(x=start+25,xend=stop-25,y=-0.05,yend=-0.05,col=NA),col="black")+
  geom_text(data=chr_labels,aes(label=chr,x=mid,y=-0.1,col=NA),
                  col="black",size=2.25,angle=0,direction="y",box.padding = 0.15,
                  segment.size=0.2) +
  facet_grid(stat ~ ., labeller=label_parsed, scales="free_y") +
  theme(strip.background = element_rect(colour="black", fill="grey100"))

#compare selective sweeps
sweeps.df <- read.csv("data/sweeps.csv")
sweeps.top <- sweeps.df[sweeps.df$mu>0.004681838,] #FPR from simulation
pi.df <- read.csv("data/pi.csv")

o <- ggplot(sweeps.top, aes(x=stat, y=mu)) + 
  theme_bw()+
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_boxplot(alpha=0.8) +
  ylab(expression(mu)) +
  xlab("Genomic subset") +
  theme(axis.text.x = element_text(size=6))

p <- ggplot(sweeps.top, aes(stat)) +
  theme_bw()+
  geom_bar(color="black",fill = c("gray70","white","gray70","white"),width=0.75) +
  ylab("Number of candidate sweeps") +
  xlab("Genomic subset") +
  theme(axis.text.x = element_text(size=6))

q <- ggplot(pi.df, aes(x=PI, fill=dataset)) +
  theme_bw() +
  geom_density() +
  scale_fill_viridis(discrete=TRUE,labels=c("megarhyncha","simulation")) +
  theme(legend.position = c(0.7, 0.8),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  xlab(~paste(pi)) +
  ylab("Count") +
  xlim(0,0.004)

top <- plot_grid(n, align = 'h', labels = c("A"), label_size = 12, nrow = 1,label_y = 0.97)
bottom <- plot_grid(o, p, q, align = 'h', labels = c("B","C","D"), label_size = 12, nrow = 1)
pdf("figures/final/Fig_4_selection.pdf",width=9,height=9,compress = TRUE)
plot_grid(top, bottom, nrow=2, rel_heights = c(0.65,0.35))
dev.off()

# look at individual chromosomes
win.5 <- subset(win.df, chr==5 | chr==11 | chr==9 | chr==23)
win.5 <- subset(win.5, stat=="F[ST]")
win.5$chr_order <- factor(win.5$chr, levels=c('5','9','11','23'))
ylab.text = expression('F'['ST'])
r <- ggplot(data=win.5,aes(x=row,y=value))+theme_bw()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1))+  
  ylab(ylab.text) +
  xlab("Position") +
  scale_y_continuous(breaks=seq(0,1.0,0.2),minor_breaks = NULL)+
  geom_point(data=subset(win.5,value>0.0 & outlier=="FALSE"),size=1.25,shape=16,alpha=0.7,color="grey30")+
  geom_point(data=subset(win.5, outlier=="TRUE"),size=1.5,shape=16,col="red")+
  geom_hline(aes(yintercept=0.409031),linetype="dashed",lwd=0.25)+
  facet_grid(. ~ chr_order, scales="free_x") +
  theme(strip.background = element_rect(colour="black", fill="grey100"))

df2 <- read.csv("data/window_correlations.csv")[,-1]
s <- ggplot(data=df2,aes(x=dxy_meg_toro,y=Fst_meg_toro)) +
  theme_bw()+
  theme(axis.line.x=element_line(),
        axis.text.x = element_text(size=6,angle=45,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=8),
        axis.ticks.x=element_line(),
        strip.background = element_blank(),
        strip.text=element_text(size=7),
        legend.background = element_blank(),
        legend.text = element_text(size=7),
        legend.title=element_text(size=7),
        axis.line.y=element_line())+
  stat_bin_hex(bins = 30,aes(fill=log(..count..)))+
  scale_fill_distiller(palette = "BrBG",direction = -1,name="log(n. windows)",guide=FALSE)+
  #geom_point(shape=21,stroke=0.3,size=0.4,alpha=0.3)+
  xlab(expression(D[XY])) +
  ylab(expression(F[ST])) +
  geom_smooth(method="lm",fill=NA,col="black",lwd=0.65)

t <- ggplot(data=df2,aes(x=pi_meg,y=Fst_meg_toro)) +
  theme_bw()+
  theme(axis.line.x=element_line(),
        axis.text.x = element_text(size=6,angle=45,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=8),
        axis.ticks.x=element_line(),
        strip.background = element_blank(),
        strip.text=element_text(size=7),
        legend.background = element_blank(),
        legend.text = element_text(size=7),
        legend.title=element_text(size=7),
        axis.line.y=element_line())+
  stat_bin_hex(bins = 30,aes(fill=log(..count..)))+
  scale_fill_distiller(palette = "BrBG",direction = -1,name="log(n. windows)",guide=FALSE)+
  #geom_point(shape=21,stroke=0.3,size=0.4,alpha=0.3)+
  xlab(~paste(pi[megarhyncha])) +
  ylab(expression(F[ST])) +
  geom_smooth(method="lm",fill=NA,col="black",lwd=0.65)

u <- ggplot(data=df2,aes(x=pi_toro,y=Fst_meg_toro)) +
  theme_bw()+
  theme(axis.line.x=element_line(),
        axis.text.x = element_text(size=6,angle=45,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=8),
        axis.ticks.x=element_line(),
        strip.background = element_blank(),
        strip.text=element_text(size=7),
        legend.background = element_blank(),
        legend.text = element_text(size=7),
        legend.title=element_text(size=7),
        axis.line.y=element_line())+
  stat_bin_hex(bins = 30,aes(fill=log(..count..)))+
  scale_fill_distiller(palette = "BrBG",direction = -1,name="log(n. windows)",guide=FALSE)+
  #geom_point(shape=21,stroke=0.3,size=0.4,alpha=0.3)+
  xlab(~paste(pi[torotoro])) +
  ylab(expression(F[ST])) +
  geom_smooth(method="lm",fill=NA,col="black",lwd=0.65)

v <- ggplot(data=df2,aes(x=pi_meg,y=pi_toro)) +
  theme_bw()+
  theme(axis.line.x=element_line(),
        axis.text.x = element_text(size=6,angle=45,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=8),
        axis.ticks.x=element_line(),
        strip.background = element_blank(),
        strip.text=element_text(size=7),
        legend.background = element_blank(),
        legend.text = element_text(size=7),
        legend.title=element_text(size=7),
        axis.line.y=element_line())+
  stat_bin_hex(bins = 30,aes(fill=log(..count..)))+
  scale_fill_distiller(palette = "BrBG",direction = -1,name="log(n. windows)",guide=FALSE)+
  #geom_point(shape=21,stroke=0.3,size=0.4,alpha=0.3)+
  xlab(~paste(pi[torotoro])) +
  ylab(~paste(pi[megarhyncha])) +
  geom_smooth(method="lm",fill=NA,col="black",lwd=0.65)

w <- ggplot(data=df2,aes(x=pi_meg,y=pi_toro)) +
  theme_bw()+
  theme(axis.line.x=element_line(),
        axis.text.x = element_text(size=6,angle=45,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=8),
        axis.ticks.x=element_line(),
        strip.background = element_blank(),
        strip.text=element_text(size=7),
        legend.background = element_blank(),
        legend.text = element_text(size=7),
        legend.title=element_text(size=7),
        axis.line.y=element_line())+
  stat_bin_hex(bins = 30,aes(fill=log(..count..)))+
  scale_fill_distiller(palette = "BrBG",direction = -1,name="log(n. windows)",guide=FALSE)+
  #geom_point(shape=21,stroke=0.3,size=0.4,alpha=0.3)+
  xlab(expression(D[XY])) +
  ylab(~paste(pi[megarhyncha])) +
  geom_smooth(method="lm",fill=NA,col="black",lwd=0.65)

x <- ggplot(data=df2,aes(x=pi_meg,y=pi_toro)) +
  theme_bw()+
  theme(axis.line.x=element_line(),
        axis.text.x = element_text(size=6,angle=45,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=8),
        axis.ticks.x=element_line(),
        strip.background = element_blank(),
        strip.text=element_text(size=7),
        legend.background = element_blank(),
        legend.text = element_text(size=7),
        legend.title=element_text(size=7),
        axis.line.y=element_line())+
  stat_bin_hex(bins = 30,aes(fill=log(..count..)))+
  scale_fill_distiller(palette = "BrBG",direction = -1,name="log(n. windows)",guide=FALSE)+
  #geom_point(shape=21,stroke=0.3,size=0.4,alpha=0.3)+
  xlab(expression(D[XY])) +
  ylab(~paste(pi[torotoro])) +
  geom_smooth(method="lm",fill=NA,col="black",lwd=0.65)

y <- ggplot(df2, aes(df2$dxy_meg_toro)) + 
  theme_bw()+
  xlim(0,0.5)+
  xlab(expression(D[XY])) +
  theme(axis.title.y=element_blank()) +
  geom_histogram(bins=100)

z <- ggplot(df2, aes(df2$Fst_meg_toro)) + 
  theme_bw()+
  xlim(0,0.5)+
  xlab(expression(F[ST])) +
  theme(axis.title.y=element_blank()) +
  geom_histogram(bins=100)

a1 <- ggplot(df2, aes(df2$pi_meg)) + 
  theme_bw()+
  xlim(0,0.5)+
  xlab(~paste(pi[megarhyncha])) +
  theme(axis.title.y=element_blank()) +
  geom_histogram(bins=100)

b1 <- ggplot(df2, aes(df2$pi_toro)) + 
  theme_bw()+
  xlim(0,0.5)+
  xlab(~paste(pi[torotoro])) +
  theme(axis.title.y=element_blank()) +
  geom_histogram(bins=100)

top <- plot_grid(q,align = 'h', label_size = 12, nrow = 1)
middle <- plot_grid(s, t, u, v, w, x, align = 'h', label_size = 12, nrow = 1)
bottom <- plot_grid(y, z, a1, a2, align = 'h', label_size = 12, nrow = 1)
pdf("figures/final/S4_genome_scans.pdf",width=9,height=7)
plot_grid(top, middle, bottom, align = 'h', labels="AUTO", label_size = 12, nrow = 3, rel_heights = c(1.25,1,1))
dev.off()
