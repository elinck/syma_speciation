### analysis of Syma spp. morphology

install.packages("xlsx");install.packages("clustvarsel");install.packages("mclust")
install.packages("ggplot2")
install.packages("rJava")
library("rJava")
library(xlsx); library(ggplot2); library(clustvarsel);library(mclust);library(gridExtra);library(grid);


library(ggplot2)
# upload data
setwd("~/Dropbox/Syma/syma_morphology")
syma <- read.csv("Syma_data.csv")
#data<-syma[syma$Sex=="m",] #filter by sex
data<-syma[syma$include=="y",] #filter by NAs

symashade <- c("#228B22", "#FFFFA8", "#C6E01B")
measurements <- data[,14:19]
species <- as.character(data$three_species_model)
measurements <- cbind(measurements,species)
measurements <- measurements[!measurements$species=="ochracea",] #if dropping ochracea
data <- data[!data$three_species_model=="ochracea",]

torotoro <- measurements[measurements$species=="torotoro",]
megarhyncha <- measurements[measurements$species=="megarhyncha",]

t.test(torotoro$bill_from_nostril, megarhyncha$bill_from_nostril) #sig
t.test(torotoro$bill_width, megarhyncha$bill_width) #sig
t.test(torotoro$bill_depth, megarhyncha$bill_depth) #sig
t.test(torotoro$tarsus, megarhyncha$tarsus) #sig
t.test(torotoro$wing_chord, megarhyncha$wing_chord) #sig
t.test(torotoro$tail_length, megarhyncha$tail_length) #sig

p1 <- ggplot(measurements, aes(species,tarsus)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean tarsus length (mm)") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank()) + 
  theme(legend.position="none") +
  theme(legend.position="none")

p2 <- ggplot(measurements, aes(species,wing_chord)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean wing chord (mm)") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.position="none")

p3 <- ggplot(measurements, aes(species, bill_from_nostril)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean bill length (mm)") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.position="none")

p4 <- ggplot(measurements, aes(species, bill_width)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean bill width (mm)") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank()) + 
  theme(legend.position="none")

p5 <- ggplot(measurements, aes(species, bill_depth)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean bill depth (mm)") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank()) + 
  theme(legend.position="none")

p6 <- ggplot(measurements, aes(species, tail_length)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean tail length (mm)") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.position="none")

#load grid / shared legend function
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(p1,p4,p5,p2,p3,p6)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position=element_blank()))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

# plot all together w/ shared legend
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)


pca<-prcomp(data[,14:19])
pca$rotation # PC loadings
summary(pca) # variance by axis
pca$x  # scores for each individual on each PC axis
data<-cbind(pca$x, data) # combining the PC scores with the old dataframe 
useful <- clustvarsel(pca$x, verbose = TRUE) # test to determine which PCs are most useful for group discrimination in NMMs

#need to examine differences in males / females
par(mfrow=c(1,1))





# evaluate the two-species model w/ clustering, real data
two_species <- data$english_name
two_species <- two_species[c(1:60)] #if dropping ochracea
table(two_species)
measurements <- measurements[,-c(7)]
clPairs(measurements, two_species)
BIC <- mclustBIC(measurements)
plot(BIC)
summary(BIC)
mod1 <- Mclust(measurements, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification",color=symashade)
table(mod1, mod1$classification)

# plot pairs of traits
par(mfrow = c(2,2))
plot(mod1, what = "uncertainty", dimens = c(2,1), main = "",col=symashade)
plot(mod1, what = "uncertainty", dimens = c(3,1), main = "",col=symashade)
plot(mod1, what = "uncertainty", dimens = c(2,3), main = "",col=symashade)
par(mfrow = c(1,1))

# evaluate the two-species model w/ clustering, recommended PC scores
PCs <- data[,1:2]
table(two_species)
clPairs(PCs, two_species, color=symashade)
BIC2 <- mclustBIC(PCs)
plot(BIC2)
summary(BIC2)
mod1pc <- Mclust(PCs, x = BIC)
summary(mod1pc, parameters = TRUE)
plot(mod1pc, what = "classification",col=symashade)
legend(4, -6, c("megarhyncha","torotoro"),fill=symashade)

table(mod1pc, mod1pc$classification)

par(mfrow = c(2,2))
plot(mod1pc, what = "uncertainty", dimens = c(2,1), main = "",col=symashade)
plot(mod1pc, what = "uncertainty", dimens = c(3,1), main = "",col=symashade)
plot(mod1pc, what = "uncertainty", dimens = c(2,3), main = "",col=symashade)
par(mfrow = c(1,1))

# violin plot
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
  

# evaluate the three-species model w/ clustering, real data
three_species <- data$three_species_model
table(three_species)
clPairs(measurements[,1:4], three_species,colors = symashade)
BIC3 <- mclustBIC(measurements)
plot(BIC3)
summary(BIC3)
mod2 <- Mclust(nums, x = BIC)
summary(mod2, parameters = TRUE)
plot(mod2, what = "classification")
table(three_species, mod2$classification)

par(mfrow = c(2,2))
plot(mod2, what = "uncertainty", dimens = c(2,1), main = "",col=symashade)
plot(mod2, what = "uncertainty", dimens = c(3,1), main = "",col=symashade)
plot(mod2, what = "uncertainty", dimens = c(2,3), main = "",col=symashade)
par(mfrow = c(1,1))

# evaluate the three-species model w/ clustering, PC scores
clPairs(PCs, three_species, colors = symashade)
BIC4 <- mclustBIC(PCs)
plot(BIC4)
summary(BIC4)
mod2pc <- Mclust(PCs, x = BIC4)
summary(mod2pc, parameters = TRUE)
plot(mod2pc, what = "classification",col=symashade)
table(three_species, mod2pc$classification)

cluster <- mod1pc$classification
ggplot(PCs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(aes(x=PC1, y=PC2, color=factor(cluster), shape=three_species), size=2) +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  scale_fill_manual(values = symashade) +
  scale_color_manual(values = symashade) +
  theme(legend.title=element_blank()) +
  guides(color=guide_legend("cluster"),fill=guide_legend("cluster"))


mod2 = MclustDA(X, class, modelType = "EDDA")
summary(mod2)

par(mfrow = c(2,2))
plot(mod2pc, what = "uncertainty", dimens = c(2,1), main = "",col=symashade)
plot(mod2pc, what = "uncertainty", dimens = c(3,1), main = "",col=symashade)
plot(mod2pc, what = "uncertainty", dimens = c(2,3), main = "",col=symashade)
par(mfrow = c(1,1))

# Density Estimates
wc <- densityMclust(data$wing_chord) 
summary(wc)
plot(wc, what = "BIC")
d1 <- plot(wc, what = "density", data = data$wing_chord, breaks = 15) #bimodal 

tar <- densityMclust(data$tarsus) 
summary(tar)
plot(tar, what = "BIC")
d2 <- plot(tar, what = "density", data = data$tarsus, breaks = 15) #normal

bfn <- densityMclust(data$bill_from_nostril) 
summary(bfn)
plot(bfn, what = "BIC")
d3 <- plot(bfn, what = "density", data = data$bill_from_nostril, breaks = 15) #normal

bw <- densityMclust(data$bill_width) 
summary(bw)
plot(bw, what = "BIC")
d4 <- plot(bw, what = "density", data = data$bill_width, breaks = 15) #bimodal

bd <- densityMclust(data$bill_depth)
summary(bd)
plot(bd, what = "BIC")
d5 <- plot(bd, what = "density", data = data$bill_depth, breaks = 15) #bimodal(ish)

tl <- densityMclust(data$tail_length)
summary(tl)
plot(tl, what = "BIC")
d6 <- plot(tl, what = "density", data = data$tail_length, breaks = 15) #unimodal

pc1d <- densityMclust(data$PC1)
summary(mpc1d)
plot(pc1d, what = "BIC")
d7 <- plot(pc1d, what = "density", data = data$PC1, breaks = 15) #bimodal

pc2d <- densityMclust(data$PC2)
summary(pc2d)
plot(pc2d, what = "BIC")
d8 <- plot(pc2d, what = "density", data = data$PC2, breaks = 15) #unimodal

par(mfrow=c(2,3))
plot(wc, what = "density", data = data$wing_chord, breaks = 15) #bimodal 
plot(tar, what = "density", data = data$tarsus, breaks = 15)
plot(bfn, what = "density", data = data$bill_from_nostril, breaks = 15)
plot(bw, what = "density", data = data$bill_width, breaks = 15) #bimodal
plot(bd, what = "density", data = data$bill_depth, breaks = 15) #bimodal(ish)
plot(tl, what = "density", data = data$tail_length, breaks = 15)
par(mfrow=c(1,1))

#multivariate, w/ PCs 
mv <- densityMclust(data[,1:2])
summary(mv)
plot(mv, what = "BIC")
plot(mv, what = "density")
plot(mv, what = "density", type = "image", col = "symashade")
plot(mv, what = "density", type = "persp")

### older clustering analysis 
scores <- pca[,1:25] # select first 3 PCs
km <- kmeans(scores, 2)# perform Kmeans clustering
km3 <- kmeans(scores, 3)# perform Kmeans clustering
scores <- cbind(scores, data2$Scientific_Name)

ggdata <- data.frame(scores, cluster=km$cluster, species=data2$English.Name) # create df
cluster <- as.factor(cluster) #transform to factors
cluster <- revalue(cluster, c("1"="Syma megarhyncha","2"="Syma torotoro")) #rename cluster levels

symashade <- c("#228B22","#C6E01B","#FFC400")
ggplot(ggdata) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(aes(x=PC1, y=PC2, color=factor(cluster),shape=species)) +
  scale_fill_manual(values = symashade) +
  scale_color_manual(values = symashade) +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  theme(legend.title=element_blank()) +
  guides(color=guide_legend("cluster"),fill=guide_legend("cluster"))
