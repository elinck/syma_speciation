### analysis of acoustic divergence between S. torotoro and S. megarhyncha

install.packages("warbleR")
install.packages("ggfortify")
install.packages("cluster")
install.packages("rgl")
library(rgl)
library(warbleR); library(ggfortify); library(ggplot2); library(cluster); library(plyr);library(mclust)
source("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R") 
setwd("~/Dropbox/Syma/song_analysis/wavs")

# download xeno canto material 
symaxc <- querxc(qword = "Syma", download = FALSE) 
View(symaxc)
symasong <- symaxc[grep("song", symaxc$Vocalization_type, ignore.case = TRUE),]
querxc(X = symasong) 
write.csv(symasong, "symasongxc.csv", row.names = FALSE)
mp32wav() # convert mp3s to wavs
checkwavs() # all files are readable
wavs <- list.files(pattern="wav$") # create list of wavs
sub <- wavs[c(27)] # subset for testing 
lspec(flist = sub, flim = c(0,3), sxrow = 10, rows = 15, ovlp = 10, it = "tiff") # test parameters
lspec(flim = c(1.5, 11), sxrow = 15, rows = 20, ovlp = 10, it = "tiff") # produce spectrograms

# manualoc(wl = 512, flim = c(0,4), seltime = 1, tdisp = NULL, reccomm =
#           FALSE, wn = "hanning", title = TRUE, selcomm = FALSE, osci = FALSE, player =
#           NULL, pal = reverse.gray.colors.2, path = NULL, flist = NULL) #manual detection params

syma.ad <- autodetec(flist = wavs, bp = c(1, 3), threshold = 4, mindur = 0.5, maxdur = 15, envt="abs",
                     ssmooth = 1500, ls = TRUE, res = 100, 
                     flim = c(1, 3), wl = 512, set =TRUE, sxrow = 15, rows = 20, 
                     redo = FALSE, it = "tiff") #best parameters 040417
write.csv(syma.ad, file = "syma_selections_all.csv")
table(syma.ad$sound.files) # view selections by individuaL 
set.seed(5) #set seed
X <- syma.ad[sample(1:nrow(syma.ad),(nrow(syma.ad)*0.1)), ] #test on 10% of selections
snrspecs(X = syma.ad, flim = c(1, 4), snrmar = 0.2, it = "tiff") # print noise margins

syma.ad <- read.csv("syma_selections.csv")
syma.snr <- sig2noise(X = syma.ad, mar = 0.2) # adjust noise to signal ratio 

syma.hisnr <- syma.snr[ave(-syma.snr$SNR, syma.snr$sound.files, FUN = rank) <= 5,] # select the 5 most high quality songs per individual 
table(syma.hisnr$sound.files) # view no. of high quality selections per individual
write.csv(syma.hisnr, file = "syma_hiquality.csv") # save as CSV
syma.hisnr <- read.csv("syma_hiquality.csv")
#params <- specan(syma.hisnr, bp = c(1,4), threshold = 4)# calculate acoustic params
params <- read.csv("params.csv")

species <- params$sp. 
duration <- params$duration
freq <- params$meanfreq
basic <- as.data.frame(stringsasfactors=F,cbind(as.character(species),as.numeric(duration),as.numeric(freq)))
#symashade <- c("#FFC400", "#337058", "#030303")
symashade <- c("#228B22","#C6E01B")

#load grid / shared legend function
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(p1,p2)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
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

p1 <- ggplot(basic, aes(species,freq)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean frequency of call element") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.position="none")

p2 <- ggplot(basic, aes(species,duration)) +
  scale_fill_manual(values = symashade) +
  geom_boxplot(aes(fill = factor(species))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean duration of call element") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.title=element_blank())+
  theme(legend.position="none")

# plot all together w/ shared legend
grid_arrange_shared_legend(ncol=2,nrow=1,p1,p2)

grid.arrange(p1,p2,ncol=2)

pca <- prcomp(x = params[, sapply(params, is.numeric)], scale. = TRUE) #perform PCA
summary(pca) #summarize PCA scores
scores <- pca$x[,1:3] # select first 3 PCs
params <- cbind(scores,params)
species <- params$sp.
table(species)
clPairs(scores,species,colors = symashade)

par(mfrow=c(1,2))

freq <- as.vector(freq)
mod1 <- densityMclust(freq)
duration <- as.vector(duration)
mod2 <- densityMclust(duration)

plot(mod1, what = "density", data = freq, breaks = 15)
plot(mod2, what = "density", data = duration, breaks = 15)
plot(mod1, what = "BIC")
plot(mod2, what = "BIC")

par(mfrow=c(1,1))


PC1 <- as.vector(params$PC1)
mod3 <- densityMclust(PC1)
plot(mod3, what = "BIC")
plot(mod3, what = "density", data = PC1, breaks = 15)


# loop to assign each individual a species binomial
names <- vector()
for (i in 1:nrow(params)){
  if (grepl("mega", params[i,1])){
    names[i] = "Syma megarhyncha"
  }
  else {names[i] = "Syma torotoro"
  }
}
params <- cbind(params,names)

pca <- prcomp(x = params[, sapply(params, is.numeric)], scale. = TRUE) #perform PCA
summary(pca) #summarize PCA scores
scores <- pca$x[,1:3] # select first 3 PCs
km <- kmeans(scores, 2) # perform Kmeans clustering
ggdata <- data.frame(scores, cluster=km$cluster, species) # create df
cluster <- as.factor(ggdata$cluster) #transform to factors
cluster <- revalue(cluster, c("1"="Syma megarhyncha","2"="Syma torotoro")) #rename cluster levels

symashade <- c("#228B22","#C6E01B")

# violin plot
m1 <- ggplot(ggdata, aes(species, PC1, fill = species, guide=F)) +
  theme_light() +
  geom_point(size=.2,col="grey",position="jitter")+
  geom_violin(draw_quantiles = c(0.5)) +
  guides(fill=guide_legend(title="Species"))+
  scale_fill_manual(labels=c("megarhyncha","torotoro"), values=symashade) + 
  theme(strip.background = element_rect(colour="black", fill="grey100")) +
  theme(axis.title.x=element_blank())

#histogram
 m2 <- ggplot(ggdata, aes(x=PC1, fill = species)) + 
  geom_histogram(aes(y=..density..), alpha=0.2, 
                 position="dodge",color="black", bins = 30)+
  theme_light() +
  coord_flip() +
  guides(fill=guide_legend(title="Species"))+
  scale_fill_manual(labels=c("megarhyncha","torotoro"), values=symashade) + 
  geom_density(alpha=0.9)+
  xlab("PC1") +
  ylab("Density")

pdf("~/Dropbox/remark_presentations/figures/vocalizations.pdf", width = 5, height = 7)
grid.arrange(m1,m2,ncol=1)
dev.off()

#plot first two PCs w/ 95% confidence ellipses 
ggplot(ggdata) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(aes(x=PC1, y=PC2, color=factor(cluster), shape=species), size=2) +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  scale_fill_manual(values = c("#C6E01B","#228B22")) +
  scale_color_manual(values = c("#C6E01B","#228B22")) +
  theme(legend.title=element_blank()) +
  guides(color=guide_legend("cluster"),fill=guide_legend("cluster"))

useful <- clustvarsel(pca$x, verbose = TRUE)
#keep <- c("PC8", "PC16", "PC11", "PC19", "PC2", "PC13", "PC1", "PC3", "PC25", "PC17", "PC18", "PC14", "PC15", "PC4", "PC7",
"PC10")
keep <- c("PC12", "PC2", "PC5")
object <- pca$x
object <- object[,keep]
species <- params$sp.
table(species)
clPairs(object, species)
BIC <- mclustBIC(object)
plot(BIC)
summary(BIC)
mod1pc <- Mclust(object, x = BIC)
summary(mod1pc, parameters = TRUE)
plot(mod1pc, what = "classification",col=symashade)
table(mod1pc, mod1pc$classification)

par(mfrow = c(1,1))

pca# set up data for PC1 / frequency plots  
meanfreq <- params$meanfreq 
scores_df <- as.data.frame(scores)
pc1 <- scores_df$PC1
pc2 <- scores_df$PC2
pc3 <- scores_df$PC3

# run linear models of PC1 scores and mean frequency of calls
pc1_freq_lm <- lm(pc1 ~ meanfreq)
summary(pc1_freq_lm)
pc2_freq_lm <- lm(pc2 ~ meanfreq)
summary(pc2_freq_lm)
pc3_freq_lm <- lm(pc3 ~ meanfreq)
summary(pc3_freq_lm)
plotdata <- cbind(meanfreq,scores_df,names)

# set syma palette pref
symashade <- c("#FFC400", "#337058", "#030303")

# plot PC1 by frequency
p1 <- ggplot(plotdata, aes(x=meanfreq, y=pc1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  geom_point(aes(color=factor(names)), size=2, shape=20) +
  xlab("Mean frequency of song element") +
  ylab("PC1") +
  theme(axis.title.y = element_text(size=10)) +
  theme(axis.title.x = element_text(size=10)) +
  theme(legend.title=element_blank()) +
  guides(color=guide_legend("names"),fill=guide_legend("names"))

# boxplots of meanfreq by species
ggplot(plotdata, aes(x=factor(names), y=meanfreq)) + 
  geom_boxplot(aes(fill = factor(names))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=FALSE) +
  xlab("Species") +
  ylab("Mean frequency of call element") +
  theme(axis.title.y = element_text(size=10)) +
  theme(axis.title.x = element_text(size=10)) +
  theme(legend.title=element_blank())

freqd <- densityMclust(as.vector(basic$freq)) 
summary(wc)
plot(wc, what = "BIC")
plot(wc, what = "density", data = data$wing_chord, breaks = 15, xlab = "Wing Chord") #bimodal 


### conventional plotting code
#pcascor <- as.data.frame(pca[[5]]) #transform to dataframe 
#pcascor <- cbind(pcascor,names)
#plot(pcascor[, 1], pcascor[, 2], pch = 20, col = pcascor$names,
#     cex = 1, xlab = "PC1", ylab = "PC2",)

