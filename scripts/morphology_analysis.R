### analysis of Syma spp. morphology
install.packages("xlsx");install.packages("clustvarsel");install.packages("mclust")
install.packages("ggplot2")
install.packages("rJava")
library("rJava")
library(xlsx); library(ggplot2); 
library(clustvarsel);
library(mclust);
library(gridExtra);
library(grid);
library(glmm)

# upload data
setwd("~/Dropbox/syma_speciation/")
syma <- read.csv("data/morphology.csv")
# data<-syma[syma$Sex=="m",] # filter by sex -- not super reliable
data <-syma[syma$include=="y",] #filter by NAs
symashade <- c("#228B22", "#FFFFA8", "#C6E01B") # for quick plots
measurements <- data[,14:19]
species <- as.character(data$three_species_model)
measurements <- cbind(measurements,species)
# data_2 <- data[!data$three_species_model=="ochracea",] 
two_species <- data$english_name
levels(two_species) <- c("megarhyncha","torotoro")
three_species <- data$three_species_model

# subset by species, t-tests (all sig even w/ bonferroni correction)
torotoro <- measurements[measurements$species=="torotoro",]
megarhyncha <- measurements[measurements$species=="megarhyncha",]
t.test(torotoro$bill_from_nostril, megarhyncha$bill_from_nostril) #sig
t.test(torotoro$bill_width, megarhyncha$bill_width) #sig
t.test(torotoro$bill_depth, megarhyncha$bill_depth) #sig
t.test(torotoro$tarsus, megarhyncha$tarsus) #sig
t.test(torotoro$wing_chord, megarhyncha$wing_chord) #sig
t.test(torotoro$tail_length, megarhyncha$tail_length) #sig

# glmms - work in progress
# mod1 <- glmm(PC1 ~ 0 + english_name, data = data, random

# perform PCA, export df
pca<-prcomp(data[,14:19])
pca$rotation # PC loadings
summary(pca) # variance by axis
pca$x  # scores for each individual on each PC axis
data<-cbind(pca$x, data) # combining the PC scores with the old dataframe 
useful <- clustvarsel(pca$x, verbose = TRUE) # test to determine which PCs are most useful for group discrimination in NMMs
write.csv(data, "data/syma_spp_morphology.csv")

# evaluate the three species model w/ clustering, recommended PC scores
PCs <- data[,1:2]
table(two_species)
clPairs(PCs, three_species, color=symashade)
BIC2 <- mclustBIC(PCs)
plot(BIC2)
summary(BIC2)
mod1pc <- Mclust(PCs, x = BIC2) #2
summary(mod1pc, parameters = TRUE) #log.likelihood = -378.0914, df = 9, BIC = -794.025
plot(mod1pc, what = "classification",col=symashade)
legend(4, -6, c("megarhyncha","torotoro"),fill=symashade)
table(three_species, mod2pc$classification)
1-(9/67) # = 86.6% accuracy in assigning to two spp. clusters