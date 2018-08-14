### analysis of acoustic divergence between S. torotoro and S. megarhyncha

library(warbleR);library(ggfortify);library(ggplot2);
library(cluster);library(plyr);library(lmerTest)

setwd("~/Dropbox/syma_speciation")

# download xeno canto material 
# symaxc <- querxc(qword = "Syma", download = FALSE) 
# View(symaxc)
# symasong <- symaxc[grep("song", symaxc$Vocalization_type, ignore.case = TRUE),]
# querxc(X = symasong) 
# write.csv(symasong, "symasongxc.csv", row.names = FALSE)
# mp32wav() # convert mp3s to wavs
# checkwavs() # all files are readable
# wavs <- list.files(pattern="wav$") # create list of wavs
# sub <- wavs[c(27)] # subset for testing 
# lspec(flist = sub, flim = c(0,3), sxrow = 10, rows = 15, ovlp = 10, it = "tiff") # test parameters
# lspec(flim = c(1.5, 11), sxrow = 15, rows = 20, ovlp = 10, it = "tiff") # produce spectrograms

# manualoc(wl = 512, flim = c(0,4), seltime = 1, tdisp = NULL, reccomm =
#           FALSE, wn = "hanning", title = TRUE, selcomm = FALSE, osci = FALSE, player =
#           NULL, pal = reverse.gray.colors.2, path = NULL, flist = NULL) #manual detection params

# syma.ad <- autodetec(flist = wavs, bp = c(1, 3), threshold = 4, mindur = 0.5, maxdur = 15, envt="abs",
#                     ssmooth = 1500, ls = TRUE, res = 100, 
#                     flim = c(1, 3), wl = 512, set =TRUE, sxrow = 15, rows = 20, 
#                     redo = FALSE, it = "tiff") #best parameters 040417
# write.csv(syma.ad, file = "syma_selections_all.csv")
# table(syma.ad$sound.files) # view selections by individuaL 
# set.seed(5) #set seed
# X <- syma.ad[sample(1:nrow(syma.ad),(nrow(syma.ad)*0.1)), ] #test on 10% of selections
# snrspecs(X = syma.ad, flim = c(1, 4), snrmar = 0.2, it = "tiff") # print noise margins

# syma.ad <- read.csv("data/syma_selections.csv")
# syma.snr <- sig2noise(X = syma.ad, mar = 0.2) # adjust noise to signal ratio 

# syma.hisnr <- syma.snr[ave(-syma.snr$SNR, syma.snr$sound.files, FUN = rank) <= 5,] # select the 5 most high quality songs per individual 
# table(syma.hisnr$sound.files) # view no. of high quality selections per individual
# write.csv(syma.hisnr, file = "syma_hiquality.csv") # save as CSV
# syma.hisnr <- read.csv("syma_hiquality.csv")
# params <- specan(syma.hisnr, bp = c(1,4), threshold = 4)# calculate acoustic params

# read in params csv, perform PCA
params <- read.csv("data/params.csv")
pca <- prcomp(x = params[, sapply(params, is.numeric)], scale. = TRUE) #perform PCA
summary(pca) #summarize PCA scores
scores <- pca$x[,1:3] # select first 3 PCs
params <- cbind(scores,params)
names(params)[names(params) == 'sp.'] <- 'two_species'
write.csv(params,"data/syma_spp_calls.csv")

# glmm -- species as fixed effect, individual as random effect
mod1 <- lmerTest::lmer(PC1 ~ species + (1 | indiv), data = df)
summary(mod1) # sig 9.43e-06 ***
