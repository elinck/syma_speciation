#d-test summaries
install.packages("ggridges", repo = "http://ftp.osuosl.org/pub/cran/"); 
library(pbapply);library(data.table);library(ggridges);library(plyr);library(ggplot2);library(viridis)
install.packages("pbapply")
setwd("~/Dropbox/syma_speciation/")
d <- fread("raw_data/out.txt",data.table = F)

# rename by geography specific subspecies
d$H1[grepl("13_toro|32_toro|8_toro",d$H1)] <- "toro_flav"
d$H2[grepl("13_toro|32_toro|8_toro",d$H2)] <- "toro_flav"
d$H3[grepl("13_toro|32_toro|8_toro",d$H3)] <- "toro_flav"
d$H1[grepl("11_toro|9_toro",d$H1)] <- "toro_meek"
d$H2[grepl("11_toro|9_toro",d$H2)] <- "toro_meek"
d$H3[grepl("11_toro|9_toro",d$H3)] <- "toro_meek"
d$H1[grepl("39_toro|21_toro",d$H1)] <- "toro_toro"
d$H2[grepl("39_toro|21_toro",d$H2)] <- "toro_toro"
d$H3[grepl("39_toro|21_toro",d$H3)] <- "toro_toro"
d$H1[grepl("10_toro",d$H1)] <- "toro_tent"
d$H2[grepl("10_toro",d$H2)] <- "toro_tent"
d$H3[grepl("10_toro",d$H3)] <- "toro_tent"
d$H1[grepl("ochr",d$H1)] <- "ochr"
d$H2[grepl("ochr",d$H2)] <- "ochr"
d$H3[grepl("ochr",d$H3)] <- "ochr"
d$H1[grepl("1_mega|24_mega|40_mega",d$H1)] <- "mega_sella"
d$H2[grepl("1_mega|24_mega|40_mega",d$H2)] <- "mega_sella"
d$H3[grepl("1_mega|24_mega|40_mega",d$H3)] <- "mega_sella"
d$H1[grepl("18_mega|19_mega|20_mega|23_mega",d$H1)] <- "mega_mega"
d$H2[grepl("18_mega|19_mega|20_mega|23_mega",d$H2)] <- "mega_mega"
d$H3[grepl("18_mega|19_mega|20_mega|23_mega",d$H3)] <- "mega_mega"
d$H1[grepl("4_mega|6_mega|27_mega",d$H1)] <- "mega_wellsi"
d$H2[grepl("4_mega|6_mega|27_mega",d$H2)] <- "mega_wellsi"
d$H3[grepl("4_mega|6_mega|27_mega",d$H3)] <- "mega_wellsi"

# sort trees
trees <- c()
trees <- apply(d,1,function(e) {
  if(grepl("toro",e[1])){
    sp1 <- "toro"
  } else if(grepl("ochr",e[1])){
    sp1 <- "ochr"
  } else if(grepl("mega_sella",e[1])){
    sp1 <- "mega_sella"
  } else if(grepl("mega_mega",e[1])){
    sp1 <- "mega_mega"
  } else if(grepl("mega_wellsi",e[1])){
    sp1 <- "mega_wellsi"
  }
  if(grepl("toro",e[2])){
    sp2 <- "toro"
  } else if(grepl("ochr",e[2])){
    sp2 <- "ochr"
  } else if(grepl("mega_sella",e[2])){
    sp2 <- "mega_sella"
  } else if(grepl("mega_mega",e[2])){
    sp2 <- "mega_mega"
  } else if(grepl("mega_wellsi",e[2])){
    sp2 <- "mega_wellsi"
  }
  if(grepl("toro",e[3])){
    sp3 <- "toro"
  } else if(grepl("ochr",e[3])){
    sp3 <- "ochr"
  } else if(grepl("mega_sella",e[3])){
    sp3 <- "mega_sella"
  } else if(grepl("mega_mega",e[3])){
    sp3 <- "mega_mega"
  } else if(grepl("mega_wellsi",e[3])){
    sp3 <- "mega_wellsi"
  }
  paste0("(",sp1,",",sp2,"),",sp3,")")
})
d <- d[,-10]
d$tree <- trees

d1 <- subset(d,tree %in% c("(mega_mega,mega_mega),toro)","(mega_sella,mega_sella),toro)",
                          "(mega_wellsi,mega_wellsi),toro)","(mega_mega,mega_sella),toro)",
                          "(mega_mega,mega_wellsi),toro)","(mega_sella,mega_mega),toro)",
                          "(mega_wellsi,mega_mega),toro)","(toro,toro),mega_mega)",
                          "(toro,toro),mega_wellsi)","(toro,toro),mega_sella)","(toro,toro),ochr)",
                          "(mega_mega,ochr),toro)","(mega_sella,ochr),toro)","(mega_wellsi,ochr),toro)",
                          "(ochr,mega_mega),toro)","(ochr,mega_sella),toro)","(ochr,mega_wellsi),toro)",
                          "(mega_mega,mega_mega),ochr)", "(mega_mega,mega_wellsi),ochr)", "(mega_mega,mega_sella),ochr)",
                          "(mega_sella,mega_sella),ochr)","(mega_sella,mega_mega),ochr)", "(mega_sella,mega_wellsi),ochr)",
                          "(mega_wellsi,mega_wellsi),ochr)","(mega_sella,mega_wellsi),ochr)",
                          "(mega_wellsi,mega_sella),ochr)","(mega_wellsi,mega_mega),ochr)"))

d1$tree <- as.factor(d1$tree)
levels(d1$tree)
levels(d1$tree)[which(levels(d1$tree)=="(toro,toro),ochr)")] <- "torotoro & ochracea"
levels(d1$tree)[which(levels(d1$tree)=="(toro,toro),mega_wellsi)")] <- "torotoro & megarhyncha W"
levels(d1$tree)[which(levels(d1$tree)=="(toro,toro),mega_sella)")] <- "torotoro & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(toro,toro),mega_mega)")] <- "torotoro & megarhyncha C"
levels(d1$tree)[which(levels(d1$tree)=="(ochr,mega_wellsi),toro)")] <- "torotoro & megarhyncha W"
levels(d1$tree)[which(levels(d1$tree)=="(ochr,mega_sella),toro)")] <- "torotoro & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_wellsi,ochr),toro)")] <- "torotoro & megarhyncha W"
levels(d1$tree)[which(levels(d1$tree)=="(mega_wellsi,mega_wellsi),toro)")] <- "torotoro & megarhyncha W"
levels(d1$tree)[which(levels(d1$tree)=="(mega_wellsi,mega_wellsi),ochr)")] <- "ochracea & megarhyncha W"
levels(d1$tree)[which(levels(d1$tree)=="(mega_wellsi,mega_sella),ochr)")] <- "ochracea & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_sella,ochr),toro)")] <- "torotoro & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_sella,mega_wellsi),ochr)")] <- "ochracea & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_sella,mega_sella),toro)")] <- "torotoro & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_sella,mega_sella),ochr)")] <- "ochracea & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_sella,mega_mega),toro)")] <- "torotoro & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_sella,mega_mega),ochr)")] <- "ochracea & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_mega,ochr),toro)")] <- "torotoro & megarhyncha C"
levels(d1$tree)[which(levels(d1$tree)=="(mega_mega,mega_wellsi),toro)")] <- "torotoro & megarhyncha C"
levels(d1$tree)[which(levels(d1$tree)=="(mega_mega,mega_wellsi),ochr)")] <- "ochracea & megarhyncha C"
levels(d1$tree)[which(levels(d1$tree)=="(mega_mega,mega_sella),toro)")] <- "torotoro & megarhyncha C"
levels(d1$tree)[which(levels(d1$tree)=="(mega_mega,mega_sella),ochr)")] <- "ochracea & megarhyncha H"
levels(d1$tree)[which(levels(d1$tree)=="(mega_mega,mega_mega),toro)")] <- "torotoro & megarhyncha C"
levels(d1$tree)[which(levels(d1$tree)=="(mega_mega,mega_mega),ochr)")] <- "ochracrea & megarhyncha C"
d1$Z <- abs(d1$Z)

write.csv(d1, "data/d_stat_results.csv")
d1 <- read.csv("~/Dropbox/syma_speciation/data/d_stat_results.csv")

# figure out zscore significance threshold for bonferroni corrected p value of 0.001
p <- 0.005/31
qnorm(p) #-3.59
d1$prop_sig <- signif(length(d1$Z[abs(d1$Z)>3.59])/length(d1$Z),3)
d1$Z <- abs(d1$Z)

# make summary table

dtable <- ddply(d1,.(tree),summarize,
                nABBA=round(median(nABBA)),
                nBABA=round(median(nBABA)),
                median_D=format(signif(median(Dstat),3),scientific=T),
                #D_CI=paste(format(signif(quantile(abs(Dstat),0.025),3),scientific=T),
                #          format(signif(quantile(abs(Dstat),0.975),3),scientific=T),sep=" - "),
                median_Z=format(signif(abs(median(Z)),3)),
                #Z_CI=paste(signif(quantile(abs(Z),0.025),3),signif(quantile(abs(Z),0.975),3),sep=" - "),
                prop_sig=signif(length(Z[abs(Z)>3.59])/length(Z),3)
)
dtable$median_Z <- as.numeric(dtable$median_Z)

# write summary for plotting
write.csv(dtable, file="data/d_summary.csv")

# test visualize
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
mag <- map2color(dtable$prop_sig,magma(31))
ggplot(data=d1,aes(x=abs(Z),y=tree,fill=tree))+
  theme_bw() +
  theme(axis.title=element_text(size=8),
        axis.text=element_text(size=8))+
  scale_fill_manual(values=mag,guide=FALSE)+
  ylab("Test")+xlab("abs(Z)")+
  xlim(-5,100)+
  #geom_vline(aes(xintercept=4.3),col="red",lwd=0.5)+
  geom_vline(aes(xintercept=3.59),col="black",linetype="dashed",lwd=0.5)+
  geom_density_ridges(scale=1.25,lwd=0.35)+
  geom_point(data=dtable,shape=8,size=.75,aes(x=100,y=as.integer(factor(tree))+.3,col=median_Z>3.59))+
  scale_color_manual(values=rep(c("white","black"),nrow(dtable)),guide=F)


