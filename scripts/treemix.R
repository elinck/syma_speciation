#convert vcftools 012 output to treemix
library(pbapply);library(magrittr);library(data.table)
input <- "~/Dropbox/syma_speciation/raw_data/syma.treemix.012"
samples <- read.table("~/Dropbox/syma_speciation/raw_data/syma.treemix.012.indv",header=F)
data <- fread(input,data.table=F)[,-1]
rownames(data) <- samples[,1]
drop <- c(26:34)
data <- data[-drop,]

#samples <- subset(samples, samples$V1 %in% subsample)
pops <- list(toro_n=c(1,8,13:14,16:17,19:22,24,27),
            toro_e=c(2:3,18,28),
            mega_h=c(6,10,15),
            mega_c=c(4:5,7,9,11,25:26),
            ochr=c(12,23))
data <- pbapply(data, 2, function(x){ # replace all -1 missing data w/ 0
  ifelse(x < 0, 0, x)})

tmix <- pbapply(data,2,function(e){
  toro_n <- e[pops$toro_n] %>% sum()
  toro_e <- e[pops$toro_e] %>% sum()
  mega_h <- e[pops$mega_h] %>% sum()
  mega_c <- e[pops$mega_c] %>% sum()
  ochr <- e[pops$ochr] %>% sum()
  out <- c(paste(28-toro_n,toro_n,sep=","),
           paste(28-toro_e,toro_e,sep=","),
           paste(28-mega_h,mega_h,sep=","),
           paste(28-mega_c,mega_c,sep=","),
           paste(28-ochr,ochr,sep=","))
})
tmix <- data.frame(t(tmix))
colnames(tmix) <- c("toro_n","toro_e","mega_h","mega_c","ochr")

write.table(tmix,"raw_data/syma.tmix",
            row.names = F,col.names = T,quote = F,sep=" ")
system("gzip raw_data/syma.tmix")