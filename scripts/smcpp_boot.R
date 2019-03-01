
install.packages("plyr", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("data.table", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("ggplot2", repo = "http://ftp.osuosl.org/pub/cran/");
library(plyr);library(magrittr);library(data.table);library(ggplot2)

contigs <- read.table("/media/burke/bigMac/ethan/keep_contigs.txt", header = FALSE, row.names = NULL)
contigs <- contigs$V1

#megarhyncha bootstraps
files <- list.files("/media/burke/bigMac/ethan/smcpp/mega/")
megaboot <- files[grepl("mega|EL40_toro",files)]
for(i in 1:10){
  boot <- sample(contigs,length(contigs),replace=T)
  set <- grep(paste(boot,collapse = "|"),megaboot,value=T)
  set <- paste0("/media/burke/bigMac/ethan/smcpp/mega/",set)
  print(set)
  write.table(set,paste0("/media/burke/bigMac/ethan/smcpp/mega_boot/",i,".txt"),row.names=F,col.names=F,quote=F)
}

#torotoro bootstraps
files <- list.files("/media/burke/bigMac/ethan/smcpp/toro/")
toroboot <- files[grepl("EL10|EL11|EL13|EL21|EL32|EL39|EL8|EL9",files)]
for(i in 1:10){
  boot <- sample(contigs,length(contigs),replace=T)
  set <- grep(paste(boot,collapse = "|"),toroboot,value=T)
  set <- paste0("/media/burke/bigMac/ethan/smcpp/toro/",set)
  print(set)
  write.table(set,paste0("/media/burke/bigMac/ethan/smcpp/toro_boot/",i,".txt"),row.names=F,col.names=F,quote=F)
}
