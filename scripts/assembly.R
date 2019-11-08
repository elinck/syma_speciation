### truseq resequencing SNP calling pipeline for degraded DNA                
### modified by E. Linck, wrappers inspired by C.J. Battey

### load libraries
setwd("/media/burke/bigMac/ethan/")
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/"); 
install.packages("rebus", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("inline", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("gam", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("Rcpp", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("ggplot2", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("RcppGSL", repo = "http://ftp.osuosl.org/pub/cran/")
library(magrittr)

### mask repeats from reference using bbmask entropy approach
system("bbmask.sh in=B10K-DU-024-03.genomic.fa out=masked.ref.fa entropy=0.7")

### run bbduk to trim adapters
R1 <- list.files("/media/burke/bigMac/ethan/syma_raw",full.names=T) %>% grep("R1",.,value=T)
R2 <- list.files("/media/burke/bigMac/ethan/syma_raw",full.names=T) %>% grep("R2",.,value=T)
commands <- c()
for(i in 1:length(R1)){
  commands[i] <- paste0("bbduk.sh in1=", R1[i],
                        " in2=",R2[i],
                        " out1=", R1[i] %>% basename() %>% strsplit(".fastq") %>% unlist() %>% .[1],".trimmed.fq",
                        " out2=", R2[i] %>% basename() %>% strsplit(".fastq") %>% unlist() %>% .[1],".trimmed.fq",
                        " ref=/media/burke/bigMac/ethan/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo t=20")
}
for(i in commands){
  system(i)
}
system("mkdir syma_trimmed; mv *trimmed* syma_trimmed")

### run bbduk to quality trim
R1 <- list.files("/media/burke/bigMac/ethan/syma_raw/syma_trimmed",full.names=T) %>% grep("R1",.,value=T)
R2 <- list.files("/media/burke/bigMac/ethan/syma_raw/syma_trimmed",full.names=T) %>% grep("R2",.,value=T)
commands <- c()
for(i in 1:length(R1)){
  commands[i] <- paste0("bbduk.sh in1=", R1[i],
                        " in2=",R2[i],
                        " out1=", R1[i] %>% basename() %>% strsplit(".trimmed.fq") %>% unlist() %>% .[1],".clean.fq",
                        " out2=", R2[i] %>% basename() %>% strsplit(".trimmed.fq") %>% unlist() %>% .[1],".clean.fq",
                        " qtrim=r trimq=10 ftr=100")
}
for(i in 1:length(R1)){
  system(commands[i])
}

system("mkdir syma_clean; mv *clean* syma_clean")

### align trimmed reads to reference sample
R1 <- list.files("/media/burke/bigMac/ethan/syma_cleaned",full.names=T) %>% grep("R1",.,value=T)
R2 <- list.files("/media/burke/bigMac/ethan/syma_cleaned",full.names=T) %>% grep("R2",.,value=T)
commands <- c()
sampleID <- c()
for(i in 1:length(R1)){
  commands[i] <- paste0("bbmap.sh", 
                        " in1=", R1[i], 
                        " in2=", R2[i],
                        " out=", R1[i] %>% basename() %>% strsplit(".clean.fq") %>% unlist() %>% .[1], ".sam",
                        " ref=/media/burke/bigMac/ethan/masked.ref.fa", 
                        " t=30 slow k=10 maxindel=200 minratio=0.1 bamscript=bs.sh >> alignment_log.txt 2>&1; sh bs.sh")
}
for(i in 1:length(R1)){
  system(commands[i])
}

### fuck I guess I have to use GATK after all
setwd("/media/burke/bigMac/ethan/alignment/sorted/") 
# can't figure out a reasonable regex for this because I suck at coding
bam_wgs <- c("EL1_mega_R1_sorted.bam", "EL10_toro_R1_sorted.bam", "EL11_toro_R1_sorted.bam", "EL13_toro_R1_sorted.bam", "EL18_mega_R1_sorted.bam",
             "EL19_mega_R1_sorted.bam", "EL20_mega_R1_sorted.bam", "EL21_toro_R1_sorted.bam", "EL23_mega_R1_sorted.bam", "EL24_mega_R1_sorted.bam",
             "EL27_mega_R1_sorted.bam", "EL29_ochr_R1_sorted.bam", "EL32_toro_R1_sorted.bam", "EL39_toro_R1_sorted.bam", "EL4_mega_R1_sorted.bam",
             "EL40_toro_R1_sorted.bam", "EL6_mega_R1_sorted.bam", "EL8_toro_R1_sorted.bam", "EL9_toro_R1_sorted.bam", "EL5_ochr_R1_sorted.bam")
bam_hyrad <- c("EL45_toro_R1_sorted.bam", "EL46_toro_R1_sorted.bam", "EL47_toro_R1_sorted.bam", "EL48_ochr_R1_sorted.bam", "EL49_toro_R1_sorted.bam",
               "EL41_toro_R1_sorted.bam",  "EL50_toro_R1_sorted.bam", "EL51_toro_R1_sorted.bam", "EL52_toro_R1_sorted.bam", "EL53_toro_R1_sorted.bam",
               "EL54_toro_R1_sorted.bam", "EL55_toro_R1_sorted.bam", "EL56_toro_R1_sorted.bam", "EL57_toro_R1_sorted.bam", "EL58_toro_R1_sorted.bam",
               "EL59_toro_R1_sorted.bam", "EL42_toro_R1_sorted.bam",  "EL60_toro_R1_sorted.bam", "EL43_toro_R1_sorted.bam", "EL44_toro_R1_sorted.bam")
commands <- c()
for(i in 1:length(bam_wgs)){
  sampleID <- basename(bam_wgs[i]) %>% strsplit("_R1_sorted.bam") %>% unlist() %>% .[1]
  commands[i] <- paste0("java -Xmx10g -jar /media/burke/bigMac/Dropbox/tools/picard.jar AddOrReplaceReadGroups",
                        " I=",bam_wgs[i],
                        " O=tmp.bam",
                        " R=/media/burke/bigMac/ethan/masked.ref.fa",
                        " RGID=1",
                        " RGLB=syma_wgs",
                        " RGPL=illumina",
                        " RGPU=pool1",
                        " RGSM=",sampleID," >> rg_log.txt 2>&1;",
                        " mv tmp.bam ",bam_wgs[i],";",
                        " samtools index ",bam_wgs[i]
  )
}
for(i in 1:length(bam_wgs)){
  system(commands[i])
}
for(i in 1:length(bam_hyrad)){
  sampleID <- basename(bam_hyrad[i]) %>% strsplit("_R1_sorted.bam") %>% unlist() %>% .[1]
  commands[i] <- paste0("java -Xmx10g -jar /media/burke/bigMac/Dropbox/tools/picard.jar AddOrReplaceReadGroups",
                        " I=",bam_hyrad[i],
                        " O=tmp.bam",
                        " R=/media/burke/bigMac/ethan/masked.ref.fa",
                        " RGID=2",
                        " RGLB=syma_hyrad",
                        " RGPL=illumina",
                        " RGPU=pool2",
                        " RGSM=",sampleID," >> rg_log.txt 2>&1;",
                        
                        "mv tmp.bam ",bam_hyrad[i],";",
                        
                        "samtools index ",bam_hyrad[i]
  )
}
for(i in 1:length(bam_hyrad)){
  system(commands[i])
}

### local realignment
bams <- list.files("/media/burke/bigMac/ethan/alignment", full.names=T) %>% grep(".bam",.,value=T)
for(i in 1:length(bams)){
  sampleID[i] <- basename(bams[i]) %>% strsplit("_R1_sorted.bam") %>% unlist() %>% .[1]
  commands[i] <- paste0("java -Xmx40g -jar /media/burke/bigMac/Dropbox/tools/GenomeAnalysisTK.jar",
                        " -T IndelRealigner",
                        " -R /media/burke/bigMac/ethan/masked.ref.fa",
                        " -I ",bams[i],
                        " -targetIntervals forIndelRealigner.intervals",
                        " -o ",sampleID[i],"_realigned.bam"
                        )
}
for(i in 1:length(bams)){
  system(commands[i])
}

### run mapDamage to recalibrate Q scores
bams <- list.files("/media/burke/bigMac/ethan/alignment/working_bams", full.names=T) %>% grep("realigned.bam$",.,value=T)
for(i in 1:length(bams)){
  sampleID[i] <- basename(bams[i]) %>% strsplit("_realigned.bam") %>% unlist() %>% .[1]
  commands[i] <- paste0("mapDamage -i ",
                        bams[i]," -r /media/burke/bigMac/ethan/masked.ref.fa --merge-reference-sequences --rescale"
  )
}
for(i in 1:length(bams)){
  system(commands[i])
}
