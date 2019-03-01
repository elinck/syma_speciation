# *Syma* bioinformatics pipeline and analysis

This digital notebook documents commands run to analyse hyRAD and whole genome sequence data from *Syma torotoro* and
*S. megarhyncha.* Downstream genotypic data manipulation, analysis of phenotypic data, and plotting in `scripts/`.
Paths and descriptions of scripts are found in `README.md.`

## Sequence alignment and variant calling pipeline

Requires `R`, [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/), [GATK](https://software.broadinstitute.org/gatk/), [Picard tools](https://broadinstitute.github.io/picard/),
[mapDamage](https://ginolhac.github.io/mapDamage/), and [MITObim](https://github.com/chrishah/MITObim). 

Full commands available in `assembly.R`. Probably best to tweak and run in discrete chunks to debug rather than via `Rscript.`

Mask repeats from reference genome:

```
bbmask.sh in=B10K-DU-024-03.genomic.fa out=masked.ref.fa entropy=0.7
```

Open R via the terminal, install and load `magrittr`, install packages for `MITObim.`

```R
setwd("/media/burke/bigMac/ethan/")
library(magrittr)
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/"); 
install.packages("rebus", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("inline", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("gam", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("Rcpp", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("ggplot2", repo = "http://ftp.osuosl.org/pub/cran/");
install.packages("RcppGSL", repo = "http://ftp.osuosl.org/pub/cran/")
```

Run `bbduk` from `bbtools` to trim adapters, make new directory and move files.

```R
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
```

Run `bbduk` for quality trimming, make new directory and move files. 

```R
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

```

Use `bbmap` to align reads to repeat masked reference genome, using parameters for high sensitivity (e.g., low DNA quality and
a divergent reference genome).

```R
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
```

Add read group headers for `GATK.` Please write me if you have a regex that works for this, lol. 

```R
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
```

Perform local realignment around indels, etc.

```R
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
```

Recalibrate quality scores for postmortem DNA damage. Bayesian, takes a while!

```R
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
```

Assembly mtDNA genomes. Warning: looping here is mysteriously memory intensive; might be better to generate commands and enter
manually. Alternatively, modify the list of commands in `mtDNA.sh`.

```R
paths <- list.files("/media/burke/bigMac/ethan/syma_raw",full.names=T) %>% grep("R1",.,value=T)
names <- list.files("/media/burke/bigMac/ethan/syma_raw",full.names=F) %>% grep("R1",.,value=T)
commands <- c()
sampleID <- c()
for(i in 1:length(names)){
  sampleID[i] <- names[i] %>% strsplit("_R1.fastq.gz") %>% unlist() %>% .[1]
  commands[i] <- paste0("mkdir ", sampleID[i], ";",
                        " cd ", sampleID[i], ";",
                        " /media/burke/bigMac/ethan/mtDNA/MITObim/MITObim.pl -start 1 -end 5 -sample ", sampleID[i],
                        " -ref syma_mtdna", i,
                        " -readpool ", paths[i],
                        " -quick /media/burke/bigMac/ethan/mtDNA/t_sanctus_ref.fasta &> log;")
}

for(i in 1:length(paths)){
  system(commands[i])
}
```

## Variant calling and VCF filtering. 

Requires [GATK](https://software.broadinstitute.org/gatk/) and [VCFtools](https://vcftools.github.io/man_latest.html). 

First, we create two .vcf files for analyses: a file with only variant sites, and a file including invariant sites for accurate calculation of windowed
summary statistics. (This could all be done from the "full" .vcf through filtering, but I ended up creating both due to the order of analyses.) A simple
shell script calls the length `GATK` command for variant sites across all samples: 

```
bash unified_genotype_variant.sh
```

This should take roughly 24 hours. To call invariant sites as well is closer to 5 days:

```
bash unified_genotype_all.sh
```

After completion, we filter the raw variant only .vcf for quality and missing data across all samples. We begin by dropping sites with more than 
0.20 missing data, and those with only a single copy of the minor allele. 

```
vcftools --vcf syma_recal_gatk.vcf --max-missing 0.80 --mac 2 --minQ 30 --recode --recode-INFO-all --out syma.80.mac2
```

(After filtering, we kept 384904 out of a possible 78882912 sites.) 


We next want to establish a maximum read depth for the entire alignment. 

```
vcftools --vcf syma.80.mac2.recode.vcf --site-depth --out depth.dist
cut -f3 depth.dist.ldepth > depth.dist
mawk '!/D/' depth.dist | mawk -v x=39 '{print $1/x}' > meandepthpersite
```

Plot as a histogram in the terminal window (rad as hell; via [dDocent](http://ddocent.com/filtering/))

```
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [0:150] 
unset label
set title "Histogram of mean depth per site"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics 5
plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF  
```  

Evaluating these data, we next drop sites covered with fewer than 3 reads and more than 120:  

```
vcftools --vcf syma.80.mac2.recode.vcf --minDP 3 --maxDP 120 --recode --recode-INFO-all --out syma.80.mac2.dp3
```

(After filtering, kept 379297 out of a possible 379297 sites.) 

Next we determine which individuals have >0.5 missing data and remove them (3, in this case) from the .vcf.

```
vcftools --vcf syma.80.mac2.dp3.recode.vcf --missing-indv
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
vcftools --vcf syma.80.mac2.dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out syma.80.mac2.dp3.rmmd
```

To generate the .vcf used in PCA and the neighbor-joining tree analyses, drop sites with more than 5% missing data, with a minor allele frequency less than 0.05:

```
vcftools --vcf syma.80.mac2.dp3.rmmd.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out syma.80.mac2.dp3.rmmd.m20 
vcftools --vcf syma.80.mac2.dp3.rmmd.recode.vcf  --recode --recode-INFO-all --out syma.80.mac2.dp3.rmmd.m20 
```

(After filtering, we kept 23634 out of a possible 379297 sites.)

To create the .vcf for `SMC++` analyses, we return to the original .vcf, drop sites with more than 25% missing data, fewer than 3 reads, more than 120, and a minimum Q score of 30:

```
vcftools --vcf syma_recal_gatk.vcf --remove hyrad.bams --max-missing 0.75 --max-alleles 2 --minDP 3 --maxDP 120 --minQ 30 --recode --recode-INFO-all --out syma.smcpp.revised

```

We'll also zip the file:

```
bgzip -c syma.smcpp.revised.recode.vcf > syma.smcpp.revised.recode.vcf.gz
tabix /media/burke/bigMac/ethan/alignment/recalibrated/syma.smcpp.revised.recode.vcf.gz
```

(After filtering, kept 1694771 out of a possible 78882912 sites).


For ABBA/BABA tests we'll do something similar, but for all WGS samples:
```
vcftools --vcf /media/burke/bigMac/ethan/alignment/recalibrated/syma_recal_gatk.vcf --remove /media/burke/bigMac/ethan/alignment/recalibrated/hyrad.bams --max-missing 0.75 --min-alleles 2 --max-alleles 2 --minDP 3 --maxDP 120 --minQ 30 --recode --recode-INFO-all --out syma.wgs.abba

```

(After filtering, kept 1385081 out of a possible 78882912 sites.)


To create the .vcf file later formatted for moments (see below), we'll take the .vcf for SMC++ and thin it for LD: 

```
vcftools --vcf syma.smcpp.revised.recode.vcf --thin 50000 --recode --recode-INFO-all --out syma.moments.revised.unlinked 

```

(After filtering, kept 19041 out of a possible 1694771 sites.) 


For RAiSD, we'll take the SMC++ file and create species-specific .vcfs, first for *megarhyncha*...

```
vcftools --vcf syma.smcpp.revised.recode.vcf --remove remove.toro --recode --recode-INFO-all --out mega.raisd
```

...and then for *torotoro*: 

```
vcftools --vcf syma.smcpp.revised.recode.vcf --remove remove.mega --recode --recode-INFO-all --out toro.raisd
```

For sliding window summary statistics, we'll call invariant sites along with variants for our largest (huge) `.vcf`, with a minQ score of 30 and a minDP of 3:

```
vcftools --gzvcf syma_recal_gatk_all.vcf.gz --remove hyrad.bams --minQ 30 --minDP 3  --recode --recode-INFO-all --out syma.wgs.all 
```

(After filtering, we kept 102606487 out of a possible 989764264 sites.) 

Finally, we'll extract Chromosome 5 for comparison w/ simulated data, then restrict our analysis to *megarhyncha*

```
vcftools --gzvcf syma.wgs.all.recode.vcf.gz --chr scaffold_13 --recode --recode-INFO-all --out syma.scaf13
vcftools --vcf syma.scaf13.recode.vcf --remove remove.toro --recode --recode-INFO-all --out syma.scaf13.mega
```

(After filtering, kept 758116 out of a possible 102606487 sites.)

## Morphology and bioacoustic analysis

We used the `R` scripts `morphology_analysis.R` and `vocalization_analysis.R` to analyse morphological data (`data/morphology.csv`)
and bioacoustic data (`data/params.csv`), respectively.

## Genotype analysis

We used the `R` script `genotype_analysis` to process `.vcf` data and prepare datasets for many of the analyses below.

## ADMIXTURE

Requires [ADMIXTURE](http://software.genetics.ucla.edu/admixture/), [PLINK](https://www.cog-genomics.org/plink2/).  
  
First, we reformat the LD-thinned SNP matrix to an "ordinary" 12 coded `PLINK` .ped file:

```
plink --vcf syma.str.recode.vcf --recode12 --out syma.admix --allow-extra-chr
```

We then run the output (`syma.admix.ped`) through ADMIXTURE for K=2 through K=5, visually inspecting the Q-matrix for population assignments: 

```
/media/burke/bigMac/ethan/alignment/recalibrated/admixture_linux-1.3.0/admixture syma.admix.ped 2 --cv
```
Loglikelihood: -106788.290035, CV error (K=2): 0.68220

```
/media/burke/bigMac/ethan/alignment/recalibrated/admixture_linux-1.3.0/admixture syma.admix.ped 3 --cv
```
Loglikelihood: -98148.979795, CV error (K=3): 0.79501

```
/media/burke/bigMac/ethan/alignment/recalibrated/admixture_linux-1.3.0/admixture syma.admix.ped 4 --cv
```
Loglikelihood: -93803.561686, CV error (K=4): 0.96999

```
/media/burke/bigMac/ethan/alignment/recalibrated/admixture_linux-1.3.0/admixture syma.admix.ped 5 --cv
```
Loglikelihood: -90818.879955, CV error (K=5): 1.11336


Across all runs, assignments are random with respect to geography, batch, and species identity, suggesting sample size is insufficient
given the level of discordance. We therefore move to supervised mode, using file `syma.admix.pop` in the same directory, with population
identities in a single column in the order of samples in the .vcf file:

```
/media/burke/bigMac/ethan/alignment/recalibrated/admixture_linux-1.3.0/admixture syma.admix.ped 3 --supervised --cv
```
Loglikelihood: -104066.811525, CV error (K=3): 0.59740

## ABBA-BABA

Requires [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD#Overview).

To calculate *D*-statistics (with "ABBA-BABA" tests), we use ANGSD's `-doAbbababa 1` [program](http://www.popgen.dk/angsd/index.php/Abbababa). We
first create a filelist with all WGS samples (`wgs.filelist`)

```
./angsd -out out -doAbbababa 1 -bam wgs.filelist -doCounts 1 -anc masked.ref.fa

Rscript R/jackKnife.R file=out.abbababa indNames=wgs.filelist outfile=out
```

We then process the outfile (`out.txt`) using the `R` script `abbababa.R`. 

## Moments

Requires [Moments](https://bitbucket.org/simongravel/moments).  
  
To perform demographic inference in moments, we'll take the `syma.moments.revised.txt` file output from `genotype_analysis.R` as input for scripts for
4 models: strict isolation (SI), isolation with migration (IM), secondary contact (SC), and initial isolation with migration (aka ancestral migration; IIM). 
   
We'll then run each model's script in `scripts/moments/*model.py`, import and process output in `genotype_analysis.R`, and then plot in `plotting.R`. We'll perform a 
likelihood ratio test between the top two models with `uncertainty.py.`

```
python SI_model.py;
python IM_model.py;
python SC_model.py;
python IIM_model.py;
uncertainty.py
```

## Treemix

Requires [Treemix](https://bitbucket.org/nygcresearch/treemix/wiki/Home).  
  
We can run Treemix by using the output from `vcftools --012` format, and an `R` script that comes with the program. 

```
vctools --vcf syma.str.recode.vcf --012 --out syma.treemix
Rscript treemix.R
treemix -i ~/Dropbox/syma_speciation/raw_data/syma.tmix.gz -o m1 -m 1 #ln(likelihood):-1289.4967 
treemix -i ~/Dropbox/syma_speciation/raw_data/syma.tmix.gz -o m2 -m 2 #ln(likelihood):-178.13221 
treemix -i ~/Dropbox/syma_speciation/raw_data/syma.tmix.gz -o m3 -m 3 #ln(likelihood):124.13581 
```

## SMC++

Requires [SMC++](https://github.com/popgenmethods/smcpp).  
  
To run `SMC++`, we'll first get a list of scaffolds over 1000000 bp in length

```
cut -f1-2 masked.ref.fa.fai masked.ref.fa.fai > scaffold_length.txt
awk '$2 >= 1000000' scaffold_length.txt > long_contigs.txt
cut -f1 long_contigs.txt > keep_contigs.txt

```

We can then run a bash script to sample over individuals and contigs and create .vcf files:

```
bash smcpp_data.sh
```

Next, bootstrap over these variables using the script `smcpp_boot.R`. (Edit for sample ID, paths, number of bootstraps)

```
Rscript smcpp_boot.R
```

We then run the program for each bootstrapped dataset:

```
bash smcpp_run.sh
```

We then plot the results as a quick-and-dirty `.pdf`, and export a dataframe for plotting 
```
cd /media/burke/bigMac/ethan/smcpp/models_mega
smc++ plot mega.plot.png 1.model.final.json 2.model.final.json 3.model.final.json 4.model.final.json 5.model.final.json 6.model.final.json 7.model.final.json 8.model.final.json 9.model.final.json 10.model.final.json -c
cd /media/burke/bigMac/ethan/smcpp/models_toro
smc++ plot toro.plot.png 1.model.final.json 2.model.final.json 3.model.final.json 4.model.final.json 5.model.final.json 6.model.final.json 7.model.final.json 8.model.final.json 9.model.final.json 10.model.final.json -c
```

## Windowed stats

Requires scripts from Simon Martin's (`genomics_general`)[https://github.com/simonhmartin/genomics_general] repository.  


First, we use the script `parseVCF.py` to convert our vcf (with invariant sites) to the `.geno` custom file format. 
```
python VCF_processing/parseVCF.py -i /media/burke/bigMac/ethan/alignment/recalibrated/syma.wgs.all.recode.vcf.gz --skipIndels | gzip > all.output.geno.gz
```

We then run `popgenWindows.py` with coordinate-based 50kb windows:
```
python popgenWindows.py -w 50000 --windType coordinate -g all.output.geno.gz -o syma_windows.csv.gz -f phased -T 24 -p meg EL19_mega,EL18_mega,EL20_mega,EL1_mega,EL23_mega,EL24_mega,EL27_mega,EL40_toro,EL4_mega,EL6_mega -p toro EL10_toro,EL11_toro,EL13_toro,EL21_toro,EL32_toro,EL39_toro,EL8_toro,EL9_toro;
```

## Selective sweeps

To determine whether divergence peaks have been shaped by positive selection, we compare a summary statistic
indicating the likelihood of a selective sweep from FST outlier regions to a random sample of nonoutlier regions. First,
we use species-specific vcf files and subset by a list of outlier and nonoutlier scaffolds generated in `genotype_analysis' 
from windowed summary statistic scans:

```
vcftools --vcf /media/burke/bigMac/ethan/alignment/recalibrated/mega.raisd.recode.vcf \
--chr scaffold_11 \
--chr scaffold_128 \
--chr scaffold_13 \
--chr scaffold_131 \
--chr scaffold_135 \
--chr scaffold_167 \
--chr scaffold_168 \
--chr scaffold_177 \
--chr scaffold_178 \
--chr scaffold_198 \
--chr scaffold_199 \
--chr scaffold_207 \
--chr scaffold_213 \
--chr scaffold_216 \
--chr scaffold_249 \
--chr scaffold_26 \
--chr scaffold_261 \
--chr scaffold_267 \
--chr scaffold_302 \
--chr scaffold_322 \
--chr scaffold_358 \
--chr scaffold_361 \
--chr scaffold_37 \
--chr scaffold_429 \
--chr scaffold_442 \
--chr scaffold_461 \
--chr scaffold_475 \
--chr scaffold_484 \
--chr scaffold_49 \
--chr scaffold_497 \
--chr scaffold_541 \
--chr scaffold_543 \
--chr scaffold_65 \
--recode --out mega.outliers

vcftools --vcf /media/burke/bigMac/ethan/alignment/recalibrated/toro.raisd.recode.vcf \
--chr scaffold_11 \
--chr scaffold_128 \
--chr scaffold_13 \
--chr scaffold_131 \
--chr scaffold_135 \
--chr scaffold_167 \
--chr scaffold_168 \
--chr scaffold_177 \
--chr scaffold_178 \
--chr scaffold_198 \
--chr scaffold_199 \
--chr scaffold_207 \
--chr scaffold_213 \
--chr scaffold_216 \
--chr scaffold_249 \
--chr scaffold_26 \
--chr scaffold_261 \
--chr scaffold_267 \
--chr scaffold_302 \
--chr scaffold_322 \
--chr scaffold_358 \
--chr scaffold_361 \
--chr scaffold_37 \
--chr scaffold_429 \
--chr scaffold_442 \
--chr scaffold_461 \
--chr scaffold_475 \
--chr scaffold_484 \
--chr scaffold_49 \
--chr scaffold_497 \
--chr scaffold_541 \
--chr scaffold_543 \
--chr scaffold_65 \
--recode --out toro.outliers

vcftools --vcf /media/burke/bigMac/ethan/alignment/recalibrated/mega.raisd.recode.vcf \
--chr scaffold_0 \
--chr scaffold_127 \
--chr scaffold_140 \
--chr scaffold_144 \
--chr scaffold_152 \
--chr scaffold_166 \
--chr scaffold_176 \
--chr scaffold_18 \
--chr scaffold_194 \
--chr scaffold_195 \
--chr scaffold_21 \
--chr scaffold_230 \
--chr scaffold_245 \
--chr scaffold_27 \
--chr scaffold_303 \
--chr scaffold_341 \
--chr scaffold_377 \
--chr scaffold_39 \
--chr scaffold_42 \
--chr scaffold_49 \
--chr scaffold_50 \
--chr scaffold_52 \
--chr scaffold_6 \
--chr scaffold_60 \
--chr scaffold_7 \
--chr scaffold_8 \
--chr scaffold_9 \
--chr scaffold_90 \
--chr scaffold_97 \
--chr scaffold_98 \
--recode --out mega.normal

vcftools --vcf /media/burke/bigMac/ethan/alignment/recalibrated/toro.raisd.recode.vcf \
--chr scaffold_0 \
--chr scaffold_127 \
--chr scaffold_140 \
--chr scaffold_144 \
--chr scaffold_152 \
--chr scaffold_166 \
--chr scaffold_176 \
--chr scaffold_18 \
--chr scaffold_194 \
--chr scaffold_195 \
--chr scaffold_21 \
--chr scaffold_230 \
--chr scaffold_245 \
--chr scaffold_27 \
--chr scaffold_303 \
--chr scaffold_341 \
--chr scaffold_377 \
--chr scaffold_39 \
--chr scaffold_42 \
--chr scaffold_49 \
--chr scaffold_50 \
--chr scaffold_52 \
--chr scaffold_6 \
--chr scaffold_60 \
--chr scaffold_7 \
--chr scaffold_8 \
--chr scaffold_9 \
--chr scaffold_90 \
--chr scaffold_97 \
--chr scaffold_98 \
--recode --out toro.normal
```

Then, we run `RAiSD` on each, outputting full reports. 

```
./RAiSD -n mega.outliers -I mega.outliers.recode.vcf -R -f -k -l -O -t -f
./RAiSD -n toro.outliers -I toro.outliers.recode.vcf -R -f -k -l -O -t -f
./RAiSD -n mega.normal -I mega.normal.recode.vcf -R -f -k -l -O -t -f
./RAiSD -n toro.normal -I toro.normal.recode.vcf -R -f -k -l -O -t -f

```
To establish a mu cutoff accounting for false positives due to background seleciton, we simulate a 1e6bp region

```
./sfs_code 1 10 -A -n 10 -N 2000 -L 1 1000000 -I -r 0.001 -W 1 0 0 0.20 -Td 0 2 -Td 0.75 0.25 -Tg 0.80 2.5 -s 182736 -TE 1 --VCF -o sim.mu
```

```
./RAiSD -n sim -I /media/burke/bigMac/ethan/sfscode_20150910/bin/sim.mu.vcf -k 0.05
```
The output in the report gives us a threshhold of 0.004681838:
```
 SORTED DATA (FPR 0.050000)
Size                    10
Highest Score           0.004681838
Lowest Score            0.002595457
FPR Threshold           0.004681838
Threshold Location      0
```

Finally, let's use vcftools to calculate per-site nucleotide diversity:
```
vcftools --vcf sim.mu.vcf --window-pi 1000 --window-pi-step 1000 --out sim.window.pi
vcftools --vcf syma.scaf13.mega.recode.vcf --window-pi 1000 --window-pi-step 1000 --out meg.window.pi
```

## Oneliners  

Various bash oneliners for coverage calculations, etc., are found in `oneliners.sh`.

## Plotting  

All plots can be reproduced from processed data files in `data/` using `plotting.R`. 

## Authorship  

Pipeline by me ((Ethan Linck)[https://elinck.org/]), with many contributions by [C.J. Battey](http://cjbattey.com/).

