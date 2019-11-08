# *Syma* bioinformatics pipeline and analysis

This digital notebook documents commands run to analyse hyRAD and whole genome sequence data from *Syma torotoro* and
*S. megarhyncha.* Downstream genotypic data manipulation, analysis of phenotypic data, and plotting in `scripts/`.
Paths and descriptions of scripts are found in `README.md.`

## Sequence alignment and variant calling pipeline

Requires `R`, [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/), [GATK](https://software.broadinstitute.org/gatk/), [Picard tools](https://broadinstitute.github.io/picard/),
[mapDamage](https://ginolhac.github.io/mapDamage/), [MitoFinder](https://github.com/RemiAllio/MitoFinder), and [MAFFT](https://www.ebi.ac.uk/Tools/msa/mafft/).

Full commands available in `assembly.R`. Probably best to tweak and run in discrete chunks to debug rather than via `Rscript.`

Mask repeats from reference genome:

```
bbmask.sh in=B10K-DU-024-03.genomic.fa out=masked.ref.fa entropy=0.7
```

Open R via the terminal an install and load `magrittr`.

```R
setwd("/media/burke/bigMac/ethan/")
library(magrittr)
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/");
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

To generate the mitochondrial DNA alignment, run the python script `scripts/msa.py`, which will loop over individuals to assemble mitochondrial genomes, extract ND2, and run MAFFT to align sequences.

```
python msa.py
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


To create the .vcf file later formatted for moments (see below), we'll drop sites with more than 10% missing data, fewer than 6 reads per individual, and a lower quality score than 30. We'll then thin it for LD, including only sites every 50 kbp.  

```
vcftools --gzvcf syma_recal_gatk_all.vcf.gz --max-missing 0.9 --remove hyrad.bams --minQ 30 --minDP 6  --recode --recode-INFO-all --out syma.wgs.6x.moments
vcftools --vcf syma.wgs.6x.moments.recode.vcf --thin 50000 --recode --recode-INFO-all --out /media/burke/bigMac/ethan/moments_revision/syma.6x.moments.unlinked
```

(After filtering, kept 4492 out of a possible 78882912 sites.)


For sliding window summary statistics, we'll call invariant sites along with variants for our largest (huge) `.vcf`, with a minQ score of 30 and a minDP of 3:

```
vcftools --gzvcf syma_recal_gatk_all.vcf.gz --remove hyrad.bams --minQ 30 --minDP 3  --recode --recode-INFO-all --out syma.wgs.all
```

To create a SNP alignment for SVD quartets, we'll drop sites with more than 10% missing data, use a minimum minor allele count of 3, a minimum depth of 5x, and a maximum depth of 120x:

```
vcftools --vcf syma_recal_gatk.vcf --remove hyrad.bams --max-missing 0.9 --mac 3 --max-alleles 2 --minDP 5 --maxDP 120 --minQ 30 --recode --recode-INFO-all --out syma.svd
```

(After filtering, kept 283485 out of a possible 78882912 sites.)

## Morphology and bioacoustic analysis

We used the `R` scripts `morphology_analysis.R` and `vocalization_analysis.R` to analyse morphological data (`data/morphology.csv`)
and bioacoustic data (`data/params.csv`), respectively.

## Genotype analysis

We used the `R` script `genotype_analysis` to process `.vcf` data and prepare datasets the analyses below.

## BEAST

We ran BEAST via the [CIPRES Science Gateaway](https://www.phylo.org/) using the .xml file `data/syma_nd2_final.xml`, which we generated in BEAUTi.

## SVDquartets

To generate the nuclear DNA phylogeny without species assignments, we use the bash script `scripts/ref_nuDNA.sh`. This script extracts the reference sequence from a .vcf, turns the the .vcf into a fasta file using [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip) (required!), and then concatenates the reference sequence with the fasta. As they come from the same source, there is no additional alignment necessary.

```
bash ref_nuDNA.sh
```

We then run SVDquartets via the [CIPRES Science Gateaway](https://www.phylo.org/) using the nexus file `data/syma_nuc.nexus`. To infer species tree using SVDquartets, we swap this file out for the file `data/syma_sptree.nexus`.


## *D*-tests

To perform bootstrapped *D*-tests (ABBA-BABA tests), we run the python script `scripts/d_tests.py`:

```
python d_tests.py
```

## Moments

Requires [Moments](https://bitbucket.org/simongravel/moments).  

To perform demographic inference in moments, we'll run each the IM model script (`scripts/IM_moments.py)`, then use the output in a second script to generate bootstrap replicates (`scripts/IM_moments.py)`, import and process the output from this in `genotype_analysis.R`, and then plot in `plotting.R`.

```
python IM_moments.py;
python IM_moments.py;
```

## Windowed stats

Requires scripts from Simon Martin's [`genomics_general`](https://github.com/simonhmartin/genomics_general) repository.  


First, we use the script `parseVCF.py` to convert our vcf (with invariant sites) to the `.geno` custom file format.
```
python VCF_processing/parseVCF.py -i /media/burke/bigMac/ethan/alignment/recalibrated/syma.wgs.all.recode.vcf.gz --skipIndels | gzip > all.output.geno.gz
```

We then run `popgenWindows.py` with coordinate-based 50kb windows with a 10kb step size:
```
python popgenWindows.py -w 50000 -s 10000 --windType coordinate -g all.output.geno.gz -o syma_windows.csv.gz -f phased -T 24 -p meg EL19_mega,EL18_mega,EL20_mega,EL1_mega,EL23_mega,EL24_mega,EL27_mega,EL40_toro,EL4_mega,EL6_mega -p toro EL10_toro,EL11_toro,EL13_toro,EL21_toro,EL32_toro,EL39_toro,EL8_toro,EL9_toro;
```

## Oneliners  

Various bash oneliners for coverage calculations, etc., are found in `oneliners.sh`.

## Plotting  

All plots can be reproduced from processed data files in `data/` using `plotting.R`.

## Authorship  

Pipeline by me ([Ethan Linck](https://elinck.org/)), with many contributions by [C.J. Battey](http://cjbattey.com/).
