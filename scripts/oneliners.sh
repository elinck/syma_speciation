# bash one-liners used in Linck et al. in prep

# determine average depth per sample
for bam_file in *.bam
do
average_depth=$(samtools depth $bam_file |  awk '{sum+=$3} END { print "Average = ",sum/NR}')
echo "$bam_file $average_depth"
done

# check if bams are valid
for bam_file in *.bam
do
gunzip -t $bam_file && echo "$bam_file VALID"
done

# get number of reads per sample
for bam_file in *.bam
do
samtools view -c $bam_file
done

#1-index for angsd 
awk '{$2+=1;$3+=1}1' OFS='\t' shared.keep

# remove headers
for bam_file in *.bam
do
samtools view -H $bam_file | grep "^@RG" | samtools reheader - $bam_file > $bam_file
done

java -Xmx40g -jar /media/burke/bigMac/Dropbox/tools/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /media/burke/bigMac/ethan/masked.ref.fa \
-I EL10_toro_R1_sorted.bam \
-I EL11_toro_R1_sorted.bam \
-I EL13_toro_R1_sorted.bam \
-I EL18_mega_R1_sorted.bam \
-I EL19_mega_R1_sorted.bam \
-I EL1_mega_R1_sorted.bam \
-I EL20_mega_R1_sorted.bam \
-I EL21_toro_R1_sorted.bam \
-I EL23_mega_R1_sorted.bam \
-I EL24_mega_R1_sorted.bam \
-I EL27_mega_R1_sorted.bam \
-I EL29_ochr_R1_sorted.bam \
-I EL32_toro_R1_sorted.bam \
-I EL39_toro_R1_sorted.bam \
-I EL40_toro_R1_sorted.bam \
-I EL41_toro_R1_sorted.bam \
-I EL42_toro_R1_sorted.bam \
-I EL43_toro_R1_sorted.bam \
-I EL44_toro_R1_sorted.bam \
-I EL45_toro_R1_sorted.bam \
-I EL46_toro_R1_sorted.bam \
-I EL47_toro_R1_sorted.bam \
-I EL48_ochr_R1_sorted.bam \
-I EL49_toro_R1_sorted.bam \
-I EL4_mega_R1_sorted.bam \
-I EL50_toro_R1_sorted.bam \
-I EL51_toro_R1_sorted.bam \
-I EL52_toro_R1_sorted.bam \
-I EL53_toro_R1_sorted.bam \
-I EL54_toro_R1_sorted.bam \
-I EL55_toro_R1_sorted.bam \
-I EL56_toro_R1_sorted.bam \
-I EL57_toro_R1_sorted.bam \
-I EL58_toro_R1_sorted.bam \
-I EL59_toro_R1_sorted.bam \
-I EL5_ochr_R1_sorted.bam \
-I EL60_toro_R1_sorted.bam \
-I EL6_mega_R1_sorted.bam \
-I EL8_toro_R1_sorted.bam \
-I EL9_toro_R1_sorted.bam \
-o forIndelRealigner.intervals

gmap_build -D gmap_ref/ -d gmap_refgenome masked.ref.fa

gmap -D gmap_ref/ -d gmap_refgenome -f samse -t 16 uces.fa | samtools view -Shb - | samtools sort - alignment
samtools index aligment.bam

bedtools bamtobed -i alignment.bam > uces.bed

./atlas task=recal bam=../alignment/EL10_toro_realigned.bam pmdFile=../alignment/EL10_toro_realigned_PMD_input_Empiric.txt regions=../alignment/uces.bed equalBaseFreq verbose
