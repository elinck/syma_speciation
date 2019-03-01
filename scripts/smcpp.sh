#!/bin/bash
cd ~/Dropbox/selasphorus/smcpp/
/media/burke/bigMac/ethan/

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
/media/burke/bigMac/ethan/alignment/recalibrated/syma.smcpp.revised.recode.vcf.gz \
/media/burke/bigMac/ethan/smcpp\_{1}\_{2}.smc.gz \
{1} \
mega:EL18_mega,EL19_mega,EL1_mega,EL20_mega,EL23_mega,EL24_mega,EL27_mega,EL40_toro,EL4_mega,EL6_mega \
::: `cat /media/burke/bigMac/ethan/keep_contigs.txt` ::: EL24_mega EL40_toro EL1_mega

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
/media/burke/bigMac/ethan/alignment/recalibrated/syma.smcpp.revised.recode.vcf.gz \
/media/burke/bigMac/ethan/smcpp\_{1}\_{2}.smc.gz \
{1} \
toro:EL10_toro,EL11_toro,EL13_toro,EL21_toro,EL32_toro,EL39_toro,EL8_toro,EL9_toro \
::: `cat /media/burke/bigMac/ethan/keep_contigs.txt` ::: EL39_toro EL8_toro EL9_toro
