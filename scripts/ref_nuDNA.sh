#!/bin/bash

# extract reference sequence
sed '/^#/d' syma.svd.recode.vcf > nohead.syma.svd.recode.vcf
awk '{print $4}' nohead.syma.svd.recode.vcf > ref.allele
tr -d '\n' < ref.allele > h_senegalensis.fasta

# run vcf2phylip on vcf; concatenate with reference sequence. 
python vcf2phylip/vcf2phylip.py --input syma.svd.recode.vcf --fasta
cat h_senegalensis.fasta syma.svd.recode.min4.fasta > syma_nuc.fasta
