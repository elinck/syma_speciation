#!/bin/bash

### sliding window popgen from simon martin (https://github.com/simonhmartin/genomics_general)

# convert vcf to .geno
python /media/burke/bigMac/ethan/genomics_general/VCF_processing/parseVCF.py -i syma_recal_gatk.vcf -o syma.windows.geno.gz --skipIndels --minQual 30 --gtf flag=DP min=3 | gzip > syma.output.geno.gz

# run sliding window popgen stats in 50kb windows (note EL40_toro = megarhyncha)
python /media/burke/bigMac/ethan/genomics_general/popgenWindows.py -w 50000 --windType coordinate -g syma.windows.geno -o syma_windows_full.csv.gz -f phased -T 20 -p meg EL19_mega,EL18_mega,EL20_mega,EL1_mega,EL23_mega,EL24_mega,EL27_mega,EL40_toro,EL4_mega,EL6_mega -p toro EL10_toro,EL11_toro,EL13_toro,EL21_toro,EL32_toro,EL39_toro,EL8_toro,EL9_toro;

# phyml 50kb gene tree estimates
python phyml_sliding_windows.py -T 5 -g /media/burke/bigMac/ethan/alignment/sorted/syma.final.geno.gz --prefix output.phyml.win50 -w 50000 -M 200 --phyml /media/burke/bigMac/ethan/genomics_general/PhyML-3.1/PhyML-3.1_linux64 --model GTR;

# topology weighting -- work in progress
python twisst.py -t /media/burke/bigMac/ethan/genomics_general/output.phyml.w50.trees.gz -w output.topology.weights.csv.gz -g meg EL19_mega_A,EL18_mega_A,EL20_mega_A,EL1_mega_A,EL23_mega_A,EL24_mega_A,EL27_mega_A,EL40_toro_A,EL4_mega_A,EL6_mega_A,EL19_mega_B,EL18_mega_B,EL20_mega_B,EL1_mega_B,EL23_mega_B,EL24_mega_B,EL27_mega_B,EL40_toro_B,EL4_mega_B,EL6_mega_B3 -g tor EL10_toro_A,EL11_toro_A,EL13_toro_A,EL21_toro_A,EL32_toro_A,EL39_toro_A,EL8_toro_A,EL9_toro_A,EL10_toro_B,EL11_toro_B,EL13_toro_B,EL21_toro_B,EL32_toro_B,EL39_toro_B,EL8_toro_B,EL9_toro_B --method complete