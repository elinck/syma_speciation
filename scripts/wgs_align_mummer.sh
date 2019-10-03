#!/bin/bash
cd /media/burke/bigMac/ethan/

refSeqs=tgut2_split/* # path to zebra finch reference genome folder

for chr in $refSeqs # loop over fastas in reference genome folder

    do
    echo "aligning to $chr"

    nucmer $chr /media/burke/bigMac/ethan/masked.ref.fa --maxgap 1000 # boost from default of 90 to account for structural rearrangements in divergent ref

    show-coords ./out.delta > /media/burke/bigMac/ethan/mummer_out/coords/${chr##*/}.coords #Check regex if changing paths

    mv out.delta /media/burke/bigMac/ethan/mummer_out/alignments/${chr##*/}.delta

done
