#!/bin/bash/
#see select_contigs.R for scripts to write bootstrapped lists of *.smc.gz files to ./boot/

cd /media/burke/bigMac/ethan/smcpp/
mkdir tmp

#megarhyncha
bootstraps=./mega_boot/*
for boot in $bootstraps
do
    smc++ estimate -o tmp/ -c 30000 --cores 30 --unfold \
    --thinning 400 2.3e-9 `cat $boot`

    mv tmp/model.final.json models_mega/`basename $boot .txt`.model.final.json
    rm tmp/*
done

#torotoro
bootstraps=./toro_boot/*
for boot in $bootstraps
do
    smc++ estimate -o tmp/ -c 30000 --cores 30 --unfold \
    --thinning 400 2.3e-9 `cat $boot`

    mv tmp/model.final.json models/sed/`basename $boot .txt`.model.final.json
    rm tmp/*
done

rm -r tmp
