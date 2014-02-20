#!/bin/bash

# trainAllModels.sh

for SP in ailMel1 anoCar2 anoGam1 apiMel2 bosTau4 calJac3 canFam2 cavPor3 ce3 ci2 danRer7 dm3 equCab2 felCat4 fr2 galGal3 gasAcu1 loxAfr3 mm10 mm9 monDom5 ornAna1 oryCun2 oryLat2 oviAri1 panTro3 ponAbe2 rheMac2 rn4 strPur2 susScr2 taeGut1 tetNig2 xenTro2
do
    echo "doing $SP"
    for WIN in 1000
    do
        echo "  using window $win"
        bin/createModelBin.pl -m $SP -w $WIN -d data2 -e chrM --keep_dw_files --no_kmer_table -v
        for KMER in 4 6 8
        do 
            echo "    using kmer $KMER"
            bin/createModelBin.pl -m $SP -w $WIN -k $KMER -d data2 -e chrM --keep_dw_files --no_repeat_table -v
        done
    done
done
