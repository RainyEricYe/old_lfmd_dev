#!/bin/bash -eux

inRef=$1
outpre=$2
copy=$3

sim=/home/yerui/src/dunovo/utils/sim2_endHighErr.py
getFrag=/home/yerui/project/LFMD/dev/getFrag_randShear.pl

perl $getFrag -n $copy $inRef $outpre.fragFile.fq

module load biopython/1.71

runit python $sim --frag-file $outpre.fragFile.fq \
    -1 $outpre\_1.fq \
    -2 $outpre\_2.fq \
    -o fastq \
    -Q 5 \
    -n $copy \
    -r 100 \
    -s 0.01 \
    -p 0.00001 \
    -c 20 \
    -B 12 \
    -I TGACT \
    -l $outpre.log \
    -q

rm -rf $outpre.fragFile.fq

#    -m $outpre.mutations \
#    -b $outpre.barcodes \

echo job-done
