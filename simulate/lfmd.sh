#!/bin/bash -eux
version="dev"
#export LD_LIBRARY_PATH=/hwfssz1/ST_MCHRI/CLINIC/SOFTWARES/lib:/ldfssz1/ST_HEALTH/Aging/qiyanwei/backup/sys/lib/gcc-5.2.0/lib:/ldfssz1/ST_HEALTH/Aging/qiyanwei/backup/sys/lib/gcc-5.2.0/lib64:/ldfssz1/ST_HEALTH/Aging/qiyanwei/backup/sys/lib:/opt/gridengine/lib/linux-x64::/opt/python/lib

### input output
#fq1=`readlink -f $1`; fq2=`readlink -f $2`; pre=$3
inbam=`readlink -f $1`; pre=$2

outd=`dirname $pre`
if [ ! -d $outd ]; then mkdir -p $outd; fi

### option
bamRdf_opt=""
bamDCS_opt=" -c -q 0 -Q 0 -t 3 -s 3 -S 3000 -f 0.001 -e 0.0001 "
RAM=7g
mut_opt=" -d 1 -c 0 -C 1 "
spacerLength=5
tagLength=12

### reference
#ref=/home/yerui/ref/GRCm38.75/Mus_musculus.GRCm38.75.dna.primary_assembly.fa
ref=/home/yerui/ref/mitochondria/rCRS.fa

### software
bin=/home/yerui/project/LFMD/$version
samtools=samtools
bwa=bwa
gatk=/psychipc01/disk2/software/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar #hk

# start
echo `date` start low-frequency mutation detection

:<<'_NOTE_'

_NOTE_

# filt mapQ
$samtools view -b -q 30 $inbam \
    | $samtools sort -n - -o $pre.sort.realign.mapQ30.sortByName.bam ; echo `date` sort by name done

### cluster family
$bin/bamRdf $bamRdf_opt -i $pre.sort.realign.mapQ30.sortByName.bam -o $pre.family.bam ; echo `date` cluster family done

rm -f $pre.sort.realign.mapQ30.sortByName.bam

### generate Double-strand Consensus Sequence (DCS)
$bin/bamDCS $bamDCS_opt $pre.family.bam $pre.dcs -o $pre.dcs.bam ; echo `date` dcs done

#rm -f $pre.family.bam* $pre.dcs.*.fq.gz

$samtools sort $pre.dcs.bam -o $pre.dcs.sort.bam
$samtools index $pre.dcs.sort.bam

rm -f $pre.dcs.bam

# Make a pileup file from the final DCS reads using
$samtools mpileup -B -A -d 5000000 \
    -f $ref $pre.dcs.sort.bam > $pre.dcs.pileup ; echo `date` pileup done

#rm -f $pre.dcs.sort.ba*

# Count the number of unique mutations present in the final DCS sequences and calculate their frequencies
python $bin/CountMuts.py $mut_opt -i $pre.dcs.pileup -o $pre.dcs.countmuts ; echo `date` count mut done

# Locate the genomic position of each mutation
python $bin/mut-position.py $mut_opt -i $pre.dcs.pileup -o $pre.dcs.mutpos ; echo `date` locate mut done

### Call muations
$bin/lhmut -s 1 -f 0.00001 -e 1e-7 -i $pre.dcs.pileup -o $pre.lh.mut ; echo `date` lhmut done
$bin/adjust_p2.pl $pre.lh.mut -f 1e-5 | grep -v F > $pre.lh.mut.adj

rm -f $pre.dcs.pileup

echo `date` job-done
