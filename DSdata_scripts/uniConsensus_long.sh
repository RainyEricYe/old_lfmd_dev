#!/bin/bash
set -e
set -u
set -x

echo `date` start unified consensus maker

fq1=$1; fq2=$2; pre=$3
#inbam=$1; pre=$2

sample=`basename $pre`
outd=`dirname $pre`
if [ ! -d $outd ]; then mkdir -p $outd; fi

picard=/home/yerui/anaconda3/share/picard-2.18.29-0/picard.jar
gatk=/psychipc01/disk2/software/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar

#ref=/home/yerui/ref/mitochondria/rCRS.fa
#ref=/home/yerui/ref/GRCm38.75/Mus_musculus.GRCm38.75.dna.primary_assembly.fa
ref=/home2/groups/pcsham/users/yerui/ref/hg19ByChr/hg19.rCRS.fa
bin=/home/yerui/project/LFMD/bin

RAM=9G
mut_opt=" -d 1 -c 0 -C 1 "
#cyclesToTrim="1-5,80-84"
cyclesToTrim="1-5,111-115"

:<<'_NOTE_'
_NOTE_

# fastq to unaligned bam
java  -Xmx7g -jar $picard FastqToSam \
    F1=$fq1 \
    F2=$fq2 \
    V=Standard \
    O=$pre.unaligned.bam \
    SM=$sample ;  echo `date` convert to bam done && \

# make dcs files
python $bin/UnifiedConsensusMaker.py \
    --input $pre.unaligned.bam \
    --taglen 12 \
    --spacerlen 5 \
    --tagstats \
    --minmem 3 \
    --maxmem 30000 \
    --cutoff 0.7 \
    --Ncutoff 0.1 \
    --rep_filt 9 \
    --prefix $pre ;  echo `date` make dcs files done && \

    rm -f $pre.unaligned.bam $pre.temp.sort.bam

# align sort realign
$bin/align.sh ${pre}_read1_dcs.fq.gz ${pre}_read2_dcs.fq.gz $pre.dcs $RAM $ref

    rm -f ${pre}_read*_dcs.fq.gz

# Perform end-trimming of DCS reads, head 5bp & tail 5bp
java -Xmx3g -jar $gatk \
    -T ClipReads \
    -R $ref \
    -I $pre.dcs.sort.realign.bam \
    -o $pre.dcs.sort.realign.trim.bam \
    --cyclesToTrim $cyclesToTrim \
    --clipRepresentation SOFTCLIP_BASES ; echo `date` dcs trim done

    rm -f $pre.dcs.sort.realign.ba*

# filt repeat of DCS
samtools view -h -q 30 $pre.dcs.sort.realign.trim.bam \
    | samtools view -bS - > $pre.dcs.final.bam ; echo `date` dcs uniq done

    rm -f $pre.dcs.sort.realign.trim.ba*

# Make a pileup file from the final DCS reads using
samtools mpileup -B -A -d 5000000 \
    -f $ref $pre.dcs.final.bam > $pre.dcs.final.pileup ; echo `date` pileup done

    rm -f $pre.dcs.final.bam

# Count the number of unique mutations present in the final DCS sequences and calculate their frequencies
python $bin/CountMuts.py $mut_opt -i $pre.dcs.final.pileup -o $pre.dcs.countmuts ; echo `date` count mut done

# Locate the genomic position of each mutation
python $bin/mut-position.py $mut_opt -i $pre.dcs.final.pileup -o $pre.dcs.mutpos ; echo `date` locate mut done

    rm -f $pre.dcs.final.pileup
# end
echo `date` job-done
