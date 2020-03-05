#!/bin/bash -eux

read inbam pre hash <<< "$@"

outd=`dirname $pre`
if [ ! -d $outd ]; then mkdir -p $outd; fi

RAM=9G
mut_opt=" -d 1 -c 0 -C 1 "
cyclesToTrim="1-5,135-139"

spacerLength=5
tagLength=12
iSize=-1
minMem=3
maxMem=200
cutOff=0.7
nCutOff=0.1
readLength=139
filtersSet=sn
readTypes=d

bin=/home/yerui/project/LFMD/bin
#ref=/home/yerui/ref/mitochondria/rCRS.fa
ref=/home2/groups/pcsham/users/yerui/ref/hg19ByChr/hg19.rCRS.fa
#ref=/home/yerui/ref/GRCm38.75/Mus_musculus.GRCm38.75.dna.primary_assembly.fa
gatk=/psychipc01/disk2/software/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar

# start
echo `date` start low-frequency mutation detection

:<<'_NOTE_'

_NOTE_
samtools view -b -q 30 $inbam > $pre.sort.realign.mapQ30.bam ; echo `date` sort done

# Run Consensus Maker
python $bin/ConsensusMaker.py \
    --infile $pre.sort.realign.mapQ30.bam \
    --tag_file $pre.pe.tagcounts \
    --outfile $pre.sscs.bam \
    --minmem $minMem \
    --maxmem $maxMem \
    --read_length $readLength \
    --cut_off $cutOff \
    --Ncut_off $nCutOff \
    --read_type $readTypes \
    --filt $filtersSet \
    --isize $iSize \
    --tag_stat $pre.pe.tagstat ; echo `date` sscs done

rm -f $pre.sort.realign.mapQ30.bam $pre.pe.tagcounts

# Sort SSCSs
samtools view -bu $pre.sscs.bam | samtools sort - -o $pre.sscs.sort.bam ; echo `date` sscs sort done

# Run Duplex Maker
python $bin/DuplexMaker.py \
    --infile $pre.sscs.sort.bam \
    --outfile $pre.dcs.bam \
    --Ncutoff $nCutOff \
    --readlength $readLength ; echo `date` dcs done

rm -f $pre.sscs*bam

# align sort realign
$bin/align.sh $pre.dcs.r1.fq $pre.dcs.r2.fq $pre.dcs $RAM $ref

rm -f $pre.dcs.r*.fq

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
