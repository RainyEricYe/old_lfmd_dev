#!/bin/bash -eux

version="v0.4"

#export LD_LIBRARY_PATH=/hwfssz1/ST_MCHRI/CLINIC/SOFTWARES/lib:/ldfssz1/ST_HEALTH/Aging/qiyanwei/backup/sys/lib/gcc-5.2.0/lib:/ldfssz1/ST_HEALTH/Aging/qiyanwei/backup/sys/lib/gcc-5.2.0/lib64:/ldfssz1/ST_HEALTH/Aging/qiyanwei/backup/sys/lib:/opt/gridengine/lib/linux-x64::/opt/python/lib

### input output
fq1=`readlink -f $1`; fq2=`readlink -f $2`; pre=$3
#inbam=`readlink -f $1`; pre=$2

outd=`dirname $pre`
if [ ! -d $outd ]; then mkdir -p $outd; fi
pre=`readlink -f $pre`

### option
bamRdf_opt=""
bamDCS_opt=" -c -q 0 -Q 0 -t 3 -s 3 -S 3000 -f 0.001 -e 0.0001 "
trim_opt="1-5,80-84"
mut_opt=" -d 1 -c 0 -C 1 "
RAM=7g

### reference
#ref=/hwfssz1/ST_HEALTH/Aging/F16ZQSB1SY2761/yerui/ref/rCRS.fa
#ref=/home2/groups/pcsham/users/yerui/ref/hg19ByChr/hg19.rCRS.fa
#ref=/home/yerui/ref/GRCm38.75/Mus_musculus.GRCm38.75.dna.primary_assembly.fa
#ref=/zfssz2/ST_HEALTH/HEALTH/yerui_5T/MacacaMtDNA/pipeline/LFMD/ref/macaca_fascicularis_mtdna.fa
ref=/home/yerui/ref/mitochondria/rCRS.fa

### software
#bin=/hwfssz1/ST_HEALTH/Aging/F16ZQSB1SY2761/yerui/project/LFMD/$version
bin=/home/yerui/project/LFMD/$version #hk

samtools=$bin/samtools
bwa=$bin/bwa
#gatk=/home/yerui/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
gatk=/psychipc01/disk2/software/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar #hk


# start
echo `date` start low-frequency mutation detection

:<<'_NOTE_'

_NOTE_

spacerLength=5
tagLength=12

# move tags to read ID
python $bin/tag_to_header.py \
    --infile1 $fq1 \
    --infile2 $fq2 \
    --outprefix $pre \
    --tagstats \
    --spacerlen $spacerLength \
    --taglen $tagLength ; echo `date` tag to head done


#ln -sf $fq1 $pre.seq1.smi.fq.gz
#ln -sf $fq2 $pre.seq2.smi.fq.gz

# alignment
$bwa aln $ref $pre.seq1.smi.fq.gz > $pre.seq1.aln
$bwa aln $ref $pre.seq2.smi.fq.gz > $pre.seq2.aln

$bwa sampe -s $ref -r '@RG\tID:foo\tSM:bar' \
    $pre.seq1.aln $pre.seq2.aln $pre.seq1.smi.fq.gz $pre.seq2.smi.fq.gz \
    | $samtools view -bS - > $pre.bam ; echo `date` bwa done

# filt unmap & sort by position
$samtools view -b -F 4 $pre.bam | $samtools sort - -o $pre.sort.bam ; echo `date` filt unmap and sort done
$samtools index $pre.sort.bam ; echo `date` index done

# identify the genome targets for local re-alignment
java -Xmx$RAM -jar $gatk \
    -T RealignerTargetCreator \
    -R $ref \
    -I $pre.sort.bam \
    -o $pre.sort.intervals ; echo `date` creat target done

# local re-alignment
java -Xmx$RAM -jar $gatk \
    -T IndelRealigner \
    -R $ref \
    --maxReadsForConsensuses 30000 \
    --maxReadsForRealignment 300000 \
    --maxReadsInMemory 300000 \
    -I $pre.sort.bam \
    -targetIntervals $pre.sort.intervals \
    -o $pre.sort.realign.bam ; echo `date` re-align done

# filt for duplex
$samtools view -b -q 30 $pre.sort.realign.bam > $pre.sort.realign.mapQ30.bam ; echo `date` sort done

# filt mapQ
$samtools view -b -q 30 $pre.sort.realign.bam \
    | $samtools sort -n - -o $pre.sort.realign.mapQ30.sortByName.bam ; echo `date` sort by name done

# cluster family
$bin/bamRdf $bamRdf_opt -i $pre.sort.realign.mapQ30.sortByName.bam -o $pre.family.bam ; echo `date` cluster family done

# generate Doubel-strand Consensus Sequence (DCS)
$bin/bamDCS $bamDCS_opt $pre.family.bam $pre.dcs -o $pre.dcs.bam ; echo `date` dcs done

$samtools sort $pre.dcs.bam -o $pre.dcs.sort.bam
$samtools index $pre.dcs.sort.bam

:<<'TMP'

TMP
# identify the genome targets for local re-alignment
java -Xmx$RAM -jar $gatk \
    -T RealignerTargetCreator \
    -R $ref \
    -I $pre.dcs.sort.bam \
    -o $pre.dcs.sort.intervals ; echo `date` dcs creat target done

# local re-alignment
java -Xmx$RAM -jar $gatk \
    -T IndelRealigner \
    -R $ref \
    --maxReadsForConsensuses 30000 \
    --maxReadsForRealignment 300000 \
    --maxReadsInMemory 300000 \
    -I $pre.dcs.sort.bam \
    -targetIntervals $pre.dcs.sort.intervals \
    -o $pre.dcs.sort.realign.bam ; echo `date` dcs re-align done

# Perform end-trimming of DCS reads
java -Xmx3g -jar $gatk \
    -T ClipReads \
    -R $ref \
    -I $pre.dcs.sort.realign.bam \
    -o $pre.dcs.final.bam \
    --cyclesToTrim $trim_opt \
    --clipRepresentation SOFTCLIP_BASES ; echo `date` dcs trim done

# Make a pileup file from the final DCS reads using
$samtools mpileup -B -A -d 5000000 -f $ref $pre.dcs.final.bam > $pre.dcs.final.pileup ; echo `date` pileup done

# Count the number of unique mutations present in the final DCS sequences and calculate their frequencies
python $bin/CountMuts.py $mut_opt -i $pre.dcs.final.pileup -o $pre.dcs.countmuts ; echo `date` count mut done

# Locate the genomic position of each mutation
python $bin/mut-position.py $mut_opt -i $pre.dcs.final.pileup -o $pre.dcs.mutpos ; echo `date` locate mut done

echo `date` job-done
