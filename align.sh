#!/bin/bash -eu
version="dev"
gatk=/psychipc01/disk2/software/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar # hk

read fq1 fq2 pre RAM ref hash <<< "$@"

if [ $# -lt 3 ]; then
    echo "Usage: $0 fq1 fq2 out_prefix [RAM] [ref]"
    exit 1
fi

if [ -z $RAM ]; then RAM=7G; fi
if [ -z $ref ]; then
    ref=/home/yerui/ref/GRCm38.75/Mus_musculus.GRCm38.75.dna.primary_assembly.fa
fi

outd=`dirname $pre`
if [ ! -d $outd ]; then mkdir -p $outd; fi

bwa aln $ref $fq1 > $pre.r1.aln
bwa aln $ref $fq2 > $pre.r2.aln

bwa sampe -s $ref -r '@RG\tID:foo\tSM:bar' $pre.r1.aln $pre.r2.aln $fq1 $fq2 \
    | samtools view -bS - > $pre.bam ; echo `date` bwa done

samtools view -b -F 4 $pre.bam | samtools sort - -o $pre.sort.bam ; echo `date` filt unmap and sort done
samtools index $pre.sort.bam ; echo `date` index done

java -Xmx$RAM -jar $gatk \
    -T RealignerTargetCreator \
    -R $ref \
    -I $pre.sort.bam \
    -o $pre.sort.intervals ; echo `date` creat target done

java -Xmx$RAM -jar $gatk \
    -T IndelRealigner \
    -R $ref \
    --maxReadsForConsensuses 30000 \
    --maxReadsForRealignment 300000 \
    --maxReadsInMemory 300000 \
    -I $pre.sort.bam \
    -targetIntervals $pre.sort.intervals \
    -o $pre.sort.realign.bam ; echo `date` re-align done

rm -f $pre.r1.aln $pre.r2.aln $pre.bam $pre.sort.bam $pre.sort.bam.bai $pre.sort.intervals
echo `date` clean tmp files done
