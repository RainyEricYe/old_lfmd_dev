#!/bin/bash -eux
# pair end
#
#source ~/anaconda3/etc/profile.d/conda.sh
#conda activate base
#module purge

read inbam pre hash <<< "$@"

echo `date` start umi_tools

# start from sort.bam
 umi_tools group --umi-separator="|" --paired \
     --unpaired-reads=discard --mapping-quality=30 \
     --chimeric-pairs=discard --soft-clip-threshold=5 \
     --edit-distance-threshold=2 \
     -I $inbam \
     --group-out=$pre.groups.tsv \
     --output-bam -S $pre.mapped_grouped.bam \
     -L $pre.group.log -E $pre.group.err ; echo `date` group UMI done

 samtools sort -n $pre.mapped_grouped.bam -o $pre.mapped_grouped.sortByName.bam

 rm -f $pre.mapped_grouped.bam $pre.groups.tsv

 perl /home/yerui/project/LFMD/dev/updateReadID.pl $pre.mapped_grouped.sortByName.bam | samtools view -bS - > $pre.umi.bam

 rm -f $pre.mapped_grouped.sortByName.bam

 samtools sort $pre.umi.bam -o $pre.umi.sort.bam
 samtools index $pre.umi.sort.bam

 rm -f $pre.umi.bam

# unbelievable request
# module load biopython/1.71
# module load java/8.0_161

 #./duplex.sh $pre.umi.sort.bam $pre.out
