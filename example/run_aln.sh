source ~/anaconda3/etc/profile.d/conda.sh
conda activate lfmd

module purge
module load biopython/1.71
module load java/8.0_161

fq1=raw_1.fq.gz
fq2=raw_2.fq.gz
pre=tmp

python ../tag_to_header.py --infile1 $fq1 --infile2 $fq2 \
    --outprefix $pre \
    --tagstats \
    --spacerlen 5 \
    --taglen 12 ; echo `date` tag to head done

../align.sh $pre.seq1.smi.fq.gz $pre.seq2.smi.fq.gz tmp/t
