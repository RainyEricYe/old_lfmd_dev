source ~/anaconda3/etc/profile.d/conda.sh
conda activate lfmd

module purge
module load biopython/1.71
module load java/8.0_161

./lfmd.sh raw_1.fq.gz raw_2.fq.gz out
