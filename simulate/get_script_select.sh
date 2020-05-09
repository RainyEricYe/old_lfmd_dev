#!/bin/bash -eux

cd ..
mkdir 0.data  1.dunovo  2.umi_tools  3.duplex  4.lfmd  5.lfmd.s2  6.uniConsensus  compare  stat
cd run
ln -sf /home/yerui/ref/mitochondria/rCRS.fa

# fake mutation
#perl ./bin/fakeMut_rand.pl rCRS.fa fMT.fa fMT.fa.mut_list -snvNum 100 -insNum 10 -delNum 10 -insMaxLen 3 -delMaxLen 3

# simulation
# cat ./bin/fragN | perl -lne '$b=1e6-$_; print "./bin/simu.randomShear_frag.sh fMT.fa ../0.data/f.$_ $_ ; ./bin/simu.randomShear_frag.sh rCRS.fa ../0.data/t.$_ $b ; cat ../0.data/f.$_\_1.fq >> ../0.data/t.$_\_1.fq ; gzip ../0.data/t.$_\_1.fq ; cat ../0.data/f.$_\_2.fq >> ../0.data/t.$_\_2.fq ; gzip ../0.data/t.$_\_2.fq ; rm -f ../0.data/f.$_\_*.fq ; echo job-done"' > run_simu.sh
# split.pl run_simu.sh -p simu -q small -m 8G -W 6:00:00

batch=$1
frac=$2

# select 0.5
  cat ./bin/fragN | while read a;do echo perl ./bin/randomSelect.pl ../../../2.simulate.$batch/0.data/t.$a\_1.fq.gz ../../../2.simulate.$batch/0.data/t.$a\_2.fq.gz ../0.data/t.$a\_1.fq.gz ../0.data/t.$a\_2.fq.gz -f $frac;done > run_select$frac.sh
 split.pl run_select$frac.sh -p se -q small -m 4G -W 2:00:00

#  aln the raw fq
 cat ./bin/fragN | perl -lne '$fq1="../0.data/t.$_\_1.fq.gz"; $fq2="../0.data/t.$_\_2.fq.gz"; $pre="../0.data/t.$_"; print "source ~/anaconda3/etc/profile.d/conda.sh; conda activate lfmd; module purge; module load biopython/1.71; module load java/8.0_161; python /home/yerui/project/LFMD/dev/tag_to_header.py --infile1 $fq1 --infile2 $fq2 --outprefix $pre --tagstats --spacerlen 5 --taglen 12 && /home/yerui/project/LFMD/dev/align.sh $pre.seq1.smi.fq.gz $pre.seq2.smi.fq.gz $pre 7G /home/yerui/ref/mitochondria/rCRS.fa && echo job-done"' > run_aln.sh
 split.pl run_aln.sh -p aln -q small -m 10G -W 2:00:00

# run dunovo from fq
 cat ./bin/fragN | perl -lne ' print "source ~/anaconda3/etc/profile.d/conda.sh; conda activate dunovo; ./bin/dunovo.sh ../0.data/t.$_\_1.fq.gz ../0.data/t.$_\_2.fq.gz ../1.dunovo/$_ && echo job-done" ' > run_dunovo.sh
 split.pl run_dunovo.sh -p dn -q small -m 10G -W 6:00:00

# run umi_tools from bam
 cat ./bin/fragN | perl -lne ' print "source ~/anaconda3/etc/profile.d/conda.sh; conda activate base; ./bin/umi_tools.sh ../0.data/t.$_.sort.realign.bam ../2.umi_tools/$_ && source ~/anaconda3/etc/profile.d/conda.sh; conda activate lfmd; module purge; module load biopython/1.71; module load java/8.0_161; ./bin/duplex.sh ../2.umi_tools/$_.umi.sort.bam ../2.umi_tools/$_.out && rm -f ../2.umi_tools/$_.umi.sort.bam* ; echo job-done " ' > run_umi_tools.sh
 split.pl run_umi_tools.sh -p umi -q medium -m 20G -W 2:00:00

# run duplex from bam
 cat ./bin/fragN | perl -lne ' print "source ~/anaconda3/etc/profile.d/conda.sh; conda activate lfmd; module purge; module load biopython/1.71; module load java/8.0_161; ./bin/duplex.sh ../0.data/t.$_.sort.realign.bam ../3.duplex/$_ && echo job-done" ' > run_duplex.sh
 split.pl run_duplex.sh -p ds -q small -m 10G -W 2:00:00

# run lfmd from bam
 cat ./bin/fragN | perl -lne ' print "source ~/anaconda3/etc/profile.d/conda.sh; conda activate lfmd; module purge; module load biopython/1.71; module load java/8.0_161; ./bin/lfmd.sh ../0.data/t.$_.sort.realign.bam ../4.lfmd/$_ && echo job-done" ' > run_lfmd.sh
 split.pl run_lfmd.sh -p lfmd -q small -m 10G -W 2:00:00

# run uniConsensus from fq
 cat ./bin/fragN | perl -lne 'print " source ~/anaconda3/etc/profile.d/conda.sh; conda activate lfmd; module purge; module load biopython/1.71; module load java/8.0_161; ./bin/uniConsensus.sh ../0.data/t.$_\_1.fq.gz ../0.data/t.$_\_2.fq.gz ../6.uniConsensus/$_ " ' > run_uniConsensus.sh
 split.pl run_uniConsensus.sh -p un -q small -m 10G -W 6:00:00

# freebayes for dunovo
# split.pl run_freebayes_after_dunovo.sh -p free -q small -m 4G -W 0:30:00

# stat
  cat ./bin/fragN | while read a;do echo sh ./bin/stat.sh ../0.data/t.$a.sort.realign.bam ../stat/$a;done > run_stat.sh
  split.pl run_stat.sh -p stat -q small -m 8G -W 2:00:00

