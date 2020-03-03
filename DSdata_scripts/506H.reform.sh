samtools sort -n Human_506H.pe.sort.bam -o Human_506H.pe.sortByName.bam
samtools fastq -1 seq1.smi.fq.gz -2 seq2.smi.fq.gz  -0 /dev/null -s /dev/null -N Human_506H.pe.sortByName.bam
perl bin/tag2read.pl seq1.smi.fq.gz seq1.fq.gz && echo job-done
perl bin/tag2read.pl seq2.smi.fq.gz seq2.fq.gz && echo job-done
