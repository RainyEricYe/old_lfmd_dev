# LFMD
Low Frequency Mutation Detector

# Preparation

* htslib(https://github.com/samtools/htslib)
- SeqLib(https://github.com/walaj/SeqLib)
* alglib(http://www.alglib.net/translator/re/alglib-3.15.0.cpp.gpl.tgz)
- bamRdf(https://github.com/RainyEricYe/bamRdf)
* bamDCS(https://github.com/RainyEricYe/bamDCS)
- lhmut(https://github.com/RainyEricYe/lhmut)

# Install
    wget https://github.com/RainyEricYe/LFMD.git

# Contact
  yerui@connect.hku.hk
  
# Citation
Rui Ye et al. LFMD: detecting low-frequency mutations in genome sequencing data without molecular tags.
  https://www.biorxiv.org/content/10.1101/617381v9
  
# Usage
    vi lfmd2.sh && ajust parameters and file paths
    sh lfmd2.sh in.fq1 in.fq2 outprefix
