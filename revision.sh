pre=$1
depth=$2
low=$3
high=$4

bin=/home/yerui/project/LFMD/v0.3

# Count the number of unique mutations present in the final DCS sequences and calculate their frequencies
cat $pre.dcs.final.pileup | python $bin/CountMuts.py -d $depth -c $low -C $high -u > $pre.countmuts.d${depth}-c${low}-${high} && echo count done && \

# Locate the genomic position of each mutation
python $bin/mut-position.py -i $pre.dcs.final.pileup -o $pre.mutpos.d${depth}-c${low}-${high} -d $depth -c $low -C $high && echo locat mutation done
