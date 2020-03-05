 perl bin/form.pl bin/smp.info > Sample.Class.Label.Dementia.Stage.Tissue.txt

 perl -e 'print "Sample\tClass\tLabel\tDementia\tStage\tTissue\tTools\tMutation\tCount\n"' > mut_info.txt

 less Sample.Class.Label.Dementia.Stage.Tissue.txt | perl -lane '@a=`perl bin/mut_spectrum.pl ../4.lfmd/$F[0]*mut.adj4.e-1`; chomp @a; for $i (@a) {print "$_\tLFMD\t$i"}' >> mut_info.txt
 less Sample.Class.Label.Dementia.Stage.Tissue.txt | perl -lane '@a=`perl bin/mut_spectrum.pl ../3.duplex/$F[0]*pos`; chomp @a; for $i (@a) {print "$_\tDS\t$i"}' >>mut_info.txt
 less Sample.Class.Label.Dementia.Stage.Tissue.txt | perl -lane '@a=`perl bin/mut_spectrum.pl ../6.uniConsensus/$F[0]*pos`; chomp @a; for $i (@a) {print "$_\tUniC\t$i"}' >> mut_info.txt
 less Sample.Class.Label.Dementia.Stage.Tissue.txt | perl -lane '@a=`perl bin/mut_spectrum.pl ../2.umi_tools/$F[0]*pos`; chomp @a; for $i (@a) {print "$_\tUMI\t$i"}' >> mut_info.txt
