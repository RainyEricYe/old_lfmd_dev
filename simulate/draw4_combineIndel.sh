perl ./bin/blast2mutList.pl ../run/rCRS.fa ../run/fMT.fa > mut.list

ls ../*/*lh.mut | while read a;do (perl ./bin/adjust_p4.pl $a -f 1e-1 > $a.adj4.e-1);done

cat ./bin/fragN | while read a;do ( perl ./bin/tpfp_combineIndel.pl mut.list ../2.umi_tools/$a.out.dcs.mutpos);done > umi.tpfp2
cat ./bin/fragN | while read a;do ( perl ./bin/tpfp_combineIndel.pl mut.list ../3.duplex/$a.dcs.mutpos);done > ds.tpfp2
cat ./bin/fragN | while read a;do ( perl ./bin/tpfp_combineIndel.pl mut.list ../6.uniConsensus/$a.dcs.mutpos);done > uni.tpfp2
cat ./bin/fragN | while read a;do ( perl ./bin/tpfp_combineIndel.pl mut.list ../4.lfmd/$a.lh.mut.adj4.e-1);done > lfmd.adj4.e-1.tpfp2
paste ./bin/fragN ds.tpfp2 umi.tpfp2 lfmd.adj4.e-1.tpfp2 uni.tpfp2 > t2.txt

#cat ./bin/fragN | while read a;do ( perl ./bin/tpfp_combineIndel.pl mut.list2 ../4.lfmd/$a.lh.mut.adj);done > lfmd.adj.tpfp
#paste ./bin/fragN ds.tpfp umi.tpfp lfmd.adj.tpfp uni.tpfp > t3.txt
