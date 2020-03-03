#perl ./bin/blast2mutList.pl ../run/rCRS.fa ../run/fMT.fa > mut.list

ls ../*/*lh.mut | while read a;do (perl ./bin/adjust_p4.pl $a -f 1e-1 > $a.adj4.e-1);done

cat ./bin/fragN | while read a;do ( perl ./bin/tpfp.pl mut.list ../2.umi_tools/$a.out.dcs.mutpos);done > umi.tpfp
cat ./bin/fragN | while read a;do ( perl ./bin/tpfp.pl mut.list ../3.duplex/$a.dcs.mutpos);done > ds.tpfp
cat ./bin/fragN | while read a;do ( perl ./bin/tpfp.pl mut.list ../6.uniConsensus/$a.dcs.mutpos);done > uni.tpfp
cat ./bin/fragN | while read a;do ( perl ./bin/tpfp.pl mut.list ../4.lfmd/$a.lh.mut.adj4.e-1);done > lfmd.adj4.e-1.tpfp
paste ./bin/fragN ds.tpfp umi.tpfp lfmd.adj4.e-1.tpfp uni.tpfp > t3.txt

#cat ./bin/fragN | while read a;do ( perl ./bin/tpfp.pl mut.list2 ../4.lfmd/$a.lh.mut.adj);done > lfmd.adj.tpfp
#paste ./bin/fragN ds.tpfp umi.tpfp lfmd.adj.tpfp uni.tpfp > t3.txt
