#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $nFrag, $skip)
 = ( 0 ,    100,    0,);
GetOptions(
    'help|?' => \$Help,
    'n=i' => \$nFrag,
    'skip=i' => \$skip,
);
die `pod2text $0` if $Help or @ARGV < 3;

my ($inbam, $ref, $pre) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my $chr = ` head -1 $ref | sed 's/^>//' | awk '{print \$1}' `;
chomp $chr;

my $cnt = 0;
open my $IN, "samtools view $inbam |" or die $!;
open my $REG, ">$pre.tmp.reg" or die $!;

while ( <$IN> ) {
    next if $skip && $. <= $skip;

    my ( $c, $p, $size ) = (split)[2,3,8];
    next if $c ne $chr && $c ne "MT2";

    if ( $size > 100 && $size < 1000 ) {
        my $reg = "$chr:$p-".($p+$size-1);
        print $REG "$reg\n";
        $cnt++;
        last if $cnt == $nFrag;
    }
}
close $IN;
close $REG;

$cnt = 1;
my $id = "";
open my $NE, "samtools faidx $ref -r $pre.tmp.reg -n 10000 |" or die $!;
open my $OUT, ">$pre" or die $!;

while ( <$NE> ) {
    chomp;
    if ( $. % 2 == 1 ) {
        s/^>//;
        my @t = split /[:-]/, $_;
        $id = join "_", @t;
    }
    else {
        my $seq = $_;
        my $qua = "I" x length($seq);

        print $OUT "\@$id\_0:0:0_0:0:0_$cnt/1\n";
        print $OUT "$seq\n+\n$qua\n";
        $cnt++;
    }
}
close $NE;
close $OUT;

`rm -f $pre.tmp.reg`;

&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900, $t[4] + 1, @t[3,2,1,0], $_[0];
}

__END__

=head1 Function

=head1 Usage
    perl $0

=head1 Options
    -h      help

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2012-

=cut
