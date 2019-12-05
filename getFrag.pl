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
die `pod2text $0` if $Help or @ARGV < 2;

my ($inbam, $ref) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my $chr = ` head -1 $ref | sed 's/^>//' | awk '{print \$1}' `;
chomp $chr;

my $cnt = 0;
open my $IN, "samtools view $inbam |" or die $!;
while ( <$IN> ) {
    next if $skip && $. <= $skip;

    my ( $c, $p, $size ) = (split)[2,3,8];
    next if $c ne $chr && $c ne "MT2";

    if ( $size > 100 && $size < 1000 ) {
        my $reg = "$chr:$p-".($p+$size-1);
        my $str = `samtools faidx $ref $reg`;
        $str =~ s/([ACGTNacgtn])\n/$1/g;

        my $seq = (split /\n/, $str)[1];
        my $qua = "I" x length($seq);

        $cnt++;
        print "\@$chr\_$p\_".($p+$size-1)."_0:0:0_0:0:0_$cnt/1\n";
        print "$seq\n+\n$qua\n";
        last if $cnt == $nFrag;
    }
}
close $IN;




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
