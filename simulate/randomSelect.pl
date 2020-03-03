#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $Frac, )
=  ( 0,     0.001, );

GetOptions(
    'help|?' => \$Help,
    'Frac=f' => \$Frac,
);
die `pod2text $0` if $Help or @ARGV < 4;

my ($fq1, $fq2, $of1, $of2) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

open my $FQ1, "gzip -cd $fq1 | " or die $!;
open my $FQ2, "gzip -cd $fq2 | " or die $!;

open my $O1, " | gzip > $of1 " or die $!;
open my $O2, " | gzip > $of2 " or die $!;

my $seed = srand();
print "Seed: $seed\nFrac: $Frac\n";

while ( my $r1 = <$FQ1> ) {
    my $r2 = <$FQ2>;

    my $s1 = <$FQ1>;
    my $s2 = <$FQ2>;
    <$FQ1>;
    <$FQ2>;

    my $q1 = <$FQ1>;
    my $q2 = <$FQ2>;

    if (rand() < $Frac) {
        print $O1 "$r1$s1+\n$q1";
        print $O2 "$r2$s2+\n$q2";
    }
}
close $FQ1;
close $FQ2;

close $O1;
close $O2;

&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900, $t[4] + 1, @t[3,2,1,0], $_[0];
}

__END__

=head1 Function
    randomly select part of fastq reads

=head1 Usage
    perl $0 fq1 fq2 out1 out2

=head1 Options
    -f [f]  fraction of reads to be selected
    -h      help

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2020-1-10

=cut
