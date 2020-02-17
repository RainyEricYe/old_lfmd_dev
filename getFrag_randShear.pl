#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;
use Math::Random;

my ( $Help, $nFrag, $Seed, $insertSize, $sd, $Phrase, )
 = ( 0 ,    100,    0,     600,         20,  "hello", );
GetOptions(
    'help|?' => \$Help,
    'n=i' => \$nFrag,
    'Seed=i' => \$Seed,
    'Phrase=s' => \$Phrase,
);
die `pod2text $0` if $Help or @ARGV < 2;

my ($ref, $pre) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

open my $IN, "<$ref.fai" or die $!;
my ($chr, $length) = split /\s+/, <$IN>;
close $IN;
die "cannot get chr and length from $ref.fai" if !defined $chr or !defined $length;

my $seed = ( $Seed == 0 ? srand() : srand($Seed) );
print "Seed: $seed\n";

open my $REG, ">$pre.tmp.reg" or die $!;
my @starts;
for (1..$nFrag) {
    my $start_max = $length - $insertSize - $sd * 4;
    my $start = int(rand($start_max));
    push @starts, $start;
}

random_set_seed_from_phrase($Phrase);
print "Phrase: $Phrase\n";

my @sizes = random_normal($nFrag, $insertSize, $sd);

for my $i ( 0..$#starts ) {
    my $s = $starts[$i];
    my $e = $s + int($sizes[$i]) - 1;
    print $REG "$chr:$s-$e\n";
}

close $REG;

my $cnt = 1;
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
    perl $0 generate randome sheared fragment from mtDNA

=head1 Options
    -h      help

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2020-02-11

=cut
