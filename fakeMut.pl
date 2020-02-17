#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $span, $indelSpan, )
=  ( 0,     200,   1000,       );
GetOptions(
    'help|?' => \$Help,
    'snvSpan=i' => \$span,
    'indelSpan=i' => \$indelSpan,
);
die `pod2text $0` if $Help or @ARGV < 2;

my ($ref, $nref) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

open my $IN, $ref or die $!;
open my $OUT, ">$nref" or die $!;

my @seq = ();
while ( <$IN> ) {
    chomp;
    if ( /^>/ ) {
        print $OUT ">fMT\n";
    }
    else {
        my @t = split //, $_;
        push @seq, @t;
    }
}
close $IN;

my %h = (
    A => "C",
    C => "G",
    G => "T",
    T => "A",
    N => "A",
);

my %indel = (
    A => "AC",
    C => "CGA",
    G => "GTCT",
    T => "",
    N => "",
);

my $o_pos = 69;
for ( my $i = 69; $i < 16569; $i += $span ) {
    if ( $i >= $o_pos + $indelSpan ) {
        print "MT\t".($i+1)."\t$seq[$i]\t$indel{ uc $seq[$i] }\tindel\n";
        $seq[$i] = $indel{ uc $seq[$i] };
        $o_pos += $indelSpan;
    }
    else {
        print "MT\t".($i+1)."\t$seq[$i]\t$h{ uc $seq[$i] }\tsnv\n";
        $seq[$i] = $h{ uc $seq[$i] };
    }
}

my $str = join "", @seq;
@seq = ();
@seq = split //, $str;

for ( my $j = 0; $j <= $#seq; $j += 50 ) {

    my $e = ($j + 49 > $#seq ? $#seq : ($j + 49) );
    print $OUT (join "", @seq[$j..$e])."\n";
}

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
