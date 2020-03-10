#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $freqCut,)
=  ( 0,     1e-8,   );

GetOptions(
    'help|?' => \$Help,
    'freqCut=f' => \$freqCut,
);
die `pod2text $0` if $Help or @ARGV < 4;

my ($ina_f, $inb_f, $out, $suf ) = @ARGV;

$suf ||= "";
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#&showLog("START");

print "ina: $ina_f   inb: $inb_f \n\n";
open my $INA, $ina_f or die $!;
open my $INB, $inb_f or die $!;

my @nt = qw/T C G A ins del N/;
my @ref_nt = qw/T C G A N/;

my @Keys = ();
for my $r (@ref_nt) {
    for my $a (@nt) {
        next if $r eq $a;
        push @Keys, "$r>$a";
    }
}

my %ha;
my %hb;
my %supA;
my %supB;
<$INA>;
<$INB>;
while ( <$INA> ) {
    my ($c, $ref, $p, $dep, $mut, @t) = split;
    $ref =~ tr/acgtn/ACGTN/;

    for my $i (0..6) {
        if ( $t[$i] / $dep >= $freqCut ) {
            my $pt = $ref.">".$nt[$i];
            my $cp = "$c\t".($p-1)."\t$p";
            $ha{$pt}{$cp} = $t[$i] / $dep;
            $supA{$pt}{$cp} = $t[$i];
        }
    }
}
close $INA;

while ( <$INB> ) {
    my ($c, $ref, $p, $dep, $mut, @t) = split;
    $ref =~ tr/acgtn/ACGTN/;

    for my $i (0..6) {
        if ( $t[$i] / $dep >= $freqCut ) {
            my $pt = $ref.">".$nt[$i];
            my $cp = "$c\t".($p-1)."\t$p";
            $hb{$pt}{$cp} = $t[$i] / $dep;
            $supB{$pt}{$cp} = $t[$i];
        }
    }
}
close $INB;

my %ov;
my %fr;
my @mutA;
my @mutB;
my @mutO;
for my $pt (@Keys) {
    next if $pt =~ /\>N/;

    $ov{$pt}{"ov"} = 0;
    $ov{$pt}{"onlyA"} = 0;
    $ov{$pt}{"onlyB"} = 0;

    for my $cp ( sort keys %{$ha{$pt}} ) {
        if ( defined $hb{$pt}{$cp} ) {
            $ov{$pt}{"ov"}++;
            push @{ $fr{$pt}{"ov"} }, $ha{$pt}{$cp}."\t".$hb{$pt}{$cp};
            push @mutO, "$cp\t$pt\t$supA{$pt}{$cp}\t$ha{$pt}{$cp}\t$supB{$pt}{$cp}\t$hb{$pt}{$cp}";
        }
        else {
            $ov{$pt}{"onlyA"}++;
            push @{ $fr{$pt}{"onlyA"} }, $ha{$pt}{$cp};
            push @mutA, "$cp\t$pt\t$supA{$pt}{$cp}\t$ha{$pt}{$cp}";
        }
    }

    for my $cp ( sort keys %{$hb{$pt}} ) {
        if ( !defined $ha{$pt}{$cp} ) {
            $ov{$pt}{"onlyB"}++;
            push @{ $fr{$pt}{"onlyB"} }, $hb{$pt}{$cp};
            push @mutB, "$cp\t$pt\t$supB{$pt}{$cp}\t$hb{$pt}{$cp}";
        }
    }
}

open my $OA, ">$out.onlyA" or die $!;
open my $OB, ">$out.onlyB" or die $!;
open my $OV, ">$out.overlap" or die $!;
open my $OOM, ">$out.overlap.mut" or die $!;
open my $OAM, ">$out.onlyA.mut" or die $!;
open my $OBM, ">$out.onlyB.mut" or die $!;

$" = "\n";
my %tot;
for my $pt ( sort keys %ov ) {
    print qq($pt\t$ov{$pt}{"onlyA"}\t$ov{$pt}{"ov"}\t$ov{$pt}{"onlyB"}\n);
    defined $fr{$pt}{"onlyA"} && print $OA qq(@{ $fr{$pt}{"onlyA"}}\n);
    defined $fr{$pt}{"onlyB"} && print $OB qq(@{ $fr{$pt}{"onlyB"}}\n);
    defined $fr{$pt}{"ov"} && print $OV qq(@{ $fr{$pt}{"ov"} }\n);

    $tot{"onlyA"} += $ov{$pt}{"onlyA"};
    $tot{"onlyB"} += $ov{$pt}{"onlyB"};
    $tot{"ov"}    += $ov{$pt}{"ov"}   ;
}

print qq(Total\t$tot{"onlyA"}\t$tot{"ov"}\t$tot{"onlyB"}\n);

print $OAM "@mutA\n";
print $OBM "@mutB\n";
print $OOM "@mutO\n";
close $OA;
close $OB;
close $OV;
close $OAM;
close $OBM;

#&showLog("END");

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
