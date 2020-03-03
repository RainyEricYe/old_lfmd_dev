#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, );
GetOptions(
    'help|?' => \$Help,
);
die `pod2text $0` if $Help or @ARGV < 1;


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my %n = (
    "T:A>C:G" => 0,
    "T:A>G:C" => 0,
    "T:A>A:T" => 0,
    "C:G>T:A" => 0,
    "C:G>G:C" => 0,
    "C:G>A:T" => 0,
);

my @nt = qw/T C G A/;
my %h;

for my $inf ( @ARGV ) {
    open my $IN, $inf or die $!;
    while ( <$IN> ) {
        next if /^Chrom/;
        my ( $chr, $ref, $pos, $depth, $mutN, @t ) = split;
        next if $chr ne "MT2";

        $ref = uc($ref);
        if ( $mutN > 0 ) {
            for my $i (0..3) {
                if ( $nt[$i] ne $ref && $t[$i] > 0 && $ref =~ /[ACGT]/) {
                    $h{"$ref>$nt[$i]"}++;
                }
            }
        }
    }
    close $IN;
}

for my $k ( keys %h ) {
    if ( $k eq "T>C" || $k eq "A>G" ) {
        $n{"T:A>C:G"} += $h{$k};
    }
    elsif ( $k eq "T>G" || $k eq "A>C" ) {
        $n{"T:A>G:C"} += $h{$k};
    }
    elsif ( $k eq "T>A" || $k eq "A>T" ) {
        $n{"T:A>A:T"} += $h{$k};
    }
    elsif ( $k eq "C>T" || $k eq "G>A" ) {
        $n{"C:G>T:A"} += $h{$k};
    }
    elsif ( $k eq "C>G" || $k eq "G>C" ) {
        $n{"C:G>G:C"} += $h{$k};
    }
    else {
        $n{"C:G>A:T"} += $h{$k};
    }
}

for my $k ( sort keys %n ) {
    print "$k\t$n{$k}\n";
}


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
