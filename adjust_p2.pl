#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $freq,)
=  ( 0,     1e-7, );
GetOptions(
    'help|?' => \$Help,
    'freq=f' => \$freq,
);
die `pod2text $0` if $Help or @ARGV < 1;

my ($inf, ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#&showLog("START");

my @ps;
open my $IN, "<$inf" or die $!;
LINE:
  while ( <$IN> ) {
    chomp;
    my @t = split;

    for my $x ( @t[12..$#t] ) {
        if ( $x =~ /^[ID]/ ) {
            next LINE;
        }
        else {
            my $p = (split /:/, $x)[1];
            next if $p == 0;
            push @ps, $p;
        }
    }
}
close $IN;

my $num = scalar @ps;
my %h;
my $rank = 1;
for my $k (sort {$a<=>$b} @ps) {
    push @{ $h{$k} }, $k * $num / $rank;
    $rank++;
}

$" = "\t";

my @nt = qw/T C G A/;

open my $RE, "<$inf" or die $!;
my $title = <$RE>;
print $title;

while ( <$RE> ) {
    chomp;
    my @t = split;

    if ( /[ID]/ && $t[4] > 0 ) {
        for my $i (0..3) {
            if ( $t[$i+5] > 0 && $nt[$i] ne $t[1] ) {
                $t[4] -= $t[$i+5];
                $t[$i+5] = 0;
            }
        }
    }

    my $out = "";
    my %indel;
    for my $x ( @t[12..$#t] ) {
        if ( $x =~ /^[ID]/ ) {
            my ($t, $num) = split /:/, $x;
            $indel{$t} = $num;
        }
        else {
            my ($nt, $p) = split /:/, $x;
            if ( !defined $p ) {
                next;
            }
            elsif ( $p == 0 ) {
                $out .= "\t$nt:$p";
            }
            elsif ( !/[ID]/ ) {
                my $adjust_p = pop @{$h{$p}};
                $out .= sprintf "\t$nt:%.2e", $adjust_p;
                $out .= "F" if $adjust_p > $freq;
            }
            else {
                1;
            }
        }
    }

    my @ks = sort { $indel{$b} <=> $indel{$a} } keys %indel;
    if ( scalar @ks > 1 ) {
        for my $k (@ks[1..$#ks]) {
            if ( $k =~ /^I/ ) {
                $t[9] -= $indel{ $k };
            }
            else {
                $t[10] -= $indel{ $k };
            }
        }
    }

    print "@t[0..11]$out";
    print "\t$ks[0]:$indel{ $ks[0] }" if defined $ks[0];
    print "\n";
}
close $RE;

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
