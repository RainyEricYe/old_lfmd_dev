#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $Tools, $MinDcsDepth ) = (0, "LFMD", 10);
GetOptions(
    'help|?' => \$Help,
    'tools=s' => \$Tools,
    'minDcsDepth=i' => \$MinDcsDepth,
);
die `pod2text $0` if $Help or @ARGV < 2;

my ($list_f, $smp_info_f) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my %info;
open my $INF, $smp_info_f or die $!;
while (<$INF>) {
    chomp;
    my $smp = (split)[0];
    $info{$smp} = $_;
}
close $INF;

open my $LIST, $list_f or die $!;
my @smpV = ();
my %h;
while (my $f = <$LIST> ) {
    my $smp = (split /\./, basename($f))[0];
    push @smpV, $smp;

    open my $IN, $f or die $!;
    <$IN>;

    while ( <$IN> ) {
        my ( $chr, $ref, $pos, $depth, $mutN, $t, $c, $g, $a, $i, $d, $n, @info ) = split;
        next if $chr ne "MT2";

        $ref = uc $ref;
        @{ $h{$pos}{$smp} } = ($ref, $depth, $mutN, $t, $c, $g, $a );
    }
    close $IN;
}
close $LIST;

my @posV = ();
for my $pos ( sort {$a<=>$b} keys %h ) {
    push @posV, $pos if scalar(keys %{ $h{$pos} }) == scalar @smpV;
}
my @nt = qw/T C G A/;


$" = "\t";

my %smp_mut;
my $pos_count = 0;
for my $pos ( @posV ) {

    my $min_depth = 1000000;
    for my $smp ( @smpV ) {
        $min_depth = $h{$pos}{$smp}->[1] if defined $h{$pos}{$smp} && $h{$pos}{$smp}->[1] < $min_depth;
    }

    next if $min_depth < $MinDcsDepth;
    $pos_count++;

    for my $smp ( @smpV ) {
        if ( defined $h{$pos}{$smp} ) {
            my $r = $h{$pos}{$smp};

            if ( $r->[2]>0 && $r->[0] ne "N" ) { # mutN > 0
                for my $i ( 3..6 ) {
                    if ( $r->[$i] > 0 && $nt[$i-3] ne $r->[0] ) {
                        my $frac = $r->[$i] / $r->[1];
                        my $prob = (
                            $frac == 1
                            ? 1
                            : 1 - 10**( $min_depth*log10(1-$frac) )
                        );

                        my $mut = classify($r->[0], $nt[$i-3]);
                        $smp_mut{$smp}{$mut} += $prob;
                    }
                }
            }
        }
    }
}

my @mutV = qw/A:T>G:C A:T>C:G A:T>T:A C:G>T:A C:G>G:C C:G>A:T/;
print "Sample\tClass\tLabel\tDementia\tStage\tTissue\tTools\tMutation\tCount\n";

for my $smp ( @smpV ) {
    for my $mut ( @mutV ) {
        print "$info{$smp}\t$Tools\t$mut\t".(defined $smp_mut{$smp}{$mut} ? "$smp_mut{$smp}{$mut}\n" : "0\n" );
    }
}

&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub log10
{
    return log($_[0]) / log(10);
}

sub classify {
    my ( $ref, $alt ) = @_;
    my $mut = "$ref>$alt";

    if ( $mut eq "T>C" || $mut eq "A>G" ) {
        return "A:T>G:C";
    }
    elsif ( $mut eq "T>G" || $mut eq "A>C" ) {
        return "A:T>C:G";
    }
    elsif ( $mut eq "T>A" || $mut eq "A>T" ) {
        return "A:T>T:A";
    }
    elsif ( $mut eq "C>T" || $mut eq "G>A" ) {
        return "C:G>T:A";
    }
    elsif ( $mut eq "C>G" || $mut eq "G>C" ) {
        return "C:G>G:C";
    }
    elsif ( $mut eq "C>A" || $mut eq "G>T" ) {
        return "C:G>A:T";
    }
    else {
        warn "unknown mutation: $ref, $alt";
    }
}

print STDERR "total pos: ".(scalar keys %h)."\tpassed pos: $pos_count\n";

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
