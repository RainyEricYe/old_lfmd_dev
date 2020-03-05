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
die `pod2text $0` if $Help or @ARGV < 2;

my ($smp_info, $inf_list, ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my %sm;
open my $SM, $smp_info or die $!;
while ( <$SM> ) {
    chomp;
    my $smp = (split)[0];
    $sm{$smp} = $_;
}
close $SM;


my %h;
my $max_p = 1;

open my $LIST, $inf_list or die $!;
while ( my $inf = <$LIST> ) {
    chomp $inf;

    my $smp = (split /\./, $inf)[0];

    open my $IN, $inf or die $!;
    while ( <$IN> ) {
        my ( $c, $p, $d ) = split;
        $h{$smp}{$p} = $d;

        $max_p = $p if $p > $max_p;
    }
    close $IN;
}
close $LIST;

print "Sample\tClass\tLabel\tDementia\tStage\tTissue\tPosition\tDepth\n";
for my $smp ( sort keys %h ) {
    for my $p ( 1..$max_p ) {
        print "$sm{$smp}\t$p\t";
        if ( defined $h{$smp}{$p} ) {
            print "$h{$smp}{$p}\n";
        }
        else {
            print "0\n";
        }
    }
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
