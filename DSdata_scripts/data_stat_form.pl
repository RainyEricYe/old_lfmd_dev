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

my ($smp_info, $inf, ) = @ARGV;

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

open my $IN, $inf or die $!;
while ( <$IN> ) {
    my ( $smp, @t ) = split;
    $h{$smp} = \@t;
}
close $IN;

$"="\t";
print "Sample\tClass\tLabel\tDementia\tStage\tTissue\tNread\tNsscs\tNdcs\tRatio\n";
for my $smp ( sort keys %h ) {
    print "$sm{$smp}\t@{$h{$smp}}\n";
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
