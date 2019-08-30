#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $type, );
GetOptions(
    'help|?' => \$Help,
    'type=s' => \$type,
);
die `pod2text $0` if $Help or @ARGV < 1;

my ($list, ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my %h;
my @labs;
open my $LIST, "<$list" or die $!;
while (my $f = <$LIST>) {
    chomp $f;
    my $lab = (split /\./, basename($f), 2)[-1];
    push @labs, $lab;

    open my $IN, $f or die $!;
    while (<$IN> ) {
        next if $type && !/$type/;

        my ( $smp, $n ) = (split /[\/\s]/, $_)[0,-1];
        $h{$smp}{$lab} = $n;
    }
    close $IN;
}
close $LIST;

$"="\t";
print "sample\t@labs\n";
for my $smp (sort keys %h) {
    print "$smp";
    for my $lab ( @labs ) {
        print "\t$h{$smp}{$lab}";
    }
    print "\n";
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
