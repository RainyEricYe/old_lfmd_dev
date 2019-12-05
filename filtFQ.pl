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

my ($inf, ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

open my $IN, "<$inf" or die $!;

while ( <$IN> ) {
    chomp;

    if ( %. % 4 == 2 )  {
        my @t = split //, $_;
        for my $i (@t) {
            if ( $i =~ /[ACGTN]/ ) {
                print "$i";
            }
            else {
                print "N";
            }
        }
        print "\n";
    }
    else {
        print "$_\n";
    }
}
close $IN;

&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900, $t[4] + 1, @t[3,2,1,0], $_[0];
}

__END__

=head1 Function
    convert fa to fq.gz file with base quality 'J'

=head1 Usage
    perl $0 in.fa | gzip > out.fq.gz

=head1 Options
    -h      help

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2019-11-23

=cut
