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
my $rid;
my $new_rid;
open my $IN, "samtools view $inf -h | " or die $!;
while ( <$IN> ) {
    if ( /^@/ ) {
        print $_;
    }
    else {
        if ( /UG:i:/ ) {

            my @t = split /\t/, $_, 2;
            $t[1] =~ /BX:Z:([ACGTN]+)/;
            my $umi = $1;
            $rid = $t[0];

            $t[0] =~ s/[ACGTN]+$/$umi/;
            $new_rid = $t[0];

            print "$t[0]\t$t[1]";
        }
        else {
            my @t = split /\t/, $_, 2;
            if ( $t[0] eq $rid ) {
                print "$new_rid\t$t[1]";
            }
        }
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
    update readID with grouped UMI

=head1 Usage
    perl $0

=head1 Options
    -h      help

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2019-11-23

=cut
