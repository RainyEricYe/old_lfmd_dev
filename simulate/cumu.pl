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

&showLog("START");

my %h;
for my $i (0..$#ARGV ) {
    open my $IN, $ARGV[$i] or die $!;
    while ( <$IN> ) {
        my ($frag, @t) = split /\s+/, $_;
        if ($i == 0) {
            $h{$frag} = \@t;
        }
        else {
            map { $h{$frag}[$_] += $t[$_] } 0..$#t;
        }
    }
    close $IN;
}

$" = "\t";
for my $k ( sort {$a<=>$b} keys %h ) {
    print "$k\t@{$h{$k}}\n";
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
