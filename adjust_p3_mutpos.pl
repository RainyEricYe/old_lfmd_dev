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

my %h;
open my $IN, "<$inf" or die $!;
my $title = <$IN>;
print $title;

while ( <$IN> ) {
    chomp;
    my @t = split;

    if ( $t[9] > 0 || $t[10] > 0 ) {
        $h{ $t[0] }{ $t[2] } = 1;
    }
}
close $IN;

my @nt = qw/T C G A/;

open my $RE, "<$inf" or die $!;
<$RE>;
while ( <$RE> ) {
    chomp;
    my @t = split;
    my $nearby = 0;

    if ( defined $h{ $t[0] } ) {
        for my $k ( keys %{ $h{ $t[0] } } ) {
            if (  ($k - $t[2])**2 <= 25 ) {
                $nearby++;

                for my $i (0..3) {
                    if ( $t[4] > 0 && $t[$i+5] > 0 && $nt[$i] ne $t[1] ) {
                        $t[4] -= $t[$i+5];
                        $t[$i+5] = 0;
                    }
                }

            }
        }
    }

    if ( $nearby ) {
        print "@t\n";
    }
    else {
        print "$_\n";
    }

}
close $RE;

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
