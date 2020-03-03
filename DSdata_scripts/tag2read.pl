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

my ($infq, $outfq ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my $tag;

open my $IN, "gzip -cd $infq |" or die $!;
open my $OUT, "|gzip > $outfq " or die $!;
while ( <$IN> ) {
    if ( $. % 4 == 1 ) {
        chomp;
        my @t = split /[|\/]/, $_;
        my $uid = $t[-2];
        my $rd = $t[-1];

        $tag = substr($uid, 0, length($uid)/2 ) if $rd eq "1";
        $tag = substr($uid, length($uid)/2    ) if $rd eq "2";

        my $name = join "", @t[0..$#t-2];
        print $OUT "$name/$rd\n";
    }
    elsif ( $. % 4 == 2 ) {
        print $OUT $tag."TGACT$_";
    }
    elsif ( $. % 4 == 0 ) {
        print $OUT "B" x 17;
        print $OUT $_;
    }
    else {
        print $OUT $_;
    }
}

close $IN;
close $OUT;





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
