#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $minSizeOfRdf, )
=  ( 0,     12,            );
GetOptions(
    'help|?' => \$Help,
    'minSizeOfRdf=i' => \$minSizeOfRdf,
);
die `pod2text $0` if $Help or @ARGV < 3;

my ($inbam, $outpre, $chrRange ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my ($chr, $start, $end) = split /[:-]/, $chrRange;
die "please set start <= end for chromosome range" if $start > $end;

open my $IN, "<$inbam.rdf" or die $!;
open my $BAM, "samtools view $inbam |" or die $!;
open my $FQ1, "| gzip > $outpre\_1.fq.gz" or die $!;
open my $FQ2, "| gzip > $outpre\_2.fq.gz" or die $!;

my $fix = "TGNCT";
while ( <$IN> ) {
    my ($c, $s, $e, $n) = split;

    if ( $c lt $chr ) {
        map { my $a = <$BAM>; } 1 .. $n;
        next;
    }

    if ($c gt $chr || ( $c eq $chr && $s > $end) ) {
        last;
    }

    if ( $n >= $minSizeOfRdf ) {
        if ( $e < $start || $s > $end ) { # not overlap
            map { my $a = <$BAM>; } 1 .. $n;
        }
        else {
            for (1 .. $n/2) {
                my (%a, %b);
                my @t = qw/id flag seq qua/;
                @a{@t} = (split /\t/, <$BAM>, 12)[0,1,9,10];
                @b{@t} = (split /\t/, <$BAM>, 12)[0,1,9,10];

                if ( $a{flag} & 0x10 ) {
                    $a{seq} = reverse_complement( $a{seq} );
                    $a{qua} = reverse $a{qua};
                }

                if ( $b{flag} & 0x10 ) {
                    $b{seq} = reverse_complement( $b{seq} );
                    $b{qua} = reverse $b{qua};
                }

                my $tag = (split /\|/, $a{id})[-1];
                my $tag_a = substr($tag, 0, 12);
                my $tag_b = substr($tag, 12   );

                if ( $a{flag} & 0x40 ) {
                    $a{id} =~ s/\|$tag$/\/1/;
                    $b{id} =~ s/\|$tag$/\/2/;

                    print $FQ1 "\@$a{id}\n$tag_a$fix$a{seq}\n+\n".("B" x 17)."$a{qua}\n";
                    print $FQ2 "\@$b{id}\n$tag_b$fix$b{seq}\n+\n".("B" x 17)."$b{qua}\n";
                }
                else {
                    $b{id} =~ s/\|$tag$/\/1/;
                    $a{id} =~ s/\|$tag$/\/2/;

                    print $FQ1 "\@$b{id}\n$tag_a$fix$b{seq}\n+\n".("B" x 17)."$b{qua}\n";
                    print $FQ2 "\@$a{id}\n$tag_b$fix$a{seq}\n+\n".("B" x 17)."$a{qua}\n";
                }

            }
        }
    }
    else {
        map { my $a = <$BAM>; } 1 .. $n;
    }
}
close $IN;
close $BAM;
close $FQ1;
close $FQ2;


&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub reverse_complement {
    my @a = split //, $_[0];
    my $re;
    while ( my $nt = pop @a ) {
        $re .= (
            $nt eq "A" ? "T"
            : $nt eq "T" ? "A"
            : $nt eq "G" ? "C"
            : $nt eq "C" ? "G"
            : $nt
        );
    }

    return $re;
}

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
