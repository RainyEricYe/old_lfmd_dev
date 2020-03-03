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

my ($ref, $query_f ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

`makeblastdb -in $ref -dbtype nucl`;
`blastn -db $ref -out $query_f.tmp -query $query_f -outfmt 3`;

my $ref_seq;
my $query_seq;

open my $IN, "<$query_f.tmp " or die $!;
while ( my $q = <$IN> ) {
    if ( $q =~ /^Query_1/ ) {
        my $r = <$IN>;
        my $space = <$IN>;

        my ( $name, $start, $seq, $end ) = split /\s+/, $q;
        my ( $r_name, $r_start, $r_seq, $r_end ) = split /\s+/, $r;

        $ref_seq .= $r_seq;
        $query_seq .= $seq;

        warn "input error for $name $r_name in line $." if $name ne "Query_1" || $r_name ne "0";
    }
}
close $IN;

my @qs = split //, $query_seq;
my @rs = split //, $ref_seq;

my $ins = "";
my $del = "";
my $cumuInsLen = 0;

for my $i ( 0..$#rs ) {
    if ( $rs[$i] eq '.' ) {
        if ( $rs[$i-1] eq '-' ) { # ins

            my $ins_start = $i - length($ins);
            print "MT\t".($ins_start - $cumuInsLen)."\t$qs[$ins_start-1]\t$qs[$ins_start-1]$ins\tINS\n";

            $cumuInsLen += length $ins;
            $ins = "";
        }
        elsif ( $qs[$i-1] eq '-' ) { # del
            my $del_start = $i - length($del);
            print "MT\t".($del_start - $cumuInsLen)."\t$qs[$del_start-1]$del\t$qs[$del_start-1]\tDEL\n";

            $del = "";
        }
        else {
            next;
        }
    }
    elsif ( $rs[$i] eq '-' ) { #insertion
        $ins .= $qs[$i];
    }
    elsif ( $qs[$i] eq '-' ) { # deletion
        $del .= $rs[$i];
    }
    elsif ( $rs[$i] =~ /[acgtnACGTN]/ ) { # snv
        print "MT\t".($i + 1 - $cumuInsLen)."\t$rs[$i]\t$qs[$i]\tSNV\n";
    }
    else {
        warn "unknown $i-th base [$rs[$i]]";
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
