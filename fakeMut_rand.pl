#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $debug, $snvNum, $insNum, $delNum, $insMaxLen, $delMaxLen, $Seed, )
 = ( 0,     0,      200,     100,     100,     6,          5,        , 0,     );
GetOptions(
    'help|?' => \$Help,
    'debug'  => \$debug,
    'snvNum=i' => \$snvNum,
    'insNum=i' => \$insNum,
    'delNum=i' => \$delNum,
    'insMaxLen=i' => \$insMaxLen,
    'delMaxLen=i' => \$delMaxLen,
    'Seed=i' => \$Seed,
);
die `pod2text $0` if $Help or @ARGV < 3;

my ($ref, $nref, $mut_list) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

open my $IN, $ref or die $!;

my @seq = ();

while ( <$IN> ) {
    chomp;
    if ( /^>/ ) {
    }
    else {
        my @t = split //, $_; 
        map { $_ = uc $_ } @t;
        push @seq, @t;
    }
}
close $IN;

#generate random positions, the shortest distance between them is $indelMaxLen+2 bp
#
my @pos;
my %flank;
my $total_mut = $snvNum + $insNum + $delNum;
my $indelMaxLen = ($insMaxLen >= $delMaxLen ? $insMaxLen : $delMaxLen);
my $refLen = scalar @seq;

my $seed = ( $Seed != 0 ? srand($Seed) : srand() );
print "seed: $seed\n";

for (1..$total_mut) {
    my $p = int( rand($refLen) ); # a random position between 0 to $refLen - 1

    while (    $p < 50
            or $p > $refLen - 50
            or defined $flank{$p}
            or $p-$indelMaxLen-2 < 0
            or $p+$indelMaxLen+2 > $refLen-1
    ) { # pos not ok
        $p = int( rand($refLen) );
    }

    push @pos, $p;
    for my $i ($p-$indelMaxLen-2 .. $p+$indelMaxLen+2) {
        $flank{$i}++;
    }
}

# generate mutations
#
my @nt = qw/A C G T/;
my @nseq = @seq;

open my $LIST, ">$mut_list" or die $!;

for (1 .. $insNum) {
    my $p = pop @pos;

    my $insLen = int( rand($insMaxLen) ) + 1;
    my $ins = "";
    map { $ins .= $nt[int(rand(4))] } 1..$insLen;

    $nseq[$p] .= $ins;

    my $np = get_new_pos(\@seq, $p, $ins);
    print $LIST "MT\t$np\t$seq[$np-1]\t$seq[$np-1]$ins\tINS\n";

    if ( $debug && $p + 1 != $np ) {
        my $reg = "NC_012920.1:".($p - 4)."-".($p + 6);
        my $info = `samtools faidx $ref $reg`;
        print $LIST "~$p\n$info";
    }
}

for (1 .. $delNum) {
    my $p = pop @pos;

    my $delLen = int( rand($delMaxLen) ) + 1;
    my $del = "";
    map { $del .= $seq[$_] } $p+1 .. $p+$delLen;
    map { $nseq[$_] = "" } $p+1 .. $p+$delLen;

    my $np = get_new_pos(\@seq, $p, $del);
    print $LIST "MT\t$np\t$seq[$np-1]$del\t$seq[$np-1]\tDEL\n";
    print $LIST "~$p\n" if $debug && $p + 1 != $np;

    if ( $debug && $p + 1 != $np ) {
        my $reg = "NC_012920.1:".($p - 4)."-".($p + 6);
        my $info = `samtools faidx $ref $reg`;
        print $LIST "~$p\n$info";
    }
}

for (1 .. $snvNum) {
    my $p = pop @pos;

    my $base = $nt[int(rand(4))];
    while ( $base eq $seq[$p] ) {
        $base = $nt[int(rand(4))];
    }
    $nseq[$p] = $base;

    my $np = $p + 1;
    print $LIST "MT\t$np\t$seq[$p]\t$nseq[$p]\tSNV\n";
}

close $LIST;
# output
#
open my $OUT, ">$nref" or die $!;
print $OUT ">fMT\n";

my $str = join "", @nseq;
@seq = (); 
@seq = split //, $str;

for ( my $j = 0; $j <= $#seq; $j += 50 ) { 

    my $e = ($j + 49 > $#seq ? $#seq : ($j + 49) );
    print $OUT (join "", @seq[$j..$e])."\n";
}
close $OUT;


&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub get_repeat_unit {
    my ($str, $n) = @_;
    if ( length($str) == $n ) {
        return $str;
    }

    my @t = split //, $str;
    my @v;

    while (@t) {
        my $s = "";
        map { $s .= shift @t } 1..$n;
        push @v, $s;
    }

    for my $i (1..$#v) {
        return "" if $v[$i] ne $v[0];
    }

    return $v[0];
}

sub max_left_shift {
    my ( $v, $p, $u ) = @_; 

    my $n = length $u; 
    my $max_shift = 0;
    return 0 if $p < $n; 

    if ( $n == 1 ) { 
        while ( $max_shift < $p ) { 
            if ( $v->[$p-$max_shift] eq $u ) { 
                $max_shift++;
            } 
            else {
                return $max_shift;
            }
        }
    }
    else {
        do {
            my $str = join "", @$v[$p-$n-$max_shift+1..$p-$max_shift];
            if ( $str eq $u ) { 
                $max_shift += $n; 
            }
            else {
                return $max_shift;
            }
        } while ($p-$n-$max_shift+1 >= 0); 

    }

    return $max_shift;
}

sub get_new_pos {
    my ($r, $p, $str) = @_;
    my $len = length $str;

    my $offset = 0;
    for my $n ( 1..$len ) {
        if ( $len % $n == 0 ) {
            my $unit = get_repeat_unit($str, $n);
            if ( length $unit > 0 ) {
                my $max_shift = max_left_shift($r, $p, $unit);
                $offset = $max_shift if $max_shift > $offset;
            }
        }
    }

    return $p - $offset + 1;
}

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900, $t[4] + 1, @t[3,2,1,0], $_[0];
}

__END__

=head1 Function
    Generate mutation list for mitochondrial genome

=head1 Usage
    perl $0

=head1 Options
    -h      help

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2020-2-6

=cut
