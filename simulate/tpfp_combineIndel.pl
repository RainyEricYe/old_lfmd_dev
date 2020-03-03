#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $debug);
GetOptions(
    'help|?' => \$Help,
    'debug' => \$debug,
);
die `pod2text $0` if $Help or @ARGV < 1;

my ($mut_list, $mutpos_f) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my %snv;
my %del;
my %ins;

open my $LIST, $mut_list or die $!;
while ( <$LIST> ) {
    chomp;
    my ( $c, $p, $ref, $mut, $type ) = split /\t/, $_;

    if ( $type eq "SNV" ) {
        $snv{$p}++;
    }
    else {
        if ( length $mut > length $ref ) { # ins
            $ins{$p}++;
        }
        else { # del
            $del{$p}++;
        }
    }
}
close $LIST;

my @tp = (0) x 2;
my @fp = @tp;

open my $IN, $mutpos_f or die $!;
<$IN>;
while ( <$IN> ) {
    my ( $c, $ref, $p, $depth, $mutN, $T,$C,$G,$A, $I,$D,$N) = split;

    if ( $mutN > 0 ) {
        if ( defined $snv{$p} ) {
            $tp[0]++;
        }
        else {
            $fp[0]++;
            print $_ if $debug;
        }
    }

    if ( $I > 0 ) {
        if ( defined $ins{$p} ) {
            $tp[1]++;
        }
        else {
            $fp[1]++;
            print $_ if $debug;
        }
    }

    if ( $D > 0 ) {
        if ( defined $del{$p} ) {
            $tp[1]++;
        }
        else {
            $fp[1]++;
            print $_ if $debug;
        }
    }
}
close $IN;

$"="\t";
print "@tp\t@fp\n";


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
