#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $trimEnd, $min_support_read_N, $Debug, )
 = ( 0,     5,        3,                   1,      );
GetOptions(
	'help|?' => \$Help,
	'trimEnd=i' => \$trimEnd,
	'min_support_read_N=i' => \$min_support_read_N,
);
die `pod2text $0` if $Help or @ARGV < 3;

my ($mut_f, $dcs_bam_f, $raw_bam_f ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my %factorial;
$factorial{0} = 1;
for my $k (1..170) {
	$factorial{$k} = $k * $factorial{$k-1};
}

open my $MUT, $mut_f or die $!;
<$MUT>;

while ( <$MUT> ) {
	my ( $chr, $ref, $pos, $depth, $mutN, $t, $c, $g, $a, $i, $d, $n, @info ) = split;

	next if $chr ne "MT2" && $chr ne "NC_012920.1" && $chr ne "chrM";
	next if $mutN == 0 or $pos < 1000 or $pos > 15569;

	$ref = uc $ref;

	# estimate number of fragments which cover this position

	my %Nf; # number of fragments,  keys => insertSize, pos
	my %Nf_total; # number of total fragments for a specific insertSize

	my $Nb; # number of bins
	my %Bin;
	my %Bin_occupied; # occupied bins for a specific insertSize
	my %pr;
	my %Size;
	my %Lam;

#	print "$chr, $ref, $pos, $depth, $mutN, $t, $c, $g, $a, $i, $d, $n, \n";

	open my $RAW, "samtools view $raw_bam_f $chr:$pos-$pos |";
	while ( my $read = <$RAW> ) {
		my ($id, $flag, $c, $p, $mapQ, $cg, $mc, $mp, $insertSize, $seq) = split /\s+/, $read, 11;

		next if ( $insertSize <= 0
				or $mc ne "="
				or $flag & 0x4 or $flag & 0x8
				or $flag & 0x100 or $flag & 0x200 or $flag & 0x400 or $flag & 0x800 );

		# skip trimEnd
		next if $p + aln_length($cg) - $trimEnd <= $pos;
		next if $p + $trimEnd > $pos;

		$Nb = length($seq) * 2 - $trimEnd * 4 if !defined $Nb;

		if ( $flag & 0x40 ) { # watson family size++
			$Bin{$insertSize}{$p}[0]++;
		}
		else { #crick family size++
			$Bin{$insertSize}{$p}[1]++;
		}

		$Bin{$insertSize}{$p}[2]++; # family size++
		$Size{$insertSize}++;
	}
	close $RAW;


	for my $is (keys %Bin) {
		my $empty = ( $Nb - scalar keys %{ $Bin{$is} } );
		my $one_read = 0;
		map { $one_read++ if $Bin{$is}{$_}[2] == 1 } keys %{ $Bin{$is} }; # count read family only has one read

		my $o = $one_read / $Nb;
		my $e = $empty / $Nb;
		my $s = $Size{$is} / $Nb;
		next if $Size{$is} < $min_support_read_N * 2;

		# SUM{ exp(-k*$s/lamda - lamda)*lamda^k/factorial(k)
		#
		print "~~~~ Nb $Nb empty $empty size $Size{$is} $s insertSize $is " if $Debug;

		my $lamda = fit_lamda($s, $e);
		my $old_lamda =  log($Nb) - log($Nb-$empty);

		my $lamda2 = $s / $lamda;
		print " lamda $lamda lamda2 $lamda2\n" if $Debug;
		$Lam{$is}[0] = $lamda;
		$Lam{$is}[1] = $lamda2;
	}

	my $mutFreq = $mutN/$depth;
	$mutFreq = 0.0007;
	my $prob_notDetect = 0;

	open my $DCS, "samtools view $dcs_bam_f $chr:$pos-$pos |" or die $!;
	while ( my $read = <$DCS> ) {
		my ($id, $flag, $c, $p, $mapQ, $cg, $mc, $mp, $insertSize, $seq, $qua, $info) = split /\s+/, $read,12;

		$insertSize *= -1 if $insertSize <= 0;
		print "insertSize $insertSize $.\n";

		next if !defined $Lam{$insertSize};

		if ( $info =~ /sp:Z:[2-9]/ ) {
			next;
		}

		my $prob_notDetect_thisFam = 0;

		if ($info =~ /fs:Z:(\d+),(\d+)/ ) {
			my ( $w_s, $c_s ) =  ($1, $2);
			my $s = $w_s + $c_s;
			if ( $s > 170 ) {
				warn "skip dcs due to large family size";
				next
			}

			my $sum = 0;
			my $target_k = int($s/$Lam{$insertSize}[1]);
			my $target_k_start = ( $target_k-30 <= 1 ? 1 : $target_k-30 );

			for my $k ($target_k_start .. $target_k+30) {
				$sum += Prob_with_fragN_FamSize($k, $s, $Lam{$insertSize}[0], $Lam{$insertSize}[1]);
			}

			for my $k ($target_k_start .. $target_k+30) { # possible fragN in this read family

				my $prob = Prob_with_fragN_FamSize($k, $s, $Lam{$insertSize}[0], $Lam{$insertSize}[1])/$sum;

				next if $prob < 1e-6;
				print "$k, $s, $Lam{$insertSize}[0], $Lam{$insertSize}[1] $prob\n";

				my $w_notDetect = notDetect(1/$k, $w_s);
				my $c_notDetect = notDetect(1/$k, $c_s);

				$prob_notDetect_thisFam += $prob * ($w_notDetect + $c_notDetect - $w_notDetect * $c_notDetect);
				print "~~$w_notDetect $c_notDetect $prob_notDetect_thisFam\n";

			}
		}
		else {
			die "NO fs:Z tag in the bam file: $dcs_bam_f";
		}

		$prob_notDetect += log(1 - $mutFreq + $mutFreq * $prob_notDetect_thisFam);
	}
	close $DCS;

	my $sens = 1 - exp($prob_notDetect);

	print "freq:sensitivity => $mutFreq\t$sens\n";

}

close $MUT;




&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
sub notDetect {
	my ( $f, $s ) = @_;

	return 0 if $f == 1;

	return (
		(1-$f)**$s
		+ $s * $f * ( (1-$f)**($s-1) )
		+ 0.5 * $s * ($s-1) * ($f**2) * ( (1-$f)**($s-2) )
		+ 0.005 * $s * ($s-1) * ($s-2) * ($f**3) * ( (1-$f)**($s-3) )
	);
}
#
sub Prob_with_fragN_FamSize {
	my ($k, $j, $lamda, $lamda2) = @_;
	#return (($k*$lamda2)**$j) * ($lamda**$k) * exp(-$k*$lamda2-$lamda) / ( $factorial{$k} * $factorial{$j} );

	return exp( $j*log($k*$lamda2) + $k*log($lamda) + (-$k*$lamda2-$lamda) - log($factorial{$k}) - log($factorial{$j} ) );
}

# SUM{ exp(-k*$s/lamda - lamda)*lamda^k/factorial(k) } - $e
#
sub diff_func {
	my ($lam, $s, $e, $max_k) = @_;

	my $r = -$e;
	for my $k (0..$max_k) {
		my $temp = exp(-$k*$s/$lam - $lam) * ($lam**$k)/$factorial{$k};
		#print "~~ $k\t$temp\n";
		$r += $temp;
		last if $temp == 0 ;
	}

	return $r;
}

sub fit_lamda {
	my ($s,$e) = @_;

	my $min = 1;

	my $best_lam = 1;

	my $step = 0.001;
	for my $x (1..10000) {
		my $lam = $x * $step;
		my $max_k = 50;
		my $f = diff_func($lam, $s, $e, $max_k);

		#	print "$lam\t$f\n";
		if ($f**2 < $min) {
			$min = $f**2;
			$best_lam = $lam;
		}
		else {
			last;
		}
	}

	return $best_lam;
}

sub aln_length {

	my @f = split /\D+/, $_[0]; # cigar field length
	my @t = split /\d+/, $_[0]; # cigar type
	shift @t;

	warn "weird cigar $_[0] " if scalar @f != scalar @t;

	my $len = 0;
	for my $i (0..$#t) {
		if ( $t[$i] =~ /[MDN=]/ ) {
			$len += $f[$i];
		}
	}
	return $len;
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
	-h	  help

=head1 Author
	yerui;	yerui@genomics.cn

=head1 Version
	v1.0;	2012-

=cut
