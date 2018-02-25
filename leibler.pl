#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

die "usage: $0 <dna.gz> <ann.gz> <type> <kmer length>" unless @ARGV == 4;
my ($dna_file, $ann_file, $type, $k) = @ARGV;

# part 1: read the dna into memory
my %dna;
my $chr;
open(my $fh1, "gunzip -c $dna_file |") or die "error reading $ARGV[0]";
while (<$fh1>) {
	chomp;
	if (/^>(\S+)/) {
		$chr = $1;
	}
	else {
		$dna{$chr} .= $_;
	}
}

# part 2: get the sequence features
my @intron_seq;
open(my $fh2, "gunzip -c $ann_file |") or die "error reading $ARGV[1]";
while (<$fh2>) {
	chomp;
	my @f = split;
	if ($_ =~ /^#/) {
		next;
	}
	if ($f[2] eq $type) {
		my $chr = $f[0];
		my $beg = $f[3] -1;
		my $len = $f[4] - $f[3] + 1;
		next if $len <= 30 or $len >= 10000;
		my $str = $f[6];
		my $seq = substr($dna{$chr}, $beg, $len); 
		if ($str eq '-') {
			$seq = reverse $seq;
			$seq =~ tr/ACGT/TGCA/;                                            
		}
		push (@intron_seq, $seq);
	}
}

@intron_seq = sort {length $a <=> length $b} @intron_seq;
print join ("\n", @intron_seq);

my $quartile = @intron_seq / 4;
my @short_type = splice (@intron_seq, 0, $quartile);
my @long_type = splice (@intron_seq, -0, $quartile);

my %intron = (				#hash (%introns) with array references 
	shortA => [],
	shortB => [],
	longA => [],
	longB => [],
);

foreach my $short_type(@short_type) {								#loop the array 
	if (rand() < .5) {push @{$intron{shortA}}, $short_type;}			#fill the arrays in the hash %introns
	else			 {push @{$intron{shortB}}, $short_type;}		
}

foreach my $long_type(@long_type) {
	if (rand() < .5) {push @{$intron{longA}}, $long_type;}
	else			 {push @{$intron{longB}}, $long_type;}		
}

# compare data sets

my @cat = keys %intron;										#cat reads as "shortA, shortB, longA, longB" in an array

for (my $i = 0; $i < @cat; $i++) {							#compares two random categories
	for (my $j = $i + 1; $j < @cat; $j++) {				
		my $dist = kl_distance($intron{$cat[$i]}, $intron{$cat[$j]});		#input intron sequences from the shortA, shortB, etc.
		print "$cat[$i] vs $cat[$j] = $dist\n";
	}
}

# print Kullback Leibler distances for __ vs. __ 

my $shA_shB = kl_distance($intron{shortA}, $intron{shortB});
my $shA_lA = kl_distance($intron{shortA}, $intron{longA});
my $shA_lB = kl_distance($intron{shortA}, $intron{longB});
my $shB_lA = kl_distance($intron{shortB}, $intron{longA});
my $shB_lB = kl_distance($intron{shortB}, $intron{longB});
my $lA_lB = kl_distance($intron{longA}, $intron{longB});

print "shortA vs shortB: $shA_shB\n";
print "shortA vs longA: $shA_lA\n";
print "shortA vs longB: $shA_lB\n";
print "shortB vs longA: $shB_lA\n";
print "shortB vs longB: $shB_lB\n";
print "longA vs longB: $lA_lB\n";
permutations();

sub kl_distance {
	my ($introns1ref, $introns2ref) = @_;	
	my @introns1 = @{$introns1ref};											#reference to arrays for intron sequences (ex: short A sequences)
	my @introns2 = @{$introns2ref};
	
	my %P = %{kmer_freq($introns1ref)};										#%P: key = kmer, value = frequency
	my %Q = %{kmer_freq($introns2ref)};
	
	my $dist = 0;
	foreach my $kmer (keys(%P)) {											#iterate through the keys of the %P hash
		if (defined $Q{$kmer}) {
			$dist += $P{$kmer} * log ($P{$kmer} / $Q{$kmer});				#$P{$kmer} = the value of the %P hash --> frequency
		} else	{
			$Q{$kmer} = 0;
			$dist += $P{$kmer} * log ($P{$kmer} / $Q{$kmer});
		}
	}
	return $dist;
}

my $total;
sub kmer_freq {
	my ($intronsref) = @_;													#input intron SEQUENCES (value from hash %introns)
	my @intron_sequences = @{$intronsref};
	my %kmer_count;															#create a hash called kmer_count															#set the total count to 0	
	my %kmer_freq;																								
	for (my $i = 0; $i < @intron_sequences; $i++) {
		for (my $j = 0; $j < length($intron_sequences[$i]) - $k + 1; $j++) {			#loop through the sequences and fill the hash (key = kmer, value = count)
			my $kmer;
			$kmer = substr($intron_sequences[$i], $j, $k);
			$kmer_count{$kmer}++;
			$total++;
		}
	}
	foreach my $kmer (keys %kmer_count) {
		$kmer_freq{$kmer} = $kmer_count{$kmer} / $total;
	}
	return \%kmer_freq;
}

sub permutation_count {
	my $permutations = 4 ** $k;
	print "$permutations permutations\n";
	
}
