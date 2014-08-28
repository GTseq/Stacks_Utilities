#!/usr/bin/perl
# Short_Scaffold_PE-fetch.pl
# by Nate
# Use results of Stacks_Summary.pl to fetch R1 and R2 sequences that overlap to generate scaffolds for targeted RAD loci.
# Useful for generating scaffolds from 250-500 bases in length for specific RAD loci.
# Generates a single .fastq file for each RAD locus in the provided stacks_summary file.
# Usage: Provide <Stacks_summary.txt> <R1.fastq> <R2.fastq>

use strict; use warnings;

die "Usage: Provide <Stacks_summary.txt> <R1.fastq> <R2.fastq>\n" unless @ARGV == 3;

my %R1_SEQ = ();
my %R2_SEQ = ();

#Read-in target sequences...
open (TARGETS, "<$ARGV[0]") or die "Error opening $ARGV[0]\n";

while (<TARGETS>) {
	chomp;
	if ($_ =~ m/^[0-9]/) {
	my @info = split "\t", $_;
	$R1_SEQ{$info[0]} = $info[2];
	$R2_SEQ{$info[0]} = $info[3];
	}
}
close TARGETS;

#Match R1 Target sequences to R1 fastq data...

foreach my $tagz (sort keys %R1_SEQ) {
	my $R1_hits = 0;
	my %Coordinate = ();
	my $MarkSEQ1 = $R1_SEQ{$tagz};
	my $MarkSEQ2 = $R1_SEQ{$tagz};
	$MarkSEQ1 =~ s/[ACGT]\]//g;
	$MarkSEQ1 =~ s/\[//g;
	$MarkSEQ2 =~ s/\[[ACGT]//g;
	$MarkSEQ2 =~ s/\]//g;
	my $MarkerLength = length($MarkSEQ1);

	open (R1, "<$ARGV[1]") or die "Error opening $ARGV[1]\n";
	open(OUT, ">$tagz.fastq")or die;

	print "TagID\t$tagz\t$R1_SEQ{$tagz}\n";
	print OUT "\@HISEQ>$tagz-Allele1\n$MarkSEQ1\n+\n";
	for (my $i = 0; $i < $MarkerLength; $i++) {
		print OUT "J";
		}
	print OUT "\n";
	print OUT "\@HISEQ>$tagz-Allele2\n$MarkSEQ2\n+\n";
	for (my $j = 0; $j < $MarkerLength; $j++) {
	print OUT "J";
		}
	print OUT "\n";

while (<R1>) {
	my @CoordinateINFO = ();
	my $R1_ID1 = $_;
	my $R1_seq = <R1>;
	my $R1_ID2 = <R1>;
	my $R1_qual = <R1>;
	chomp($R1_seq);
	chomp($R1_qual);
	my $read_length = length($R1_seq) - 6;
	my $trimmed_R1_seq = substr($R1_seq,6,$read_length);
	my $trimmed_R1_qual = substr($R1_qual,6,$read_length);

	if ($R1_seq =~ m/$R1_SEQ{$tagz}/) {
		@CoordinateINFO = split(/ /, $R1_ID1);
		$Coordinate{$CoordinateINFO[0]} = $CoordinateINFO[0];
		$R1_hits++;
		if ($R1_hits < 40) 
			{
		print OUT "\@HISEQ>$tagz-F-$CoordinateINFO[0]\n$trimmed_R1_seq\n+\n$trimmed_R1_qual\n";
			}
		}
	}
close R1;
print "Hits in read 1 = $R1_hits\nMax R1 reads written = 40\n...Done with Read 1...\n";

#Use read 1 matches to gather read 2 sequences...

open(R2, "<$ARGV[2]") or die "Error opening $ARGV[2]\n";
print "Reading R2 file...\n";
my $R2_matches = 0;
my $R2_NoMatches = 0;

	while (<R2>) {
		my @CoordinateINFO2 = ();
		my $R2_ID1 = $_;
		my $R2_seq = <R2>;
		my $R2_ID2 = <R2>;
		my $R2_qual = <R2>;
		chomp ($R2_seq);
		chomp ($R2_qual);
		@CoordinateINFO2 = split(/ /, $R2_ID1);

		if ((exists $Coordinate{$CoordinateINFO2[0]}) && ($R2_seq =~ m/$R2_SEQ{$tagz}/) && ($R2_matches < 30))
			{
		$R2_matches++;
		print OUT "\@HISEQ>$tagz-R-$CoordinateINFO2[0]\n$R2_seq\n+\n$R2_qual\n";
			}
		elsif ((exists $Coordinate{$CoordinateINFO2[0]}) && ($R2_NoMatches < 30))
			{
		$R2_NoMatches++;
		print OUT "\@HISEQ>$tagz-R-$CoordinateINFO2[0]\n$R2_seq\n+\n$R2_qual\n";
			}		
		}
	print "R2 output reads spanning gap = $R2_matches\nR2 output reads not spanning gap = $R2_NoMatches\n";
	close R2;
	close OUT;
}

