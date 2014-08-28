#!/usr/bin/perl
# filterGENEPOP.pl
# by Nate Campbell
# filter genepop file for selected loci.
# $ARGV[0] is a .txt file containing a single column of selected loci and $ARGV[1] is the genepop file you wish to filter.
# Some simple edits are required following filtering.  *Remove line ending at end of locus list *Insert line ending at EOF
# *Change line endings to windows format before importing into windows programs.

use strict; use warnings;

die "provide txt file of selected loci and genepop file\n" unless @ARGV == 2;

my @loci = ();
my %ArrayPos = ();

open (LIST, "<$ARGV[0]") or die "Cannot open locus list\n";

while (<LIST>) {
	chomp;
	push @loci, $_;
	}
close LIST;

foreach my $SNPs (@loci) {
	open (GENEPOP, "<$ARGV[1]") or die "Cannot open GENEPOP file\n";
		while (<GENEPOP>) {
			if ($. == 2) {my @markers = split ",", $_;
			my $length = @markers;
			for (my $i = 0; $i < $length; $i++) {
				if ($SNPs =~ m/$markers[$i]/) {$ArrayPos{$SNPs} = $i + 1}
					}
				}
		}
	close GENEPOP;
}

print ">$ARGV[0]\n";

foreach my $selected (sort keys %ArrayPos) {
	print "$selected\n";
	}

open (GENEPOP2, "<$ARGV[1]") or die;
	while (<GENEPOP2>) {
		chomp;
		#if ($_ =~ m/POP/) {print "\n$_"}
		if ($_ =~ m/,/) {
			print "\n";
			my @info = split "\t", $_;
			print "$info[0]";
			foreach my $genos (sort keys %ArrayPos) {
				print "\t$info[$ArrayPos{$genos}]";
				}
			}
	}
close GENEPOP2;
