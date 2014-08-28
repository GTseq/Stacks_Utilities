#!/usr/bin/perl
# Stacks_Summary.pl
# by Nate Campbell
# Uses a list of STACKS tags and gathers consensus sequence from tags.tsv file and SNP position info from sumstats.tsv.
# Exports Summary data for SNP loci that can be used to gather paired end info for each locus later.
# works for RAD loci with up to 4 SNP sites.

use strict; use warnings;

die "usage: Provide a single column list of tags, a tags.tsv file, and a sumstats.tsv file" unless @ARGV == 3;

my @tags = ();

open (LIST, "<$ARGV[0]") or die "error reading $ARGV[0]\n";
	while (<LIST>) {
		if ($_ =~ m/^[0-9]/) 
			{
		chomp;
		push @tags, $_;
			}
		}

my %polys = ();
my %Positions = ();
my %SNP_ARR1 = ();
my %SNP_ARR2 = ();
my %SNP_ARR3 = ();
my %SNP_ARR4 = ();
my %SNP_Pos1 = ();
my %SNP_Pos2 = ();
my %SNP_Pos3 = ();
my %SNP_Pos4 = ();
my %SNP_Pos1_INS = ();
my %SNP_Pos2_INS = ();
my %SNP_Pos3_INS = ();
my %SNP_Pos4_INS = ();
my %ConsensusSEQ = ();
my %Search_String = ();
my %Read2_Search_String = ();

open (TAGS, "<$ARGV[1]") or die "error reading $ARGV[1]\n";
	while (<TAGS>) {
		chomp;
		my @info1 = split "\t", $_;
			foreach my $loci (@tags) {
				$polys{$loci} = 0;
				if ($info1[2] == $loci) {$ConsensusSEQ{$loci} = $info1[9]}
			}
		}

open (SUMSTATS, "<$ARGV[2]") or die "error reading $ARGV[2]\n";
	while (<SUMSTATS>) {
		chomp;
		if ($_ =~ m/^[0-9]/) {
			my @info2 = split "\t", $_;
				foreach my $loci2 (@tags) {
					if ($info2[1] == $loci2) {
						$polys{$loci2}++;
						if ($polys{$loci2} == 1) {
							$SNP_ARR1{$loci2} = $info2[4]; 
							$SNP_Pos1{$loci2} = $info2[4] + 1;
							$SNP_Pos1_INS{$loci2} = "[$info2[6]$info2[7]]";
							$Positions{$loci2} = "$SNP_Pos1{$loci2}";
							} 
						elsif ($polys{$loci2} == 2) {
							$SNP_ARR2{$loci2} = $info2[4]; 
							$SNP_Pos2{$loci2} = $info2[4] + 1;
							$SNP_Pos2_INS{$loci2} = "[$info2[6]$info2[7]]";
							$Positions{$loci2} = "$Positions{$loci2},$SNP_Pos2{$loci2}";
							}
						elsif ($polys{$loci2} == 3) {
							$SNP_ARR3{$loci2} = $info2[4]; 
							$SNP_Pos3{$loci2} = $info2[4] + 1;
							$SNP_Pos3_INS{$loci2} = "[$info2[6]$info2[7]]";
							$Positions{$loci2} = "$Positions{$loci2},$SNP_Pos3{$loci2}";
							}
						elsif ($polys{$loci2} == 4) {
							$SNP_ARR4{$loci2} = $info2[4]; 
							$SNP_Pos4{$loci2} = $info2[4] + 1;
							$SNP_Pos4_INS{$loci2} = "[$info2[6]$info2[7]]";
							$Positions{$loci2} = "$Positions{$loci2},$SNP_Pos4{$loci2}";
							}
					}
				}
			}
		}

foreach my $tagz (@tags) {
	my @bases = split "", $ConsensusSEQ{$tagz};
	splice @bases, $SNP_ARR1{$tagz}, 1, $SNP_Pos1_INS{$tagz};
	if (exists $SNP_Pos2_INS{$tagz}) {
	splice @bases, $SNP_ARR2{$tagz}, 1, $SNP_Pos2_INS{$tagz};
		}
	if (exists $SNP_Pos3_INS{$tagz}) {
	splice @bases, $SNP_ARR3{$tagz}, 1, $SNP_Pos3_INS{$tagz};
		}
	if (exists $SNP_Pos4_INS{$tagz}) {
	splice @bases, $SNP_ARR4{$tagz}, 1, $SNP_Pos4_INS{$tagz};
		}
	$Search_String{$tagz} = join ( '', @bases );
	splice @bases, 0, 69;
	$Read2_Search_String{$tagz} = join ( '', @bases );
	$Read2_Search_String{$tagz} = reverse $Read2_Search_String{$tagz};
	$Read2_Search_String{$tagz} =~ tr/ACGT][/TGCA[]/;
	}

#Print Headers...
print "TagID\tSNP-Position(s)\tTag-Summary\tRead2_Search-SEQ\n";

foreach my $loci3 (sort keys %ConsensusSEQ) {
	print "$loci3\t$Positions{$loci3}\t$Search_String{$loci3}\t$Read2_Search_String{$loci3}\n";
	}
