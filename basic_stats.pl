#!/usr/bin/perl
# ------------------------------------------------------------------------------
# This program reads a fasta file, then computes and reports
# number of sequences or contigs, the average sequence length,
# the total number of bases, the minimum and the maximum sequence length, 
# and the %GC.
# AUTHOR: Bruno Gomez-Gil (modified from script by Joseph Fass <- Brad Sickler)
# LAST REVISED: Agst 2016
# USAGE: basic_stats.pl your_fasta_file
#------------------------------------------------------------------------------- 
# Read in sequences from one or more fasta file

my @data_files;
for(my $i = 0; $i < ($#ARGV + 1); $i++){
	$data_files[$i] = $ARGV[$i];
}
my $Id;
my %seq;
foreach my $file (@data_files){
	open(FASTA, $file) or die"Can't open file $file\n";
	while (<FASTA>) {
#         if (/^>(.+)\s/)  { $Id = $1; }  # overwrites if id up to 1st whitespace is non-unique in file
        if (/^>(.*)$/)  { $Id = $1; }
# 		if (/(\w+)/) {
		elsif (/^(\S+)$/) 	{ $seq{$Id} .= $1 if $Id; }
	}
	close (FASTA);
}
#---------------------------------------------------------------------------------
# Calculate basic fasta statistics

my $count = 0;
my $totalLength = 0;
my $gcCount = 0;
my @seqLengths;
my $min = 9999;
my $max = 0;
foreach my $id (keys %seq) {
	push @seqLengths, length($seq{$id}); # record length for N50 calc's
	$count++; # count the number of contigs or sequences
	$totalLength += length($seq{$id}); # determines the total number of bases
	$gcCount += ($seq{$id} =~ tr/gGcC/gGcC/); #determines %GC
	if (length($seq{$id}) < $min) {
		$min = length($seq{$id});
	}
	if (length($seq{$id}) > $max) {
		$max = length($seq{$id});
	}
}
#-------------------------------------------------------------------------------------
# Calculate N25, N50, and N75 and counts number of sequences

my $N25; my $N50; my $N75;
my $N25count=0; my $N50count=0; my $N75count=0;
my $frac_covered = $totalLength;
@seqLengths = reverse sort { $a <=> $b } @seqLengths;
$N25 = $seqLengths[0];
while ($frac_covered > $totalLength*3/4) {
	$N25 = shift(@seqLengths);
	$N25count++; $N50count++; $N75count++;
	$frac_covered -= $N25;
}
$N50 = $N25;
while ($frac_covered > $totalLength/2) {
	$N50 = shift(@seqLengths);
	$N50count++; $N75count++;
	$frac_covered -= $N50;
}
$N75 = $N50;
while ($frac_covered > $totalLength/4) {
	$N75 = shift(@seqLengths);
	$N75count++;
	$frac_covered -= $N75;
}
#-------------------------------------------------------------------------------------
# Add commas in large numbers
  # total length
sub commify {
    my $totalLength = reverse $_[0];
    $totalLength =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $totalLength;
}
  # Number of sequences
sub commify {
    my $count = reverse $_[0];
    $count =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $count;
}
#-------------------------------------------------------------------------------------
# Print out the results
use Term::ANSIColor;
print color 'cyan';
print "\n ------------- BASIC FASTA STATISTICS --------------------\n";
print color 'reset';
print "\nTotal number of sequences:";
print color 'bold';
print "\t$count\n";
print color 'reset';
print commify ("Total number of bases: \t\t$totalLength bp ");
#print "Total number of bases:\t\t";
printf '(%.2f', ($totalLength/1000000);
print " Mb)";
printf "\nAverage sequence length: \t\%6.2f\n", ($totalLength/$count);
print "Minimum sequence length: \t$min bp\n";
print "Maximum sequence length: \t$max bp\n";
print "N50:\t\t\t\t$N50 bp\n";
printf "GC %%:\t\t\t\t%.2f %%\n", ($gcCount/$totalLength * 100);
print "50% of total sequence length is contained in ";
print color 'bold';
print "$N50count";
print color 'reset';
print " sequences.\n";
#print "\n";
print color 'cyan';
print "----------------------------------------------------------\n";
print color 'reset';
# this is the end
