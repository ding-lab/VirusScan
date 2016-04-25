#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts, "i:s", "o:s", "n=i", "p=s", "h");

die "usage:\t$0 <-i fasta_file> <-o out_put_dir> <-n number_of_file> <-p prefix_of_seq>[-h]\n" if(defined($opts{h}) || !defined($opts{i}) || !defined($opts{o}) || !defined($opts{p}));

# get total number of sequence in input file
my $count_seq= &count_num_of_seq($opts{i}); 

# calculate how many sequences in each file
my $size = $count_seq/$opts{n};
#print "$count_seq $size\n";

# start spliting
my $count = 1;
my $count_seq_each=0;
open(OUT, ">$opts{o}/$opts{p}${count}".".fa")||die $!;

open(SEQ, $opts{i}) || die "cannot open file : $opts{i}\n";
$/='>'; 
<SEQ>;
while(<SEQ>) {
	chomp;
	if($count_seq_each > $size) {
		close OUT;
		$count++; $count_seq_each=0; 
		open(OUT, ">$opts{o}/$opts{p}${count}".".fa")||die $!;
	}
	print OUT $/.$_; 

	$count_seq_each++; 
}
close OUT;
close SEQ;

############################################################################
sub count_num_of_seq () {
	my ($fastaFile) = @_;
	my $count = 0;
	
	open (FastaFile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";
	while (my $line = <FastaFile>){
		if ($line =~ ">") {
			$count++;
		}
	}
	close FastaFile;

	return $count;
}

