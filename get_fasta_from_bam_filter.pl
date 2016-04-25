#!/usr/bin/perl
use strict;
my $usage = '
perl $file_in $file_out
';
die $usage unless scalar @ARGV == 2;
my ( $file_in, $file_out) = @ARGV;
open(IN,"<$file_in");
open(OUT,">$file_out"); 
while(<IN>)
	{
		my $line=$_; 
		chomp($line); 
		my @ss=split("\t",$line);
		if(!($ss[2]=~/gi\|548558394/) && !($ss[2]=~/gi\|9626372/))
		{
		chomp($ss[0]); 
		chomp($ss[9]); 
		print OUT ">",$ss[0],"\n";
		print OUT $ss[9],"\n";}
	}
close IN; 
close OUT;
