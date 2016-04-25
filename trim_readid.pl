#!/usr/bin/perl -w
use strict;

my $usage = "
This script will trim the name of read id

perl $0 <file_in> <file_out>

";

die $usage unless scalar @ARGV == 2;
my ( $filein, $fileout ) = @ARGV;

open(IN,"<$filein"); 
open(OUT,">$fileout"); 

my $cc=0; 

while(<IN>)
 {

  my $line=$_; 
  if($line=~/^\>/) { $cc++; print OUT ">read".$cc,"\n"; }
  else { print OUT $line; }

 }

close IN;
close OUT; 
