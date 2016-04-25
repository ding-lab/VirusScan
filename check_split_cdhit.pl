#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample folder>
<sample folder> = full path of the folder holding files for a sample

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

my $finished = &check_split_output($dir);
#print $finished; 
exit ($finished);

###########################################################################
sub check_split_output {
	my ( $dir ) = @_;
	my $tot_cdhit_seq=0;
	my $tot_seq = 0;

	opendir(DH, $dir) or return 10;
	foreach my $name (readdir DH) {
		if ($name =~/\.cdhit_out$/) {
			my $cdFile = $dir."/".$name;
		    open (IN, $cdFile) or return 10;
		    while (my $line = <IN>){
				if ($line =~ ">") {
					$tot_cdhit_seq++;
				}
			}
			close IN;
		}
		if ($name =~ /\.cdhit_out_RepeatMasker$/) { # RepeatMasker directory
			my $full_path = $dir."/".$name;
			opendir(SubDH, $full_path) or return 10;
			foreach my $file (readdir SubDH) {
				if ($file =~ /\.fa$/) { 
					my $faFile = $full_path."/".$file;
					my $count = 0;
					open (IN, $faFile) or return 10;
					while (my $line = <IN>){
						if ($line =~ ">") {
							$count++;
						}
					}
					close IN;
					$tot_seq += $count;                    
				}
			}
			close SubDH; 
		}
	}
	close DH; 

#	print "$tot_cdhit_seq\n";
#	print "$tot_seq\n"; 

	if($tot_cdhit_seq==$tot_seq) { return 1; } 
	else { 
		opendir(DH, $dir) or return 10;
		foreach my $name (readdir DH) {
			if ($name =~ /\.cdhit_out_RepeatMasker$/) {
				my $full_path = $dir."/".$name;
				my $com1 = "rm -rf $full_path";
				my $com2 ="mkdir $full_path"; 
				# print "com is $com\n";
				system ( $com1 );
				system ( $com2 ); 
			}
  		}
  		close DH;  
  		return 10; 
	}
}

