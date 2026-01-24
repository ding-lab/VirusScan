#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample dir>
<sample dir> = full path of the folder holding files for this sample
               without last "/"
';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

my $finished = &check_split_output($dir);
#print $finished;
exit ($finished);

##############################################################
sub check_split_output {
	my ( $dir ) = @_;
	my $tot_BN_seq=0;
	my $tot_seq = 0;

	opendir(DH, $dir) or return 10;
	foreach my $name (readdir DH) {
		if ($name =~/\.RefGfiltered\.fa$/) {
			my $RefGFile = $dir."/".$name;
			open (IN, $RefGFile) or return 10;
			while (my $line = <IN>){
				if ($line =~ ">") {
					$tot_BN_seq++;
				}
			}
			close IN;
		}
		if ($name =~ /\.RefGfiltered_BLASTN$/) { # BLASTN directory
			my $full_path = $dir."/".$name;
			opendir(SubDH, $full_path) or return 10;
			foreach my $file (readdir SubDH) {
				if ($file =~ /\.fa$/ && !($file=~/\.BNfiltered\.fa/)) { 
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

#	print "$tot_BN_seq\n";
#	print "$tot_seq\n"; 

	if($tot_BN_seq==$tot_seq) { return 1; } 
	else { 
		opendir(DH, $dir) or return 10;
		foreach my $name (readdir DH) {
#   	    print "$name\n";
			if ($name =~ /\.RefGfiltered_BLASTN$/) { 
				my $full_path = $dir."/".$name;
				opendir(SubDH, $full_path) or return 10;
				foreach my $file (readdir SubDH) {	
					my $faFile = $full_path."/".$file;
    	            # print "$faFile\n"; 
					unlink $faFile; 
				}
				close SubDH; 
			}
		}
	  	close DH;  
	  	return 10; 
	}
}
