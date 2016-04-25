#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample dir>
<sample dir> = full path of the folder holding files for this sample
               without last "/"
';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

my $finished = &check_RefG_split($dir);
#print $finished; 
exit ($finished);

###########################################################################
sub check_RefG_split {
	my ( $dir ) = @_;
	my $tot_RefG_seq=0;
	my $tot_seq = 0;

	opendir(DH, $dir) or return 10;
	foreach my $name (readdir DH) {
		if ($name =~/\.goodSeq$/) {
			my $RefGFile = $dir."/".$name;
		    open (IN, $RefGFile) or return 10;
		    while (my $line = <IN>){
				if ($line =~ ">") {
					$tot_RefG_seq++;
				}
			}
			close IN;
 		}
	
        if ($name =~ /\.goodSeq_RefGblast$/) { # Blast RefG directory
			my $full_path = $dir."/".$name;
			opendir(SubDH, $full_path) or return 10;
			foreach my $file (readdir SubDH) {
 				if ($file =~ /\.fa$/ && !($file=~/\.RefGfiltered\.fa/)) { 
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

#	print "total seq after spliting: $tot_RefG_seq\n";
#	print "total input seq in .goodSeq: $tot_seq\n"; 

	if($tot_RefG_seq==$tot_seq) { return 1; } 
	else { 
		opendir(DH, $dir) or return 10;
		foreach my $name (readdir DH) {
 			if ($name =~ /\.goodSeq_RefGblast$/) { 
				my $full_path = $dir."/".$name;
				opendir(SubDH, $full_path) or return 10;
				foreach my $file (readdir SubDH) {
					my $File = $full_path."/".$file;
#            	    unlink $File; 
				}
 				close SubDH; 
#           	`rmdir $full_path`;
			}
 		}
  		close DH;  
  		return 10; 
	}
}
