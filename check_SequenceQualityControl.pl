#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample dir>
<sample dir> = full path of the folder holding files for a sample

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

my $finished = &check_QC_read_number($dir);
#print $finished; 
exit ($finished);

##########################################################################
sub check_QC_read_number {
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
		if ($name =~ /\.badSeq$/) { 
			my $full_path = $dir."/".$name;
			open(IN, $full_path) or return 10;
			while(<IN>){ 
				if($_=~/total unique seq = (\d+)/) {$tot_seq=$1;} 
			}
            close IN; 
		}
	}
	close DH; 

	print "total unique seq in CD-HIT output file: $tot_cdhit_seq\n";
	print "total unique sequence in QC output file: $tot_seq\n"; 

	if(abs($tot_cdhit_seq-$tot_seq)/$tot_cdhit_seq<=0.00001) { print "1","\n"; return 1; } 
	else {
		 
		return 10; 
	}
}

