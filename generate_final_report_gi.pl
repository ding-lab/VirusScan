
#!/usr/bin/perl
use strict;

my $usage = "
This script will read corresponding files in the given director and 
generate a report which contains SampleDescription, SequenceReport,
AssignmentSummary, InterestingReads.

perl $0 <run folder> <program version>
<run folder> = full path of the folder holding files for this sequence run

";
die $usage unless scalar @ARGV == 2;
my ( $dir, $version ) = @ARGV;

my @temp = split("/", $dir);
my $run_name = pop @temp;
my $outFile = $dir."/Analysis_Report_".$run_name.".tsv";
open (OUT, ">$outFile") or die "can not open file $outFile!\n";

my ($wkday,$month,$day,$time,$year) = split(/\s+/, localtime);
print OUT "Viruscan V${version}; Processing date: $day-$month-$year\n";

my $c = "**************************************************************************\n";
my $c2 = "#########################################################################\n\n";
print OUT $c;


print OUT "Summary:\n\n";
&generate_SampleDescription( $dir );
print OUT "End of Summary\n\n";
#print OUT $c ;

#print OUT "\n\nSequence Report\n\n";
#&generate_SequenceReport( $dir );
#print OUT "End of Sequence Report\n\n";
#print OUT $c ;

#print OUT "\n\nTaxonomy Assignment:\n\n";
#&generate_AssignmentSummary( $dir );
#print OUT "End of Assignment\n\n";
#print OUT $c ;

#print OUT "\n\nInteresting Reads\n\n";
#&generate_InterestingReads( $dir );
#print OUT "End of Interesting Reads\n\n";
#print "\n";

#print OUT "# Finished\n";

exit;

############################################################################
sub generate_SampleDescription {
	my ($dir) = @_;

	# sample name => num of total sequence in the sample
	my %total_seq = ();

	print OUT $dir,"\n";
	printf OUT "%10s\t", " ";
	printf OUT  "%5s\t%15s\t%10s\t%40s\n", "NoHumanRead", "PercentIDrange", "gi", "IdentifiedNoHuman";

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		if (!($name =~ /\./)) {
		# name is either file name or sample name (directory)
			my $full_path = $dir."/".$name;
			if (-d $full_path) { # is a directory, sample directory
				# get total number of sequences in the sample
				my $tempF = $full_path."/".$name.".fa";
				$total_seq{$name} = &count_num_of_seq($tempF);

				# print out report for this sample
				printf OUT "%30s\t%8d\n", $name,  $total_seq{$name};
				my $Summary_file = $full_path."/".$name.".gi.AssignmentSummary";
				if (-e $Summary_file) {
					open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
					foreach (1..17) {
						<IN>;
					}
					while (<IN>) {
						if ($_ =~ /^\s*$/) { # empty line
							next;
						}
						elsif ($_ =~ /# Finished Assignment Summary/) { 
							next;
						}
						else {
							chomp $_;
							my $number_reads = 0;
							my $range = "";
							my @temp = split(/\t/, $_);
							my $range = pop @temp;
							my $info = pop @temp;
							my $virus_info = pop @temp;
							my $virus = "";
                                                        my $gi=$temp[0]; 
							if ($info =~ /total number of reads: (\d+)/) {
								$number_reads = $1;
							}
	
							if ($virus_info =~ /hit does not have taxonomy entry/) {
								my @temp2 = split (",", $virus_info);
								$virus = shift @temp2;
							}
							else {
								my @temp2 = split(";", $virus_info);
								$virus = pop @temp2;
							}

						        if($number_reads>=1) 
							{
							printf OUT "%10s\t", " ";
							printf OUT "%5d\t%20s\t%10s\t%40s\n", $number_reads, $range, $gi, $virus; }
						}
					}
				}
				else {
					print OUT "$Summary_file does not exist!\n";
				}
			}
		}
	}
}

#####################################################################
# Assignment Summary
sub generate_AssignmentSummary {
	my ( $dir ) = @_;

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $full_path = $dir."/".$name;
		if (!($name =~ /\./)) {
			if (-d $full_path) { # is a directory
				my $Summary_file = $full_path."/".$name.".gi.AssignmentSummary";
				if (-e $Summary_file) {
					open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
					while (<IN>) {
						if ($_ =~ /# Finished Assignment Summary/) { 
							next;
						}

						print OUT $_;
					}
				}		
				print OUT $c2 ;
			}
		}
	}
}

##########################################################################
sub generate_SequenceReport {
	my ( $dir ) = @_;

	# sample name => num of total sequence in the sample
	my %total_seq = ();

	# sample name => num of unique sequence in the sample
	my %unique_seq = ();
	my %unique_seq_percent = ();

	# sample name => num of Filtered sequence in the libary
	my %bad_seq = ();

	# sample name => percentage of Filtered seq in the lib
	my %bad_percent = (); 

	# sample name => num of Filtered sequence in the libary
	my %lowComplex_seq = ();

	# sample name => percentage of Filtered seq in the lib
	my %lowComplex_percent = (); 

	# libary name => num of good sequenc in the sample
	my %good_seq = (); 

	# sample name => percentage of Filtered seq in the lib
	my %good_percent = ();

	# sample name => num of sequence assigned by BLASTN
	my %blastn_assigned = ();

	# sample name => percentage of sequences assigned by blastn
	my %blastn_assigned_percent = ();

	# sample name => num of sequence assigned by BLASTN
	my %blastx_assigned = ();

	# sample name => percentage of sequences assigned by blastn 
	my %blastx_assigned_percent = ();

	print OUT $dir,"\n";
	printf OUT "%30s\t", "sampleName";
	print OUT "total\tuniq\t\%\t Filtered\t\%\tLowComplex\t\%\tgood\t\%\tBNassign\t\%\tBXassign\t\%\n";
	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $full_path = $dir."/".$name;
		if (!($name =~ /\./)) {
			if (-d $full_path) { # is a directory
				# get total number of sequences in the sample
				my $tempF = $full_path."/".$name.".fa";
				$total_seq{$name} = &count_num_of_seq($tempF);

				# get number of unique sequence in the sample 
				$tempF = $full_path."/".$name.".fa.cdhit_out";
				if (-e $tempF) {
					$unique_seq{$name} = &count_num_of_seq($tempF);
					$unique_seq_percent{$name} = $unique_seq{$name}*100/$total_seq{$name};
					print "total # seq = ", $total_seq{$name}, " unique # seq: ", $unique_seq{$name}, "\n";
				}
				else {
					print OUT "$full_path does not have cdhit_out file!\n";
					return;
				}

				# get number of Filtered and good sequences
				##############################################################################
				# need to change here if seg masker enabled
				$tempF = $full_path."/".$name.".fa.cdhit_out.masked.badSeq";
				open (IN, $tempF) or die "can not open file $tempF!\n";
				while (<IN>) {
					if ($_ =~ /good seq = (\d+)/) {
#						print "num of good seq: $1, percentage: $2 (percentage of unique sequences\n";
						$good_seq{$name} = $1;
						$good_percent{$name} = $1*100/$total_seq{$name};
					}
					if ($_ =~ /bad seq = (\d+)/) {
#						print "num of Filtered seq: $1, percentage: $2 percentage of unique sequences\n";
						$bad_seq{$name} = $1;
						$bad_percent{$name} = $1*100/$total_seq{$name};
					}
					if ($_ =~ /Repeat and Low complexicity seq = (\d+)/) {
#						print "num of Filtered seq: $1, percentage: $2 percentage of unique sequences\n";
						$lowComplex_seq{$name} = $1;
						$lowComplex_percent{$name} = $1*100/$total_seq{$name};
					}
				}
				
				# get number of sequences assigned by BLASTn and number of sequences saved for BLASTX 
				my $total_saved = 0;
				my $total_BNassigned = 0;
				$tempF = $full_path."/".$name.".BNFiltered.fa";
				my $BNFiltered;
				if (-e $tempF) {
					$BNFiltered = &count_num_of_seq($tempF);
					$blastn_assigned{$name} = $good_seq{$name} - $BNFiltered;
					$blastn_assigned_percent{$name} = $blastn_assigned{$name}*100/$total_seq{$name};
				}
				else {
					$BNFiltered = 0;
					$blastn_assigned{$name} = $good_seq{$name} - $BNFiltered;
					$blastn_assigned_percent{$name} = $blastn_assigned{$name}*100/$total_seq{$name};
				}

				my $total_BXassigned = 0;
				$tempF = $full_path."/".$name.".gi.unassigned.fa";
				my $unassigned; 
				if (-e $tempF) {
					$unassigned = &count_num_of_seq($tempF);
				}
				else {
					$unassigned = 0;
				}
				$blastx_assigned{$name} = $BNFiltered - $unassigned;
				$blastx_assigned_percent{$name} = $blastx_assigned{$name}*100/$total_seq{$name};

				# print out report for this sample
				printf OUT "%30s\t%5d\t%5d\t%5.1f\t", $name, $total_seq{$name}, $unique_seq{$name}, $unique_seq_percent{$name};
				printf OUT "%5d\t%5.1f\t%5d\t%5.1f\t%5d\t%5.1f\t", $bad_seq{$name}, $bad_percent{$name}, $lowComplex_seq{$name}, $lowComplex_percent{$name}, $good_seq{$name}, $good_percent{$name};
				printf OUT "%5d\t%9.1f\t%5d\t%5.1f\n", $blastn_assigned{$name}, $blastn_assigned_percent{$name}, $blastx_assigned{$name},  $blastx_assigned_percent{$name};
			}
		}
	}


			# caclculate and print statistics for this run
			my $total = 0;
			my $unique = 0;
			my $bad = 0;
			my $good = 0;
			my $BNassign = 0;
			my $BXassign = 0;
			foreach my $key (keys %total_seq) {
				$total += $total_seq{$key};
				$unique += $unique_seq{$key};
				$bad += $bad_seq{$key};
				$good += $good_seq{$key};
				$BNassign += $blastn_assigned{$key};
				$BXassign += $blastx_assigned{$key};
			}
			$total_seq{"total"} = $total;
			$unique_seq{"total"} = $unique;
			$unique_seq_percent{"total"} = $unique*100/$total;
			$bad_seq{"total"} = $bad;
			$bad_percent{"total"} = $bad*100/$total;
			$lowComplex_seq{"total"} = $bad;
			$lowComplex_percent{"total"} = $bad*100/$total;
			$good_seq{"total"} = $good;
			$good_percent{"total"} = $good*100/$total;
			$blastn_assigned{"total"} = $BNassign;
			$blastn_assigned_percent{"total"} = $BNassign*100/$total;
			$blastx_assigned{"total"} = $BXassign;
			$blastx_assigned_percent{"total"} = $BXassign*100/$total;
				
			printf OUT "%30s\t%5d\t%5d\t%5.1f\t", "total", $total_seq{"total"}, $unique_seq{"total"}, $unique_seq_percent{"total"};
			printf OUT "%5d\t%5.1f\t%5d\t%5.1f\t%5d\t%5.1f\t", $bad_seq{"total"}, $bad_percent{"total"}, $lowComplex_seq{"total"}, $lowComplex_percent{"total"},  $good_seq{"total"}, $good_percent{"total"};
			printf OUT "%5d\t%9.1f\t%5d\t%5.1f\n", $blastn_assigned{"total"}, $blastn_assigned_percent{"total"}, $blastx_assigned{"total"},  $blastx_assigned_percent{"total"};

}

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

####################################################################################
# Assignment Summary
sub generate_InterestingReads {
	my ( $dir ) = @_;

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $full_path = $dir."/".$name;
		if (!($name =~ /\./)) {
			if (-d $full_path) { # is a directory
				print OUT $name, "\n";
				my $tempF = $full_path."/".$name.".gi.InterestingReads";
				if ( -e $tempF ) {
					open (IN, $tempF) or die "can not open file $tempF!\n";
					while (<IN>) {
						print OUT $_;
					}
					close IN;
				}
				else {
					print OUT "$name does not have .InteresingReads file!\n";
				}		
				print OUT $c2;
			}
		}
	}
}
