
#!/usr/bin/perl
use strict;

my $usage = '
This script will read the assignment report files in the given 
directory and generate a summary report for a given library. It will report 
in each library, for each category, how many total sequence were 
assigned to this category, how many were assigned by BLASTN, how many
were assigned by BLASTX.

It will also filter the virus lineage, leave out virus that are phage.
It will rank the virus lineage by range of percent ID from low to high. 

It will also generate a .InterestingReads report about the details of each lineage.

perl script <sample folder> <bad seq of repeatmasker (full path)>
<sample folder> = full path to the folder holding files for a given sample 

';

die $usage unless scalar @ARGV == 2;
my ( $dir, $bad_seq ) = @ARGV;

# cutoff for sequences to be interesting, we choose to report everything.
my $percentID_cutoff = 100;

my @temp = split("\/", $dir);
my $lib_name = pop @temp;

my $out = $dir."/".$lib_name.".gi.AssignmentSummary";
open (OUT, ">$out") or die "can not open file $out!\n";
my $out2 = $dir."/".$lib_name.".gi.InterestingReads";
open (OUT2, ">$out2") or die "can not open file $out2!\n";

my $seq_file = $dir."/".$lib_name.".fa";
my %sequences = &read_FASTA_data($seq_file); # read_ID => sequence

my %ID_low = ();  # lineage => lowest percent identity to hits
my %ID_high = (); # lineage => highest percent identity to hits
my %num_reads = (); 
my $C = "##############################################\n\n";

print OUT "$dir\n\n";
# get sequence statistics
my @nums = &get_SequenceInfo_OneSample($dir);
#print "@nums\n";

print OUT "#total\tuniq\ttotal\%\tFiltered\ttotal\%\tLowComplex\ttotal\%\tgood\ttotal\%\tBNRefG\ttotal\%\tBNNT\ttotal\%\tBXNR\ttotal\%\n";
printf OUT ("%d\t%d\t%5.1f\t%d\t%5.1f\t%d\t%5.1f\t%d\t%5.1f\t%d\t%5.1f\t%d\t%5.1f\t%d\t%5.1f\n", $nums[0], $nums[1], $nums[1]*100/$nums[0], $nums[2], $nums[2]*100/$nums[0], $nums[3], $nums[3]*100/$nums[0], $nums[4], $nums[4]*100/$nums[0], $nums[5], $nums[5]*100/$nums[0], $nums[6], $nums[6]*100/$nums[0], $nums[7], $nums[7]*100/$nums[0]);
print OUT "\n\n";

  
my $oldSeperator = $/;
$/ = "###########\n";
my $AssignmentReport_file = $dir."/".$lib_name.".gi.AssignmentReport";
open (IN, $AssignmentReport_file) or die "can not open file $AssignmentReport_file!\n";
my $line = <IN>;
$line =~ s/#//g;
my @temps = split("\n", $line);
shift @temps;
foreach my $temp (@temps) {
	print OUT $temp, "\n";
}
print OUT "\n\n";
		
while (<IN>) {
	if ($_ =~ /^\s*$/) { # skip blank line
		next;
	}
	elsif ($_ =~ /Finished Assignment Report/) { next; }

	my @lines = split("\n", $_);
	my $lineage = shift @lines;
	$lineage = shift @lines;
        #$lineage = shift @lines;
        #print $lineage,"\n";
	#<STDIN>; 
	my $high = 0;
	my $low = 100;
	my %readID_Identity = (); # readID => percent ID
	my %readID_desc = (); # readID => description of the read
	foreach my $l (@lines) {
		if ($l =~ /^\s*$/) { next; }
		elsif ($l =~ /QueryName/) { next; }
		elsif ($l =~ /reads from/) { next; }
		elsif ($l =~ /#+/) { next; }
		my ($read_ID, $Qlength, $hitName, $hitLen, $hitDesc, $alnLen, $ID, $hitS, $hitE, $e) = split("\t", $l);
		if($ID > $high) { $high = $ID;}
		if($ID < $low) { $low = $ID;}

		if (defined ($readID_Identity{$read_ID})) {
			if ($ID > $readID_Identity{$read_ID}) { 
				$readID_Identity{$read_ID} = $ID;
				$readID_desc{$read_ID} = $l;	
			}
		}
		else {
			$readID_Identity{$read_ID} = $ID;	
			$readID_desc{$read_ID} = $l;
		}
	}
	if ($high == 0) {
		$high = 100;
	}

	$ID_low{$lineage} = $low;
	$ID_high{$lineage} = $high;

        if($lineage=~/total number of reads: (\d+)/) { $num_reads{$lineage}=$1; }

	if ($low <= $percentID_cutoff) { 
		print OUT2 $lineage, "\t[$low, $high]\n\n";
		foreach my $key (sort {$readID_Identity{$a} <=> $readID_Identity{$b}} keys %readID_Identity) {
			print OUT2  $readID_desc{$key}, "\n";
		}
		print OUT2 "\n";
		foreach my $key (sort {$readID_Identity{$a} <=> $readID_Identity{$b}} keys %readID_Identity) {
			print OUT2 ">$key\n";
			print OUT2 "$sequences{$key}\n\n";
		}
	}
}
close IN;


foreach my $key (sort {$num_reads{$b} <=> $num_reads{$a}} keys %num_reads) {
	printf OUT  ("%s\t[%4.1f, %4.1f]%\n", $key, $ID_low{$key}, $ID_high{$key});
}
print OUT "# Finished Assignment Summary\n";

$/ = $oldSeperator;

close OUT;
close OUT2;

exit;

#####################################################################
sub read_FASTA_data () {
    my $fastaFile = shift @_;

    #keep old read seperator and set new read seperator to ">"
    my $oldseperator = $/;
    $/ = ">";
	 
    my %fastaSeq;	 
    open (FastaFile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";

    while (my $line = <FastaFile>){
		# Discard blank lines
        if ($line =~ /^\s*$/) {
			next;
		}
		# discard comment lines
		elsif ($line =~ /^\s*#/) {
	       next;
		}
		# discard the first line which only has ">", keep the rest
		elsif ($line ne ">") {
		    chomp $line;
		    my @rows = ();
		    @rows = split (/\n/m, $line);	
		    my $temp = shift @rows;
			my @temp_arr = split(/\s/, $temp);
			my $contigName = shift @temp_arr;
		    my $contigSeq = join("", @rows);
		    $contigSeq =~ s/\s//g; #remove white space
		    $fastaSeq{$contigName} = $contigSeq;
#			print " name = \\$contigName\\, seq  = \\$contigSeq\\\n\n";
		}
    }

    # check
#	 foreach my $key (keys %fastaSeq){
#	      print "Here is the key for fasta seq: $key \t $fastaSeq{$key}\n";
#	 }

    #reset the read seperator
    $/ = $oldseperator;
    
    return %fastaSeq;
}

 
##########################################################################
sub get_SequenceInfo_OneSample {
	my ( $dir ) = @_;

	my $total_seq = 0;
	my $unique_seq = 0;
	my $good_seq = 0;
	my $filtered_seq = 0;	
	my $RepeatLowComplex_seq = 0;
	my $blast_RefG_assigned = 0;
        my $blastn_assigned = 0;
	my $blastx_NR_assigned = 0;
#	my $tblastx_NTVS_assigned = 0;

	# get directory path
	my @fields = split(/\//, $dir);
	my $libName = $fields[$#fields];

	# get total number of sequences in the sample
	my $tempF = $dir."/".$libName.".fa";
	$total_seq = &count_num_of_seq($tempF);

	# get number of unique sequence in the sample 
	$tempF = $dir."/".$libName.".fa.cdhit_out";
	if (-e $tempF) {
		$unique_seq = &count_num_of_seq($tempF);
#		print "total # seq = ", $total_seq, " unique # seq: ", $unique_seq, "\n";
	}		
	
	# get number of Filtered and good sequences 
	#$tempF = $dir."/".$libName.".fa.cdhit_out.masked.badSeq";
	$tempF = $bad_seq; #cai changed, added segmasker
	if (-e $tempF) {
		open (IN, $tempF) or die "can not open file $tempF!\n";
	}	
	while (<IN>) {
		if ($_ =~ /good seq = (\d+)/) {
#			print "num of good seq: $1, percentage: $2 (percentage of unique sequences\n";
			$good_seq = $1;
		}
		if ($_ =~ /bad seq = (\d+)/) {
#			print "num of Filtered seq: $1, percentage: $2 percentage of unique sequences\n";
			$filtered_seq = $1;
		}
		if ($_ =~ /Repeat and Low complexicity seq = (\d+)/) {
#			print "num of Filtered seq: $1, percentage: $2 percentage of unique sequences\n";
			$RepeatLowComplex_seq = $1;
		}

	}
	close IN;

	# get number of sequences assigned by BLAST ReferenceGenome
	my $RefGfiltered = 0;
        my $tempF = $dir."/".$libName.".RefGfiltered.fa";
        if (-e $tempF) {
                $RefGfiltered = &count_num_of_seq($tempF);
        }
        else {
                $RefGfiltered = 0;
        }
        $blast_RefG_assigned = $good_seq - $RefGfiltered;

        # get number of sequences assigned by BLASTN
		my $BNFiltered = 0;
        my $tempF = $dir."/".$libName.".BNFiltered.fa";
        if (-e $tempF) {
                $BNFiltered = &count_num_of_seq($tempF);
        }
        else {
                $BNFiltered = 0;
        }
        $blastn_assigned = $RefGfiltered - $BNFiltered;

	# get number of sequences assigned by BLASTX NR  
	my $unassigned_num = 0;
	my $tempF = $dir."/".$libName.".unassigned.fa";
	if (-e $tempF) {
		$unassigned_num = &count_num_of_seq($tempF);
	}
	else {
		$unassigned_num = 0;
	}
	$blastx_NR_assigned = $BNFiltered - $unassigned_num;


	return ($total_seq, $unique_seq,  $filtered_seq, $RepeatLowComplex_seq, $good_seq, $blast_RefG_assigned, $blastn_assigned,  $blastx_NR_assigned);
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
