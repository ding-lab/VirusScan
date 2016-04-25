
#!/usr/bin/perl -w

use strict;
use Bio::SearchIO;

my $Usage = '
This script accepts a BLASTn output file that were blasted against Reference 
genome, find out whether the best hit has a e value lower than the cutoff. If 
yes, output query information. If no, the sequence will be kept for further analysis.

perl script <dir><blast output file><Ref genome Taxonomy> 
<dir> = directory that blast output file resides in, without last "/"
<blast output file> = name of the blastn output file
<Ref genome Taxonomy> = Bacteria, Homo, Phage, Fungi, Mus, other

';

die $Usage unless scalar @ARGV == 3;
my ($dir, $blastout, $RefGenomeTaxonomy) = @ARGV;

# cutoff value for having a good hit, 1e-10 is a value that gives reasonable confidence
my $E_cutoff = 1e-10;

# create ouput file
my $outFile = $blastout;
$outFile =~ s/RefGblast\.out/RefGblast.parsed/;
$outFile = $dir."/".$outFile;
open (OUT, ">$outFile") or die "can not open file $outFile!\n";

my @keep = (); # query should be kept for further analysis
my @known = (); # queries that are significantly similar to Reference sequences
my $total_records = 0;

#print "parsing blast output files...\n\n";

my $input_file = $dir."/".$blastout;
my $report = new Bio::SearchIO(-format => 'blast', -file => $input_file, -report_type => 'blastn');

# Go through BLAST reports one by one        
while(my $result = $report->next_result) {# next query output
#	print "\\", $result->query_name, "\\ input\n\n";
	if ($result->query_name eq "") { # deals with situation where blast 1st report is empty
		next;
	}

	$total_records++;
	my $haveHit = 0;
	my $keep = 1;
	while(my $hit = $result->next_hit) {
		$haveHit = 1;
		# check whether the query should be kept for further analysis
		if ($hit->significance <= $E_cutoff) {
			$keep = 0;	
#			print $result->query_name, " similar to known, output information!\n\n";
			print OUT $result->query_name, "\t", $result->query_length, "\t$RefGenomeTaxonomy\t$RefGenomeTaxonomy\t".$hit->name."\t".$hit->significance,"\n";
		}
		last; # only need to look at the first hit
	}

	if ($haveHit) {
		if ($keep) {
			push @keep, $result->query_name;
#			print $result->query_name, " keep!\n\n";
		}	
		else {		
			push @known, $result->query_name;
		}
	}
	else { # does not have a hit, keep for further analysis
		push @keep, $result->query_name;
#		print $result->query_name, " keep!\n\n";
	}	

}
print OUT "# Summary: ", scalar @keep, " out of $total_records ", scalar @keep/$total_records, " is saved for BLASTN analysis.\n";

close OUT;

# generate a fasta file that contains all the non-Reference sequences
# read in blastn input sequences
my $file = $blastout;
$file =~ s/\.RefGblast\.out//;
$file = $dir."/".$file.".fa";
my %seq = &read_FASTA_data($file);

$outFile = $blastout;
$outFile =~ s/\.RefGblast\.out//;
$outFile = $dir."/".$outFile.".RefGfiltered.fa";
open (OUT2, ">$outFile") or die "can not open file $outFile!\n";
foreach my $seq_name (@keep) {
	if ($seq_name eq "") { # deals with situation where blast 1st report is empty
		next;
	}
	print OUT2 ">$seq_name\n";
	print OUT2 $seq{$seq_name}, "\n";
}
close OUT2;


exit;


############################################################################
# subroutines
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
	    	@rows = split (/\s/, $line);	
	    	my $contigName = shift @rows;
	    	my $contigSeq = join("", @rows);
	    	$contigSeq =~ s/\s//g; #remove white space
	    	$fastaSeq{$contigName} = $contigSeq;
		}
    }

    # check for correctness
#	 foreach my $key (keys %fastaSeq){
#	      print "$key \t $fastaSeq{$key}\n";
#	 }

    #reset the read seperator
    $/ = $oldseperator;
	   
 	close FastaFile;
    return %fastaSeq;
}


