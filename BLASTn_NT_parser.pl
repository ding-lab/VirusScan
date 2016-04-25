
#!/usr/bin/perl -w

use strict;
use Bio::SearchIO;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use DBI();

my $Usage = '
This script accepts a blastn output file and parse the information

perl script <dir><blast output file> 
<dir> = directory that blast output file resides in, without last "/"
<blast output file> = name of the blastn output file

';

die $Usage unless scalar @ARGV == 2;
my ($dir, $blastout) = @ARGV;

###################################################################################
# This section needs to be modified to use local configuration
my $database_dir = "/gscmnt/gc3027/info/medseq/taxdump_2014_01_08";

# open a connection to mysql database
my $dbh_mysql = DBI->connect("DBI:mysql:database=scao_taxondb;host=mysql1","scao", "asdf1234",{'RaiseError'=>1}) or die "Unable to connect $DBI::errstr\n";

###################################################################################
# Everhting below should not need modification
my $HOME =$ENV{HOME};

my %assignment = (); 

# cutoff value for having a good hit
my $E_cutoff = 1e-10;
my $havedefinedHit=0; # Song added
# create ouput file
my $outFile = $blastout;
$outFile =~ s/blastn\.out/blastn.parsed/;
$outFile = $dir."/".$outFile;
open (OUT, ">$outFile") or die "can not open file $outFile!\n";

# create a tmp directory in the home directory if tmp does not exist
if (! -d $HOME."/taxo") {
        `mkdir $HOME"/taxo"`;
}

# get a Taxon from a Bio::DB::Taxonomy object
my $dbh = Bio::DB::Taxonomy->new(-source => 'flatfile',
		-directory=> "$HOME/taxo",
		-nodesfile=> "$database_dir/nodes.dmp",
		-namesfile=> "$database_dir/names.dmp",
		);

my @keep_for_tblastx = (); # query should be kept for further analysis
my @known = (); # queries that are significantly similar to known sequences
my $total_records = 0;

print "parsing blast output files...\n\n";

my $input_file = $dir."/".$blastout;
my $report = new Bio::SearchIO(-format => 'blast', -file => $input_file, -report_type => 'blastn');

# Go through BLAST reports one by one        
while(my $result = $report->next_result) {# next query output
	$total_records++;
	my $haveHit = 0;
	my $keep_for_tblastx = 1;
	%assignment = ();

	# only take the best hits
	my $best_e = 100;
	my $hit_count = 0;
        $havedefinedHit=0; #song added; 	
	while(my $hit = $result->next_hit) {
		# from hit name get hit gi number
		my $hit_name = $hit->name; # gi|num|database|accessionNum|
		my @temp_arr = split(/\|/, $hit_name);
		my $gi = $temp_arr[1];
		#print $gi,"\n";
		if ($temp_arr[2] eq "pdb") { # skip data from pdb database
			next;
		}
		$haveHit = 1;
		$hit_count++;
		if ($hit_count == 1) {
			$best_e = $hit->significance;
		}

		# check whether the hit should be kept
		if ($best_e <= $E_cutoff) { # similar to known, need Phylotyped
			$keep_for_tblastx = 0;

#			print $result->query_name, " similar to known, output information!\n\n";
#			print "the $hit_count hit, $best_e \n"; 
			if ($hit->significance == $best_e || ($hit->significance <= $E_cutoff && $havedefinedHit==1)) { # only get best hits #song changed
				# from gi get taxonomy lineage
				my $sth = $dbh_mysql->prepare("SELECT * FROM gi_taxid_nucl where gi = $gi");
				$sth->execute();
				my $ref = $sth->fetchrow_hashref();
#				print "gi = $ref->{'gi'}, taxid = $ref->{'tax_id'}\n";
				
				$sth->finish();
				my $taxID = $ref->{'tax_id'};
				if ($taxID) { # some gi don't have record in gi_taxid_nucl
#					print "taxID is $taxID\n";
					my $taxon_obj = $dbh->get_taxon(-taxonid => $taxID);

					if (!(defined $taxon_obj)) {
#						die "unable to get taxon_obj object\n";
						my $description = "undefined taxon ".$hit->description."\t".$hit->name."\t".$hit->significance;
						$assignment{"other"} = $description; 
					}

					else {
						my $tree_function = Bio::Tree::Tree->new();
						my @lineage = $tree_function->get_lineage_nodes($taxon_obj);
						# each lineage node is a Bio::Tree::NodeI object

						#if($gi eq "61741475") {
						#print "hit gi is $gi\n";
						#print "id is ", $taxon_obj->id, "\n";
						#print "rank is ", $taxon_obj->rank, "\n";
						#print "divison is ", $taxon_obj->division, "\n\n";
						#print "lineage is @lineage\n"; 
						#<STDIN>;
						#}

						if (scalar @lineage) {				
#							print "PhyloTyped, don't save for further analysis\n";
							&PhyloType(\@lineage,$hit, $best_e, $dbh_mysql, $dbh, \%assignment);
						}
						#}
					}
				}	
				else { # for situations that gi does not have corresponding taxid
#					print $result->query_name, " ", $hit->name, "\n";
#					print "gi = $ref->{'gi'}, taxid = $ref->{'tax_id'}\n";
#					print "hit gi is $gi\n";
					my $desc = $hit->description."\t".$hit->name."\t".$hit->significance;
#					print $result->query_name, "\t", $desc, "\n";
					$assignment{"other"} = $desc;
				} 
			}
			else {
				last;
			}
		} # finish phylotype for given hit
	}  # finish all hits

#	foreach my $key (keys %assignment) {
#		print "after parsing ", $key, "\t", $assignment{$key},"\n";
#	}
	# consolidate assignment
	# If a query is assigned both Homo and Primates, it will be reported as Homo only
	# If a query is assigned a real taxon name and "other" for reason like"other sequences;
	# artificial sequences", or no taxon id in taxon database it will be reported only as 
	# the real taxon name
	my $num_assignment = keys %assignment;
	if ($num_assignment > 1) { # have multiple assignment
		# handle the situation that assigned both a specific category and "other"
		# only specific category will be save.
		my $has_specific = 0;
		my $has_other = 0;
		if ((defined $assignment{"Bacteria"}) || (defined $assignment{"Artificial"}) || (defined $assignment{"Fungi"}) || (defined $assignment{"Homo"}) || (defined $assignment{"Mus"}) || (defined $assignment{"Phage"}) || (defined $assignment{"Viruses"})) {
			$has_specific = 1;
		}
		if (defined $assignment{"other"}) {
			$has_other = 1;
		}
		#################################################################
		# If a sequence hits virus and any other species with the same e value, 
		# the sequence is assigned to "Ambiguous" category. cai added 12/2010
		#remove human since we have done extensive filtering for human sequence 10/12/2014

		if (((defined $assignment{"Bacteria"}) || (defined $assignment{"Fungi"}) || (defined $assignment{"Mus"}) || (defined $assignment{"Phage"}) || (defined $assignment{"other"})) && (defined $assignment{"Viruses"})) {
			$assignment{"Ambiguous"} = $assignment{"Viruses"};
			delete $assignment{"Viruses"};
		}
		if (((defined $assignment{"Viruses"}) || (defined $assignment{"Fungi"}) || (defined $assignment{"Mus"}) || (defined $assignment{"Phage"}) || (defined $assignment{"other"})) && (defined $assignment{"Bacteria"})) {
                        $assignment{"Ambiguous"} = $assignment{"Bacteria"};
                        delete $assignment{"Bacteria"};
                }
		#################################
		if ($has_specific && $has_other) {
			delete $assignment{"other"}; 
		}

	}

#	foreach my $key (keys %assignment) {
#		print "after consolidateion ", $key, "\t", $assignment{$key},"\n";
#	}

	# print out assignment for this query
	foreach my $assign (keys %assignment) {
		print OUT $result->query_name, "\t", $result->query_length, "\t", $assign, "\t", $assignment{$assign}, "\n";
#		print $result->query_name, "\t", $result->query_length, "\t", $assign, "\t", $assignment{$assign}, "\n";

	}
	
	if ($keep_for_tblastx) {
		push @keep_for_tblastx, $result->query_name;
#		print $result->query_name, " keep_for_tblastx!\n\n";
	}	
	else {		
		push @known, $result->query_name;
	}
#}
}
print OUT "# Summary: ", scalar @keep_for_tblastx, " out of $total_records ", scalar @keep_for_tblastx/$total_records, " is saved for next step analysis.\n";

close OUT;

# generate a fasta file that contains all the sequences that will be kept for further analysis
# read in blast input sequences
my $file = $blastout;
$file =~ s/\.blastn\.out//;
$file = $dir."/".$file.".fa";
my %seq = &read_FASTA_data($file);

$outFile = $blastout;
$outFile =~ s/\.blastn\.out//;
$outFile = $dir."/".$outFile.".BNfiltered.fa";
open (OUT2, ">$outFile") or die "can not open file $outFile!\n";
foreach my $seq_name (@keep_for_tblastx) {
	print OUT2 ">$seq_name\n";
	print OUT2 $seq{$seq_name}, "\n";
}
close OUT2;

$dbh_mysql->disconnect();

exit;


############################################################################
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

    # to check the correctness 
#	 foreach my $key (keys %fastaSeq){
#	      print "Here is the key for fasta seq: $key \t $fastaSeq{$key}\n";
#	 }

    #reset the read seperator
    $/ = $oldseperator;
    close FastaFile;

    return %fastaSeq;
}

	
###############################################################################
# subroutine to determine the taxonomy lineage for a given blast hit
sub PhyloType {
	my ($lineage_ref, $hit_ref, $best_e, $dbh_mysql, $dbh_taxonomy, $assignment_ref) = @_;
	my $description = "";
	my $node_id; 
	my $obj;
	my $name;
	my $assigned = 0;

	my $Lineage = "";
	for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
		my $temp_node_id = $lineage_ref->[$i]->id;
		my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
		my $temp_name = $temp_obj->scientific_name;
		$Lineage .= $temp_name.";";
            }					
	#print "linease is $Lineage\n";
   
        if($Lineage =~/Mimiviridae/i || $Lineage =~/Phycodnaviridae/i || $Lineage =~/marseillevirus/i || $Lineage =~/Iridoviridae/i) { $havedefinedHit=1; } #song added;  
		
	# check to see if it is a human sequence
	if (scalar @{$lineage_ref} >= 4) {
		$node_id = $lineage_ref->[3]->id;
		$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
		$name = $obj->scientific_name;
		if ($name eq "Metazoa") {
			# make assignment
			for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
				my $temp_node_id = $lineage_ref->[$i]->id;
				my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
				my $temp_name = $temp_obj->scientific_name;
				#print "name = $temp_name\n";
				#<STDIN>;
				if ($temp_name eq "Homo") {
                    if(!defined $assignment_ref->{"Homo"}) {  # only keep the first best hit description, song added 1/7/2012
#						print "assigned to Homo\n\n";
						$description .= "Homo\t".$hit_ref->name."\t".$hit_ref->significance;
						$assignment_ref->{"Homo"} = $description; 
					}
					$assigned = 1;
					last; 
				}
			}
			if (!$assigned) {
				for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
					my $temp_node_id = $lineage_ref->[$i]->id;
					my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
					my $temp_name = $temp_obj->scientific_name;
#					print "name = $temp_name\n";
	
					if ($temp_name eq "Mus") {
                       if(!defined $assignment_ref->{"Mus"}) {  # only keep the first best hit description, song added 1/7/2012
#							print "assigned to Mus\n\n";
							$description .= "Mus\t".$hit_ref->name."\t".$hit_ref->significance;
							$assignment_ref->{"Mus"} = $description; 
						}
						$assigned = 1;
						last; 
					}
				}
			}
			if (!$assigned) {
				if(!defined $assignment_ref->{"other"})  { # only take the first best hit description
					$description .= $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
#					print "assigned to other\n\n";
					$assignment_ref->{"other"} = $description; 
				}
				$assigned = 1; 
			}
		}
	}

	# check to see if it is bacteria sequence
	if ((scalar @{$lineage_ref} >= 2)&&(!$assigned)) {
		$node_id = $lineage_ref->[1]->id;
		#print $node_id,"\n";
		$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
		$name = $obj->scientific_name;
		#print $name,"\n";

		if($name=~/artificial sequences/i) {
			if(!defined $assignment_ref->{"Artificial"})
			 {
			  $description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
              $assignment_ref->{"Artificial"} = $description;}
			 $assigned=1; 
		}

		if ($name eq "Bacteria") {
			if(!defined $assignment_ref->{"Bacteria"})  { # take the first best hit description   
				$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
				$assignment_ref->{"Bacteria"} = $description;
			}
			$assigned = 1; 
		}
	}

	# check to see if it is a phage virus sequence
	if (!$assigned) {
		$node_id = $lineage_ref->[0]->id;
		$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
		$name = $obj->scientific_name;
		if ($name eq "Viruses") {
			for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
				my $temp_node_id = $lineage_ref->[$i]->id;
				my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
				my $temp_name = $temp_obj->scientific_name;
				$description .= $temp_name.";";
				if (($temp_name eq "Lipothrixviridae")||($temp_name eq "Caudovirales")||($temp_name eq "Corticoviridae")||($temp_name eq "Cystoviridae")||($temp_name eq "Inoviridae")||($temp_name eq "Leviviridae")||($temp_name eq "Microviridae")||($temp_name eq "Tectiviridae")||($temp_name =~ /phage/i)) {
#					print "assigned to phage\n\n";
					if(!defined $assignment_ref->{"Phage"}) { # take the first best hit description
						$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
						$assignment_ref->{"Phage"} = $description;
					}
					$assigned = 1;
					last;
				}
			}
		}
	}

	# check to see if it is a virus sequence
	$description = "";
	if (!$assigned) {
		$node_id = $lineage_ref->[0]->id;
		$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
		$name = $obj->scientific_name;
		if ($name eq "Viruses") {
			if(!defined $assignment_ref->{"Viruses"}) { # take the first best hit description
				$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
				$assignment_ref->{"Viruses"} = $description;
			}
			$assigned = 1;
		}
	}

	# check to see if it is a fungi sequence
	if ((scalar @{$lineage_ref} >= 4)&&(!$assigned)) {
		$node_id = $lineage_ref->[3]->id;
		$obj = $dbh->get_taxon(-taxonid=>$node_id);
		$name = $obj->scientific_name;
		if ($name eq "Fungi") {
			if(!defined $assignment_ref->{"Fungi"}) { # take the first best hit description
				$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
				$assignment_ref->{"Fungi"} = $description;
			}
			$assigned = 1;
		}
	}

	# if still not assigned, assigned to "other" category
	if (!$assigned) {
		if(!defined $assignment_ref->{"other"}) {
			$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
			$assignment_ref->{"other"} = $description;
		}
		$assigned = 1;
	}
	
	return $assigned;
}


