
#!/usr/bin/perl
use strict;
use Switch;
use Bio::SearchIO;

my $usage = '
perl $full_path $blast_file
';

die $usage unless scalar @ARGV == 2;
my ( $full_path, $blast_file) = @ARGV;

# get all the viral read sequences
my %viral_reads_blastn = ();
my %viral_reads_blastx = ();
my %best_e_blastn = (); # viral_read_ID => best_e value for this read in blastn
my %best_e_blastx = (); # viral_read_ID => best_e value for this read in blastx
my @blast_files_blastn = (); # all blastn.out files
my @blast_files_blastx = (); # all blastx.out files
my @unassigned_reads = ();
####################################
my @ambiguous_reads = (); #cai added 12/2010
####################################

# read in original sequences
#my @temp = split("\/", $dir);
#my $lib_name = pop @temp;
# print "lib is $lib_name\n";
#my $fasta_file = $dir."/".$lib_name.".fa.cdhit_out.masked.goodSeq";
#my $fasta_file = $input_good_seq_fasta_file; #cai changed, added segmasker

#my %seq = &read_FASTA_data($fasta_file);

#my $out1 = $dir."/".$lib_name.".gi.AssignmentReport";
#open (OUT1, ">$out1") or die "can not open file $out1!\n";
#my $OUT2 = $dir."/".$lib_name.".gi.ViralReads_all.fa";
#open (OUT2, ">$OUT2") or die "can not open file $OUT2!\n";
#my $OUT3 = $dir."/".$lib_name.".gi.unassigned.fa";
#open (OUT3, ">$OUT3") or die "can not open file $OUT3!\n";

##################################cai added 12/2010
#my $out4 = $dir."/".$lib_name.".gi.AmbiguousReads_all.fa";
#open (OUT4, ">$out4") or die "can not open file $out4!\n";
##################################


# category => num of sequence assigned to this category by blastn
my %blastn = ( 
	"Bacteria" => 0,
	"Fungi" => 0,
	"Homo" => 0,
	"Mus" => 0,
	"Phage" => 0,
	"Viruses" => 0,
	"other" => 0,
	"unassigned" => 0,
	##################################cai added 12/2010
	"Ambiguous" => 0,
	##################################cai added
);

# category => num of sequence assigned to this category by blastn of Reference genome
my %blastn_RefG = ();
foreach my $key (keys %blastn) {
        $blastn_RefG{$key} = 0;
}

# category => num of sequence assigned to this category by tblastx of viral genome
my %blastx = ();
foreach my $key (keys %blastn) {
	$blastx{$key} = 0;
}

# viral_lineage => number of reads assigned to this lineage in the library
my %num_reads = ();
my %blast_readinfo =();    # readID => information about this read
my %lineage_blastn = ();    # lineage => [read ID]
my %lineage_gi = (); 
my %lineage_blastx = ();    # lineage => [read ID]
#if ($blast_file =~ /blastn\.parsed$/) {
my $blast_out = $blast_file;
$blast_out =~ s/\.blastn\.parsed/\.blastn\.out/;
$blast_out = $full_path."/".$blast_out;
my $blast_s = $blast_file;
print $blast_s,"\n";
$blast_s =~ s/\.blastn\.parsed/\.blastn\.summary/;
$blast_s = $full_path."/".$blast_s;
#print $blast_file,"\n";
#print $blast_out,"\n";
#print $blast_s,"\n";

#open(OUT,">$blast_s"); 

#foreach my $id (keys %blast_readinfo) 
 #{
  # print OUT $id,"\t",$blast_readinfo{$id},"\n"; 
 #}

push @blast_files_blastn, $blast_out;
my $parsed = $full_path."/".$blast_file;
&collect_information($parsed, \%blastn, \%viral_reads_blastn, \%best_e_blastn, \%lineage_blastn, \%lineage_gi, \%num_reads, \@unassigned_reads, \@ambiguous_reads);

&get_viral_read_info( \@blast_files_blastn, "blastn", \%viral_reads_blastn, \%best_e_blastn, \%blast_readinfo);

open(OUT,">$blast_s");

foreach my $id (keys %blast_readinfo)
 {
   print OUT $blast_readinfo{$id};
 }

print OUT "Finished summary\n"; 
close OUT; 
#####################################################################################
# collecte information from given directory
sub collect_information {
	##################################cai changed 12/2010
	my ($infile, $category_hash_ref, $viral_reads_hash_ref, $best_e_hash_ref, $lineage_hash_ref, $lineage_hash_gi, $num_reads_hash_ref, $unassigned_reads_arr_ref, $ambiguous_reads_arr_ref) = @_;
	##################################
	open (IN, $infile) or die "can not open file $infile!\n";
	while (<IN>) {
		if ($_ =~ /#/) { # skip comment line
			next;
		}
		chomp;
		my ($read_ID, $length, $category, $lineage, $hit_name, $e_value) = split("\t", $_);
#		print "readID = $read_ID, length = $length, category = $category, lineage = $lineage, hit name = $hit_name, e = $e_value\n";
		my $gid=0; 
                if($hit_name=~/gi\|(\d+)\|/) { $gid=$1; $lineage_hash_gi->{$gid}=$lineage; } 
                
		switch ($category ) {
			case "Bacteria" { $category_hash_ref->{"Bacteria"}++	}
			case "Fungi" { $category_hash_ref->{"Fungi"}++ }
			case "Homo" { $category_hash_ref->{"Homo"}++ }
			case "Mus" { $category_hash_ref->{"Mus"}++ }
			case "Phage" {$category_hash_ref->{"Phage"}++ }
			case "Viruses" { $category_hash_ref->{"Viruses"}++ }
			case "other" {$category_hash_ref->{"other"}++ }
			case "unassigned" {$category_hash_ref->{"unassigned"}++}
			case "Ambiguous" {$category_hash_ref->{"Ambiguous"}++ } #cai added
		}
	
		if (($category eq "Viruses" || $category eq "Bacteria") && $gid!=0) {

			$viral_reads_hash_ref->{$read_ID} = 1;

			$best_e_hash_ref->{$read_ID} = $e_value;

			if (!(defined $lineage_hash_ref->{$gid})) {
				$lineage_hash_ref->{$gid} = [$read_ID];
			}
			else {
				push @{$lineage_hash_ref->{$gid}}, $read_ID;
			}

			if (defined $num_reads_hash_ref->{$gid}) {
				$num_reads_hash_ref->{$gid}++;
			}
			else {
				$num_reads_hash_ref->{$gid} = 1;
			}

		##################################cai added 12/2010
		}elsif ($category eq "Ambiguous"){
			push @{$ambiguous_reads_arr_ref}, $read_ID;
		##################################
		}elsif ($category eq "unassigned") {
			push @{$unassigned_reads_arr_ref}, $read_ID;
		}
	}
	close IN;
}

#############################################################################
# get detailed information about each viral read
sub get_viral_read_info {
	my ($report_file_ref, $report_type, $viral_reads_hash_ref, $best_e_hash_ref, $blast_readinfo_hash_ref) = @_;
	my $report; # blast report object
	foreach my $file (@{$report_file_ref}) {
		$report = new Bio::SearchIO(-format => 'blast', -file => $file, -report_type => $report_type);
		# Go through BLAST reports one by one        
		while(my $result = $report->next_result) {# next query output
			my $read_ID = $result->query_name;
			if (defined $viral_reads_hash_ref->{$read_ID}) {
				my $desc = "";
				my $hit_count = 0;
				while (my $hit = $result->next_hit()) {
					if ($hit->significance() == $best_e_hash_ref->{$read_ID}) {
						$hit_count++;
						# for those with hundreads hits, only take the first 100
						if ($hit_count == 2) {
							last; 
						}
						$desc .= $result->query_name()."\t";
						$desc .= $result->query_length()."\t";
						$desc .= $hit->name()."\t";
						$desc .= $hit->length()."\t";
						$desc .= $hit->description(60)."\t";
						while (my $hsp = $hit->next_hsp()) {
							$desc .= $hsp->length('hit')."\t";
							my $percent_id = sprintf("%4.1f", $hsp->percent_identity());
							$desc .= $percent_id."\%\t[";
							$desc .= $hsp->start('hit')."\t";
							$desc .= $hsp->end('hit')."]\t";
							$desc .= $hsp->evalue()."\n";
							last;
						}
					}
				}
				$blast_readinfo_hash_ref->{$read_ID} = $desc;
			}
		}
	}
}
