#!/usr/bin/perl -w

use strict;
use Bio::SearchIO;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;

my $Usage = '
Usage:
perl script <HOME> <dir> <blast output file>

<HOME> = directory used for temporary taxonomy files
<dir>  = directory that blast output file resides in (no trailing "/")
<blast output file> = blastn output file name
';

die $Usage unless scalar @ARGV == 3;
my ($HOME, $dir, $blastout) = @ARGV;

###################################################################################
# Local configuration (NO SQL needed)
my $database_dir = "/storage1/fs1/songcao/Active/Database/taxonomy";

# BLAST database prefix (NOT .nal; do NOT include extension)
# Example: if you have nt.nal and nt.00.nsq in /storage1/fs1/songcao/Active/Database,
# set this to /storage1/fs1/songcao/Active/Database/nt
my $blastdb_nt = "/storage1/fs1/songcao/Active/Database/nt/nt";

###################################################################################
# my $HOME = $ENV{HOME};

my %assignment = ();

# cutoff value for having a good hit
my $E_cutoff = 1e-10;
my $havedefinedHit = 0; # Song added

# create output file
my $outFile = $blastout;
$outFile =~ s/blastn\.out/blastn.parsed/;
$outFile = $dir . "/" . $outFile;
open(OUT, ">$outFile") or die "can not open file $outFile!\n";

# create a tmp directory in the home directory if tmp does not exist
if (! -d $HOME . "/taxo") {
    `mkdir $HOME"/taxo"`;
}

# get a Taxon from a Bio::DB::Taxonomy object (flatfile nodes/names)
my $dbh = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -directory => "$HOME/taxo",
    -nodesfile => "$database_dir/nodes.dmp",
    -namesfile => "$database_dir/names.dmp",
);

# Cache taxid lookups to avoid repeated blastdbcmd calls
my %ENTRY2TAX = ();

my @keep_for_tblastx = (); # query should be kept for further analysis
my @known = ();            # queries that are significantly similar to known sequences
my $total_records = 0;

print "parsing blast output files...\n\n";

my $input_file = $dir . "/" . $blastout;
my $report = Bio::SearchIO->new(
    -format      => 'blast',
    -file        => $input_file,
    -report_type => 'blastn'
);

# Go through BLAST reports one by one
while (my $result = $report->next_result) { # next query output
    $total_records++;
    my $haveHit = 0;
    my $keep_for_tblastx = 1;
    %assignment = ();

    # only take the best hits
    my $best_e = 100;
    my $hit_count = 0;
    $havedefinedHit = 0; # song added

    while (my $hit = $result->next_hit) {

        my $hit_name = $hit->name; # often: gi|num|db|acc| or ref|acc| ...
        my @temp_arr = split(/\|/, $hit_name);

        # skip data from pdb database
        if (defined $temp_arr[2] && $temp_arr[2] eq "pdb") {
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

            # only get best hits (same logic as original, with Song tweak)
            if ($hit->significance == $best_e || ($hit->significance <= $E_cutoff && $havedefinedHit == 1)) {

                # Determine an entry identifier to query blastdbcmd:
                # Prefer GI if present, otherwise use accession
                my $entry = parse_blast_entry_id($hit_name);

                # Lookup taxid from local BLAST database (no SQL)
                my $taxID = get_taxid_from_blastdb($blastdb_nt, $entry);

                if ($taxID) {
                    my $taxon_obj = $dbh->get_taxon(-taxonid => $taxID);

                    if (!(defined $taxon_obj)) {
                        my $description = "undefined taxon " . $hit->description . "\t" . $hit->name . "\t" . $hit->significance;
                        $assignment{"other"} = $description;
                    } else {
                        my $tree_function = Bio::Tree::Tree->new();
                        my @lineage = $tree_function->get_lineage_nodes($taxon_obj);

                        if (scalar @lineage) {
                            &PhyloType(\@lineage, $hit, $best_e, $dbh, \%assignment);
                        }
                    }
                } else {
                    # situations where entry does not have a corresponding taxid
                    my $desc = $hit->description . "\t" . $hit->name . "\t" . $hit->significance;
                    $assignment{"other"} = $desc;
                }
            } else {
                last;
            }
        } # finish phylotype for given hit
    }  # finish all hits

    # consolidate assignment
    my $num_assignment = keys %assignment;
    if ($num_assignment > 1) { # have multiple assignment
        my $has_specific = 0;
        my $has_other = 0;
        if ((defined $assignment{"Bacteria"}) || (defined $assignment{"Artificial"}) || (defined $assignment{"Fungi"}) || (defined $assignment{"Homo"}) || (defined $assignment{"Mus"}) || (defined $assignment{"Phage"}) || (defined $assignment{"Viruses"})) {
            $has_specific = 1;
        }
        if (defined $assignment{"other"}) {
            $has_other = 1;
        }

        # If a sequence hits virus and any other species with the same e value, assign "Ambiguous"
        if (((defined $assignment{"Bacteria"}) || (defined $assignment{"Fungi"}) || (defined $assignment{"Mus"}) || (defined $assignment{"Phage"}) || (defined $assignment{"other"})) && (defined $assignment{"Viruses"})) {
            $assignment{"Ambiguous"} = $assignment{"Viruses"};
            delete $assignment{"Viruses"};
        }
        if (((defined $assignment{"Viruses"}) || (defined $assignment{"Fungi"}) || (defined $assignment{"Mus"}) || (defined $assignment{"Phage"}) || (defined $assignment{"other"})) && (defined $assignment{"Bacteria"})) {
            $assignment{"Ambiguous"} = $assignment{"Bacteria"};
            delete $assignment{"Bacteria"};
        }

        # If assigned both a specific category and "other", keep only specific
        if ($has_specific && $has_other) {
            delete $assignment{"other"};
        }
    }

    # print out assignment for this query
    foreach my $assign (keys %assignment) {
        print OUT $result->query_name, "\t", $result->query_length, "\t", $assign, "\t", $assignment{$assign}, "\n";
    }

    if ($keep_for_tblastx) {
        push @keep_for_tblastx, $result->query_name;
    } else {
        push @known, $result->query_name;
    }
}

print OUT "# Summary: ", scalar @keep_for_tblastx, " out of $total_records ", scalar @keep_for_tblastx/$total_records, " is saved for next step analysis.\n";
close OUT;

# generate a fasta file that contains all the sequences that will be kept for further analysis
my $file = $blastout;
$file =~ s/\.blastn\.out//;
$file = $dir . "/" . $file . ".fa";
my %seq = &read_FASTA_data($file);

$outFile = $blastout;
$outFile =~ s/\.blastn\.out//;
$outFile = $dir . "/" . $outFile . ".BNfiltered.fa";
open(OUT2, ">$outFile") or die "can not open file $outFile!\n";

foreach my $seq_name (@keep_for_tblastx) {
    print OUT2 ">$seq_name\n";
    print OUT2 $seq{$seq_name}, "\n";
}
close OUT2;

exit;

############################################################################
sub read_FASTA_data {
    my $fastaFile = shift @_;

    my $oldseperator = $/;
    $/ = ">";

    my %fastaSeq;
    open(FastaFile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";

    while (my $line = <FastaFile>) {
        if ($line =~ /^\s*$/) {
            next;
        } elsif ($line =~ /^\s*#/) {
            next;
        } elsif ($line ne ">") {
            chomp $line;
            my @rows = split(/\s/, $line);
            my $contigName = shift @rows;
            my $contigSeq  = join("", @rows);
            $contigSeq =~ s/\s//g;
            $fastaSeq{$contigName} = $contigSeq;
        }
    }

    $/ = $oldseperator;
    close FastaFile;
    return %fastaSeq;
}

############################################################################
# Parse an identifier to query blastdbcmd.
# Prefer GI if present; otherwise use accession (ref|ACC| etc.).
sub parse_blast_entry_id {
    my ($hit_name) = @_;
    return undef if (!defined $hit_name || $hit_name eq "");

    # Common legacy format: gi|12345|db|ACC|
    if ($hit_name =~ /^gi\|(\d+)\|/) {
        return $1;
    }

    # Common formats: ref|ACC|  gb|ACC|  emb|ACC|  dbj|ACC|
    # Take the second field as accession
    my @t = split(/\|/, $hit_name);
    if (scalar @t >= 2 && defined $t[1] && $t[1] ne "") {
        return $t[1];
    }

    # Fallback: use whole name (sometimes BLAST stores accession directly)
    return $hit_name;
}

############################################################################
# Get taxid from BLAST database without SQL using blastdbcmd.
# Requires BLAST+ 'blastdbcmd' to be in PATH.
sub get_taxid_from_blastdb {
    my ($db_prefix, $entry) = @_;
    return undef if (!defined $db_prefix || $db_prefix eq "");
    return undef if (!defined $entry || $entry eq "");

    if (exists $ENTRY2TAX{$entry}) {
        return $ENTRY2TAX{$entry};
    }

    # %T returns taxid; suppress stderr
    my $cmd = "blastdbcmd -db '$db_prefix' -entry '$entry' -outfmt '%T' 2>/dev/null";
    my $taxid = `$cmd`;
    chomp($taxid);

    if (!defined $taxid || $taxid eq "" || $taxid !~ /^\d+$/) {
        $taxid = undef;
    }

    $ENTRY2TAX{$entry} = $taxid;
    return $taxid;
}

###############################################################################
# subroutine to determine the taxonomy lineage for a given blast hit
sub PhyloType {
    my ($lineage_ref, $hit_ref, $best_e, $dbh_taxonomy, $assignment_ref) = @_;
    my $description = "";
    my $node_id;
    my $obj;
    my $name;
    my $assigned = 0;

    my $Lineage = "";
    for (my $i = 0; $i <= $#$lineage_ref; $i++) {
        my $temp_node_id = $lineage_ref->[$i]->id;
        my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid => $temp_node_id);
        my $temp_name = $temp_obj->scientific_name;
        $Lineage .= $temp_name . ";";
    }

    if ($Lineage =~ /Mimiviridae/i || $Lineage =~ /Phycodnaviridae/i || $Lineage =~ /marseillevirus/i || $Lineage =~ /Iridoviridae/i) {
        $havedefinedHit = 1; # song added
    }

    # check to see if it is a human sequence
    if (scalar @{$lineage_ref} >= 4) {
        $node_id = $lineage_ref->[3]->id;
        $obj = $dbh_taxonomy->get_taxon(-taxonid => $node_id);
        $name = $obj->scientific_name;

        if ($name eq "Metazoa") {
            # make assignment
            for (my $i = 0; $i <= $#$lineage_ref; $i++) {
                my $temp_node_id = $lineage_ref->[$i]->id;
                my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid => $temp_node_id);
                my $temp_name = $temp_obj->scientific_name;

                if ($temp_name eq "Homo") {
                    if (!defined $assignment_ref->{"Homo"}) {
                        $description .= "Homo\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                        $assignment_ref->{"Homo"} = $description;
                    }
                    $assigned = 1;
                    last;
                }
            }

            if (!$assigned) {
                for (my $i = 0; $i <= $#$lineage_ref; $i++) {
                    my $temp_node_id = $lineage_ref->[$i]->id;
                    my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid => $temp_node_id);
                    my $temp_name = $temp_obj->scientific_name;

                    if ($temp_name eq "Mus") {
                        if (!defined $assignment_ref->{"Mus"}) {
                            $description .= "Mus\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                            $assignment_ref->{"Mus"} = $description;
                        }
                        $assigned = 1;
                        last;
                    }
                }
            }

            if (!$assigned) {
                if (!defined $assignment_ref->{"other"}) {
                    $description .= $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                    $assignment_ref->{"other"} = $description;
                }
                $assigned = 1;
            }
        }
    }

    # check to see if it is bacteria sequence
    if ((scalar @{$lineage_ref} >= 2) && (!$assigned)) {
        $node_id = $lineage_ref->[1]->id;
        $obj = $dbh_taxonomy->get_taxon(-taxonid => $node_id);
        $name = $obj->scientific_name;

        if ($name =~ /artificial sequences/i) {
            if (!defined $assignment_ref->{"Artificial"}) {
                $description = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                $assignment_ref->{"Artificial"} = $description;
            }
            $assigned = 1;
        }

        if ($name eq "Bacteria") {
            if (!defined $assignment_ref->{"Bacteria"}) {
                $description = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                $assignment_ref->{"Bacteria"} = $description;
            }
            $assigned = 1;
        }
    }

    # check to see if it is a phage virus sequence
    if (!$assigned) {
        $node_id = $lineage_ref->[0]->id;
        $obj = $dbh_taxonomy->get_taxon(-taxonid => $node_id);
        $name = $obj->scientific_name;

        if ($name eq "Viruses") {
            for (my $i = 0; $i <= $#$lineage_ref; $i++) {
                my $temp_node_id = $lineage_ref->[$i]->id;
                my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid => $temp_node_id);
                my $temp_name = $temp_obj->scientific_name;

                $description .= $temp_name . ";";
                if (($temp_name eq "Lipothrixviridae") || ($temp_name eq "Caudovirales") || ($temp_name eq "Corticoviridae") || ($temp_name eq "Cystoviridae") || ($temp_name eq "Inoviridae") || ($temp_name eq "Leviviridae") || ($temp_name eq "Microviridae") || ($temp_name eq "Tectiviridae") || ($temp_name =~ /phage/i)) {
                    if (!defined $assignment_ref->{"Phage"}) {
                        $description = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
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
        $obj = $dbh_taxonomy->get_taxon(-taxonid => $node_id);
        $name = $obj->scientific_name;

        if ($name eq "Viruses") {
            if (!defined $assignment_ref->{"Viruses"}) {
                $description = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                $assignment_ref->{"Viruses"} = $description;
            }
            $assigned = 1;
        }
    }

    # check to see if it is a fungi sequence
    if ((scalar @{$lineage_ref} >= 4) && (!$assigned)) {
        $node_id = $lineage_ref->[3]->id;
        $obj = $dbh_taxonomy->get_taxon(-taxonid => $node_id);
        $name = $obj->scientific_name;

        if ($name eq "Fungi") {
            if (!defined $assignment_ref->{"Fungi"}) {
                $description = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                $assignment_ref->{"Fungi"} = $description;
            }
            $assigned = 1;
        }
    }

    # if still not assigned, assigned to "other" category
    if (!$assigned) {
        if (!defined $assignment_ref->{"other"}) {
            $description = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
            $assignment_ref->{"other"} = $description;
        }
        $assigned = 1;
    }

    return $assigned;
}