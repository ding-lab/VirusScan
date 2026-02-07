#!/usr/bin/perl -w
use strict;

use Bio::SearchIO;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use File::Path qw(make_path);
use File::Spec;

my $Usage = '
Usage:
perl script <HOME> <dir> <blast output file>

<HOME> = directory used for temporary taxonomy files
<dir>  = directory that blast output file resides in (no trailing "/")
<blast output file> = blastn output file name
';

die $Usage unless @ARGV == 3;
my ($HOME, $dir, $blastout) = @ARGV;

###################################################################################
# Local configuration (NO SQL needed)
my $database_dir = "/storage1/fs1/songcao/Active/Database/taxonomy";

# BLAST database prefix (NOT .nal; do NOT include extension)
my $blastdb_nt   = "/storage1/fs1/songcao/Active/Database/nt/nt";
###################################################################################

# cutoff value for having a good hit
my $E_cutoff = 1e-10;
my $havedefinedHit = 0;

# create output file
my $outFile = $blastout;
$outFile =~ s/blastn\.out/blastn.parsed/;
$outFile = File::Spec->catfile($dir, $outFile);
open(my $OUT, ">", $outFile) or die "cannot open $outFile: $!\n";

# create a tmp taxonomy cache directory
my $taxo_dir = File::Spec->catdir($HOME, "taxo");
if (! -d $taxo_dir) {
    make_path($taxo_dir) or die "cannot create $taxo_dir: $!\n";
}

# get a Taxon from a Bio::DB::Taxonomy object (flatfile nodes/names)
my $dbh = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -directory => $taxo_dir,
    -nodesfile => File::Spec->catfile($database_dir, "nodes.dmp"),
    -namesfile => File::Spec->catfile($database_dir, "names.dmp"),
);

# Cache taxid lookups to avoid repeated blastdbcmd calls
my %ENTRY2TAX;

my @keep_for_tblastx; # query should be kept for further analysis
my @known;            # queries that are significantly similar to known sequences
my $total_records = 0;

print "parsing blast output files...\n\n";

my $input_file = File::Spec->catfile($dir, $blastout);
my $report = Bio::SearchIO->new(
    -format      => 'blast',
    -file        => $input_file,
    -report_type => 'blastn'
);

# Go through BLAST reports one by one
while (my $result = $report->next_result) {
    $total_records++;

    my %assignment;
    my $keep_for_tblastx = 1;

    # only take the best hits
    my $best_e = 100;
    my $hit_count = 0;
    $havedefinedHit = 0;

    while (my $hit = $result->next_hit) {

        my $hit_name = $hit->name; # often: gi|num|db|acc| or ref|acc| ...
        my @temp_arr = split(/\|/, $hit_name);

        # skip data from pdb database
        if (defined $temp_arr[2] && $temp_arr[2] eq "pdb") {
            next;
        }

        $hit_count++;
        if ($hit_count == 1) {
            $best_e = $hit->significance;
        }

        # check whether the hit should be kept
        if ($best_e <= $E_cutoff) {
            $keep_for_tblastx = 0;

            # keep only best hits (Song logic kept)
            if ($hit->significance == $best_e || ($hit->significance <= $E_cutoff && $havedefinedHit == 1)) {

                my $entry = parse_blast_entry_id($hit_name);

                # Lookup taxid from local BLAST database (no SQL)
                my $taxID = get_taxid_from_blastdb($blastdb_nt, $entry, \%ENTRY2TAX);

                if ($taxID) {
                    my $taxon_obj = $dbh->get_taxon(-taxonid => $taxID);

                    if (!defined $taxon_obj) {
                        my $desc = "undefined taxon " . $hit->description . "\t" . $hit->name . "\t" . $hit->significance;
                        $assignment{"other"} = $desc;
                    } else {
                        my @lineage_names = build_clean_lineage_names($taxon_obj, $dbh);
                        if (@lineage_names) {
                            PhyloType(\@lineage_names, $hit, $best_e, \%assignment);
                        }
                    }
                } else {
                    my $desc = $hit->description . "\t" . $hit->name . "\t" . $hit->significance;
                    $assignment{"other"} = $desc;
                }

            } else {
                last;
            }
        }
    }

    # consolidate assignment
    my $num_assignment = scalar keys %assignment;
    if ($num_assignment > 1) {
        my $has_specific = 0;
        my $has_other    = 0;

        if (defined $assignment{"Bacteria"} ||
            defined $assignment{"Artificial"} ||
            defined $assignment{"Fungi"} ||
            defined $assignment{"Homo"} ||
            defined $assignment{"Mus"} ||
            defined $assignment{"Phage"} ||
            defined $assignment{"Viruses"}) {
            $has_specific = 1;
        }
        $has_other = 1 if defined $assignment{"other"};

        # If a sequence hits virus and any other species with the same e value, assign "Ambiguous"
        if ((defined $assignment{"Viruses"}) &&
            (defined $assignment{"Bacteria"} || defined $assignment{"Fungi"} || defined $assignment{"Mus"} || defined $assignment{"Phage"} || defined $assignment{"other"})) {
            $assignment{"Ambiguous"} = $assignment{"Viruses"};
            delete $assignment{"Viruses"};
        }

        if ((defined $assignment{"Bacteria"}) &&
            (defined $assignment{"Viruses"} || defined $assignment{"Fungi"} || defined $assignment{"Mus"} || defined $assignment{"Phage"} || defined $assignment{"other"})) {
            $assignment{"Ambiguous"} = $assignment{"Bacteria"};
            delete $assignment{"Bacteria"};
        }

        # If assigned both a specific category and "other", keep only specific
        if ($has_specific && $has_other) {
            delete $assignment{"other"};
        }
    }

    # print out assignment for this query
    for my $assign (keys %assignment) {
        print $OUT $result->query_name, "\t", $result->query_length, "\t", $assign, "\t", $assignment{$assign}, "\n";
    }

    if ($keep_for_tblastx) {
        push @keep_for_tblastx, $result->query_name;
    } else {
        push @known, $result->query_name;
    }
}

if ($total_records > 0) {
    print $OUT "# Summary: ", scalar(@keep_for_tblastx), " out of $total_records ",
               (scalar(@keep_for_tblastx)/$total_records), " is saved for next step analysis.\n";
}
close $OUT;

# generate a fasta file that contains all the sequences that will be kept for further analysis
my $fa_file = $blastout;
$fa_file =~ s/\.blastn\.out$//;
$fa_file = File::Spec->catfile($dir, $fa_file . ".fa");

my %seq = read_FASTA_data($fa_file);

my $filtered = $blastout;
$filtered =~ s/\.blastn\.out$//;
$filtered = File::Spec->catfile($dir, $filtered . ".BNfiltered.fa");

open(my $OUT2, ">", $filtered) or die "cannot open $filtered: $!\n";
for my $seq_name (@keep_for_tblastx) {
    next unless exists $seq{$seq_name};
    print $OUT2 ">$seq_name\n", $seq{$seq_name}, "\n";
}
close $OUT2;

exit;

############################################################################
sub read_FASTA_data {
    my ($fastaFile) = @_;

    my $oldsep = $/;
    $/ = ">";

    my %fastaSeq;
    open(my $FH, "<", $fastaFile) or die "Can't open FASTA file: $fastaFile\n";

    while (my $chunk = <$FH>) {
        next if $chunk =~ /^\s*$/;
        next if $chunk =~ /^\s*#/;

        chomp $chunk;
        next if $chunk eq ">";

        my @rows = split(/\s+/, $chunk);
        my $name = shift @rows;
        my $seq  = join("", @rows);
        $seq =~ s/\s//g;

        $fastaSeq{$name} = $seq if defined $name && $name ne "";
    }

    close $FH;
    $/ = $oldsep;
    return %fastaSeq;
}

############################################################################
# Parse an identifier to query blastdbcmd.
# Prefer GI if present; otherwise use accession (ref|ACC| etc.).
sub parse_blast_entry_id {
    my ($hit_name) = @_;
    return undef if (!defined $hit_name || $hit_name eq "");

    # Legacy: gi|12345|...
    if ($hit_name =~ /^gi\|(\d+)\|/) {
        return $1;
    }

    # ref|ACC|  gb|ACC|  emb|ACC|  dbj|ACC|
    my @t = split(/\|/, $hit_name);
    if (@t >= 2 && defined $t[1] && $t[1] ne "") {
        return $t[1];
    }

    return $hit_name;
}

############################################################################
# Get taxid from BLAST database without SQL using blastdbcmd.
# Requires BLAST+ 'blastdbcmd' to be in PATH.
sub get_taxid_from_blastdb {
    my ($db_prefix, $entry, $cache_ref) = @_;
    return undef if (!defined $db_prefix || $db_prefix eq "");
    return undef if (!defined $entry || $entry eq "");

    if (exists $cache_ref->{$entry}) {
        return $cache_ref->{$entry};
    }

    # %T returns taxid; suppress stderr
    my $cmd = "blastdbcmd -db '$db_prefix' -entry '$entry' -outfmt '%T' 2>/dev/null";
    my $taxid = `$cmd`;
    chomp($taxid);

    if (!defined $taxid || $taxid eq "" || $taxid !~ /^\d+$/) {
        $taxid = undef;
    }

    $cache_ref->{$entry} = $taxid;
    return $taxid;
}

############################################################################
# Build a clean taxonomy lineage (root -> ... -> leaf) for one taxon_obj.
# Uses parent_id chain so we get a TRUE lineage, not a mixed set of nodes.
sub build_clean_lineage_names {
    my ($taxon_obj, $dbh_taxonomy) = @_;
    return () if (!defined $taxon_obj);

    my @names;
    my %seen;

    my $tax = $taxon_obj;
    while (defined $tax) {
        my $tid  = $tax->id;
        my $name = $tax->scientific_name;

        last if (!defined $tid || $tid eq "");
        last if ($seen{$tid}++);  # avoid cycles

        push @names, $name if defined $name && $name ne "";

        my $pid = eval { $tax->parent_id };
        last if (!defined $pid || $pid eq "" || $pid == $tid);

        my $parent = $dbh_taxonomy->get_taxon(-taxonid => $pid);
        last if (!defined $parent);

        $tax = $parent;
    }

    @names = reverse @names; # root -> leaf
    return @names;
}

############################################################################
sub lineage_contains {
    my ($lineage_ref, $target) = @_;
    for my $n (@{$lineage_ref}) {
        return 1 if defined $n && $n eq $target;
    }
    return 0;
}

sub lineage_contains_regex {
    my ($lineage_ref, $re) = @_;
    for my $n (@{$lineage_ref}) {
        next unless defined $n;
        return 1 if $n =~ $re;
    }
    return 0;
}

###############################################################################
# Determine taxonomy class for a blast hit based on CLEAN lineage names.
sub PhyloType {
    my ($lineage_ref, $hit_ref, $best_e, $assignment_ref) = @_;

    my $Lineage = join(";", @{$lineage_ref}) . ";";

    # Song tweak: set "definedHit" when large DNA virus groups appear
    if ($Lineage =~ /Mimiviridae/i || $Lineage =~ /Phycodnaviridae/i || $Lineage =~ /marseillevirus/i || $Lineage =~ /Iridoviridae/i) {
        $havedefinedHit = 1;
    }

    # Artificial
    if (lineage_contains_regex($lineage_ref, qr/^artificial sequences$/i) ||
        lineage_contains_regex($lineage_ref, qr/artificial sequences/i)) {
        if (!defined $assignment_ref->{"Artificial"}) {
            $assignment_ref->{"Artificial"} = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
        }
        return 1;
    }

    # Human / Mouse (most specific)
    if (lineage_contains($lineage_ref, "Homo")) {
        if (!defined $assignment_ref->{"Homo"}) {
            $assignment_ref->{"Homo"} = "Homo\t" . $hit_ref->name . "\t" . $hit_ref->significance;
        }
        return 1;
    }
    if (lineage_contains($lineage_ref, "Mus")) {
        if (!defined $assignment_ref->{"Mus"}) {
            $assignment_ref->{"Mus"} = "Mus\t" . $hit_ref->name . "\t" . $hit_ref->significance;
        }
        return 1;
    }

    # Fungi
    if (lineage_contains($lineage_ref, "Fungi")) {
        if (!defined $assignment_ref->{"Fungi"}) {
            $assignment_ref->{"Fungi"} = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
        }
        return 1;
    }

    # Bacteria
    if (lineage_contains($lineage_ref, "Bacteria")) {
        if (!defined $assignment_ref->{"Bacteria"}) {
            $assignment_ref->{"Bacteria"} = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
        }
        return 1;
    }

    # Viruses / Phage
    if (@{$lineage_ref} && $lineage_ref->[0] && $lineage_ref->[0] eq "Viruses") {

        # Phage clades (keep list explicit; DO NOT use /phage/i loose match)
        my @phage_markers = qw(
            Caudovirales
            Caudoviricetes
            Caudoviridae
            Siphoviridae
            Myoviridae
            Podoviridae
            Microviridae
            Inoviridae
            Leviviridae
            Cystoviridae
            Tectiviridae
            Corticoviridae
            Lipothrixviridae
        );

        for my $m (@phage_markers) {
            if (lineage_contains($lineage_ref, $m)) {
                if (!defined $assignment_ref->{"Phage"}) {
                    $assignment_ref->{"Phage"} = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
                }
                return 1;
            }
        }

        if (!defined $assignment_ref->{"Viruses"}) {
            $assignment_ref->{"Viruses"} = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
        }
        return 1;
    }

    # other
    if (!defined $assignment_ref->{"other"}) {
        $assignment_ref->{"other"} = $Lineage . "\t" . $hit_ref->name . "\t" . $hit_ref->significance;
    }
    return 1;
}