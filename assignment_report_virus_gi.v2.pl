#!/usr/bin/perl
use strict;
# use Switch;   # REMOVED
use Bio::SearchIO;

my $usage = '
This script will read corresponding files in the given director and
generate a report. It will report in each library, for each category,
how many total sequence were assigned to this category, how many were
assigned by BLASTN, how many were assigned by TBLASTX, the range of
percent identity. It will also generate four fasta format files which
contain viral reads from blastn, tblastx, all viral reads and reads
that can not be assigned to any category.

perl script <sample dir> <input repeatmasked good seq(full path)> <ref genome>
<sample dir> = full path to the directory holding files for the given
               library
               e.g. .../S21_Rota_other
<ref genome> = 1. Human
               2. Mouse
               3. Worm (C. elegans, C. briggsae)
               4. Mouse lemur (Microcebus_murinus)
               5. sand fly (Lutzomyia longipalpis)
';

die $usage unless scalar @ARGV == 3;
my ( $dir, $input_good_seq_fasta_file, $ref_genome_choice ) = @ARGV;

# get all the viral read sequences
my %viral_reads_blastn = ();
my %viral_reads_blastx = ();

my %best_e_blastn = (); # viral_read_ID => best_e value for this read in blastn
my %best_e_blastx = (); # viral_read_ID => best_e value for this read in blastx

my @blast_files_blastn = (); # all blastn.out files
my @blast_files_blastx = (); # all blastx.out files

my @unassigned_reads = ();
my @ambiguous_reads = (); #cai added 12/2010

# read in original sequences
my @temp = split("\/", $dir);
my $lib_name = pop @temp;

my $fasta_file = $input_good_seq_fasta_file; #cai changed, added segmasker
my %seq = &read_FASTA_data($fasta_file);

my $out1 = $dir."/".$lib_name.".gi.AssignmentReport";
open (OUT1, ">$out1") or die "can not open file $out1!\n";
my $OUT2 = $dir."/".$lib_name.".gi.ViralReads_all.fa";
open (OUT2, ">$OUT2") or die "can not open file $OUT2!\n";
my $OUT3 = $dir."/".$lib_name.".gi.unassigned.fa";
open (OUT3, ">$OUT3") or die "can not open file $OUT3!\n";

my $out4 = $dir."/".$lib_name.".gi.AmbiguousReads_all.fa";
open (OUT4, ">$out4") or die "can not open file $out4!\n";

# category => num of sequence assigned to this category by blastn
my %blastn = (
    "Bacteria"   => 0,
    "Fungi"      => 0,
    "Homo"       => 0,
    "Mus"        => 0,
    "Phage"      => 0,
    "Viruses"    => 0,
    "other"      => 0,
    "unassigned" => 0,
    "Ambiguous"  => 0,
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
my %blast_readinfo =();     # readID => information about this read
my %lineage_blastn = ();    # lineage => [read ID]
my %lineage_gi = ();
my %lineage_blastx = ();    # lineage => [read ID]

opendir(DH, $dir) or die "Can not open dir $dir!\n";
foreach my $name (readdir DH) {
    my $full_path = $dir."/".$name;

    # Reference genome blast result
    if ($name =~ /goodSeq_RefGblast$/) {
        opendir (RefGDIR, $full_path) or die "can not open dir $full_path!\n";
        foreach my $blast_file (readdir RefGDIR) {
            if ($blast_file =~ /RefGblast\.parsed$/) {
                my $parsed = $full_path."/".$blast_file;
                open (IN, $parsed) or die "can not open file $parsed!\n";
                while (<IN>) {
                    if ($_ =~ /#/) { next; }
                    chomp;
                    my ($read_ID, $length, $category, $lineage, $hit_name, $e_value) = split("\t", $_);
                    $blastn_RefG{$ref_genome_choice}++;
                }
                close IN;
            }
        }
        closedir RefGDIR;
    }

    # RefGfiltered BLASTN
    if ($name =~ /RefGfiltered_BLASTN$/) {
        opendir (BNDIR, $full_path) or die "can not open dir $full_path!\n";
        foreach my $blast_file (readdir BNDIR) {
            if ($blast_file =~ /blastn\.parsed$/) {
                my $blast_out = $blast_file;
                $blast_out =~ s/\.blastn\.parsed/\.blastn\.out/;
                $blast_out = $full_path."/".$blast_out;

                my $blast_s = $blast_file;
                $blast_s =~ s/\.blastn\.parsed/\.blastn\.summary/;
                $blast_s = $full_path."/".$blast_s;

                push @blast_files_blastn, $blast_s;

                my $parsed = $full_path."/".$blast_file;
                &collect_information(
                    $parsed,
                    \%blastn,
                    \%viral_reads_blastn,
                    \%best_e_blastn,
                    \%lineage_blastn,
                    \%lineage_gi,
                    \%num_reads,
                    \@unassigned_reads,
                    \@ambiguous_reads
                );
            }
        }
        closedir BNDIR;
    }
}
closedir DH;

# get detailed information about each viral read
&get_viral_read_info(\@blast_files_blastn, \%blast_readinfo);

# print out report for this library
print OUT1 $dir, "\n";
printf OUT1 "%12s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\n", "category", "total", "BN_RefG", "BN", "BX_NR";

foreach my $key (sort {$a cmp $b } keys %blastx) {
    printf OUT1 "%12s\t%7d\t%7d\t%7d\t%7d\n",
        $key,
        $blastn_RefG{$key} + $blastn{$key} + $blastx{$key},
        $blastn_RefG{$key},
        $blastn{$key},
        $blastx{$key};
}

print OUT1 "\n###########################################################\n\n";

foreach my $gi (sort {$num_reads{$a} <=> $num_reads{$b}} keys %num_reads) {
    print OUT1 $gi, "\t", $lineage_gi{$gi}, "\ttotal number of reads: ", $num_reads{$gi},  "\n\n";
    print OUT1 "QueryName\tQuerylength\t         HitName       \tHitLen\t                             HitDesc                       \tAlnLen\t%ID\tHitStart\tHitEnd\te\n";

    if (defined $lineage_blastn{$gi}) {
        if (scalar @{$lineage_blastn{$gi}}) {
            print OUT1 "reads from blastn:\n";
            foreach my $read (sort {$a cmp $b} @{$lineage_blastn{$gi}}) {
                print OUT1 $blast_readinfo{$read};
            }
        }
    }

    print OUT1 "\n##################################################\n\n";
}

# get all the viral reads and put into output file:
foreach my $gi (keys %num_reads) {
    foreach my $read (@{$lineage_blastn{$gi}}) {
        print OUT2 ">$read\n";
        print OUT2 $seq{$read}, "\n";
    }
}

print OUT1 "# Finished Assignment Report\n";

exit;

#####################################################################################
# collect information from given directory
sub collect_information {
    my ($infile, $category_hash_ref, $viral_reads_hash_ref, $best_e_hash_ref,
        $lineage_hash_ref, $lineage_hash_gi, $num_reads_hash_ref,
        $unassigned_reads_arr_ref, $ambiguous_reads_arr_ref) = @_;

    open (IN, $infile) or die "can not open file $infile!\n";
    while (<IN>) {
        if ($_ =~ /#/) { next; }
        chomp;

        my ($read_ID, $length, $category, $lineage, $hit_name, $e_value) = split("\t", $_);

        my $gid = 0;
        if ($hit_name =~ /gi\|(\d+)\|/) {
            $gid = $1;
            $lineage_hash_gi->{$gid} = $lineage;
        }

        # NO Switch: increment category counter safely
        if (exists $category_hash_ref->{$category}) {
            $category_hash_ref->{$category}++;
        } else {
            $category_hash_ref->{"other"}++;
        }

        if (($category eq "Viruses") && $gid != 0) {
            $viral_reads_hash_ref->{$read_ID} = 1;
            $best_e_hash_ref->{$read_ID} = $e_value;

            if (!(defined $lineage_hash_ref->{$gid})) {
                $lineage_hash_ref->{$gid} = [$read_ID];
            } else {
                push @{$lineage_hash_ref->{$gid}}, $read_ID;
            }

            if (defined $num_reads_hash_ref->{$gid}) {
                $num_reads_hash_ref->{$gid}++;
            } else {
                $num_reads_hash_ref->{$gid} = 1;
            }

        } elsif ($category eq "Ambiguous") {
            push @{$ambiguous_reads_arr_ref}, $read_ID;

        } elsif ($category eq "unassigned") {
            push @{$unassigned_reads_arr_ref}, $read_ID;
        }
    }
    close IN;
}

############################################################################
sub read_FASTA_data {
    my $fastaFile = shift @_;

    my $oldseperator = $/;
    $/ = ">";

    my %fastaSeq;
    open (FAfile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";
    while (my $line = <FAfile>) {
        if ($line =~ /^\s*$/) {
            next;
        } elsif ($line =~ /^\s*#/) {
            next;
        } elsif ($line ne ">") {
            chomp $line;
            my @rows = split (/\n/, $line);
            my $temp = shift @rows;
            my @temp = split(/\s+/, $temp);
            my $name = shift @temp;

            my $Seq = join("", @rows);
            $Seq =~ s/\s//g;
            $fastaSeq{$name} = $Seq;
        }
    }

    $/ = $oldseperator;
    close FAfile;

    return %fastaSeq;
}

#############################################################################
# get detailed information about each viral read (reads *.blastn.summary)
sub get_viral_read_info {
    my ($report_file_ref, $blast_readinfo_hash_ref) = @_;

    foreach my $file (@{$report_file_ref}) {
        foreach my $line (`cat $file`) {
            if ($line =~ /Finished summary/) { next; }
            my @ss = split("\t", $line);
            $blast_readinfo_hash_ref->{$ss[0]} = $line;
        }
    }
}