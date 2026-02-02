 
#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;

#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m"; 
my $version = "1.2";
#usage information
(my $usage = <<OUT) =~ s/\t+//g;
This script will run the virus discovery pipeline on LSF cluster.
Pipeline version: $version
$yellow		Usage: perl $0 --sre --rdir --log --q --groupname --users --step

 $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$green		[1]  Run bwa
$red		[2, or <=22]  Split files for RepeatMasker
			[3 or <=23]  Submit RepeatMasker job array
$yellow		[4 or <=24]  Sequence Qulity Control
$green		[5 or <=25]  Split files for Blast Reference Genome
			[6 or <=26]  Submit Blast Reference Genome job array
			[7 or <=27]  Parse Reference Genome Blast result
$gray		[8 or <=28]  Pool and split files for BlastN
			[9 or <=29]  Submit BlastN job array
			[10 or <=30] Parse BlastN result
			[11 or <=31] Get summary of BlastN
$purple		[12 or <=32] Assignment report for each sample
			[13 or <=33] Assignment summary for each sample
			[14 or <=34] Generate report for the run
$normal
OUT


#__DEFAULT NUMBER OF BINS IE (MUST BE INTEGER)
my $step_number = -1;
my $status_rerun=0; 

#__HELP (BOOLEAN, DEFAULTS TO NO-HELP)
my $help = 0;
my $q_name="";
#__FILE NAME (STRING, NO DEFAULT)
my $run_dir="";
my $log_dir="";

my $compute_username="";
my $group_name="";

#__PARSE COMMAND LINE
my $status = &GetOptions (
      "step=i" => \$step_number,
      "sre=i" => \$status_rerun,	
      "groupname=s" => \$group_name,
      "users=s" => \$compute_username,	
      "rdir=s" => \$run_dir,
	  "log=s"  => \$log_dir,
	  "q=s" => \$q_name,
   	  "help" => \$help
	);
 
#print $status,"\n";

if ($help || $run_dir eq "" || $log_dir eq "" || $group_name eq "" || $compute_username eq "" || $step_number<0 ) {
	print "wrong option\n";
	print $usage;
    exit;
   }

# unless ($step_number >=0)&&(($step_number <= 17) || ($step_number >= 22));



#####################################################################################
# values need to be modified to adapt to local environment
my $email = "scao\@wustl\.edu";

# software path
#my $cd_hit = "/gscuser/mboolcha/software/cdhit/cd-hit-est";
# my RepeatMasker = "RepeatMasker";
# my blastn = "/gscuser/scao/tools/ncbi-blast+/bin/blastn";
# my $bwa="/storage1/fs1/songcao/Active/Software/anaconda3/bin/bwa"
#my $blastx = "/gscuser/scao/tools/software/ncbi-blast+/bin/blastx";

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# path and name of databases
#my $db_BN = "/gscuser/scao/gc3027/nt/nt";
#my $db_BX = "/gscuser/scao/gc3027/nr/nr";
#my $bwa_ref = "/gscuser/scao/gc3027/fasta/virus/virusdb_082414.fa";

my $db_BN = "/storage1/fs1/songcao/Active/Database/nt/nt";
my $db_BX = "/storage1/fs1/songcao/Active/Database/nr";
my $bwa_ref = "/storage1/fs1/songcao/Active/Database/nt012414_RE_Split/nt012414_virus_abbr_cdhit98.fa";

# reference genome taxonomy classification and database location.
# It's better to change $refrence_genome_taxonomy and $reference_genome based on the data being analyzed.
my $refrence_genome_taxonomy = "";
my $reference_genome = "";

#if ($ref_genome_choice == 1) {
#	$refrence_genome_taxonomy = "Homo"; # use Bacteria, Homo, Phage, Fungi, Mus, other

	# path to the reference genome	
#	$reference_genome = "/gscmnt/gc3027/dinglab/medseq/human70.37/humandnacdna.fa";
#}

$refrence_genome_taxonomy = "Homo";

$reference_genome = "/storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa";

#####################################################################################
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];

# To run jobs faster, split large fasta files to small ones. Split to specific number of 
# files instead of specific sequences in each small file, because the number of job array 
# cannot be determined if spliting to specific number of sequences in each file. Job 
# number is required by qsub ${SGE_TASK_ID}. The minimum size of each file is 4kb. 
# The number of files should be determined accourding to CPUs available in the computer
# cluster.

# The number of small fasta files to split to from a large file for RepeatMasker
my $file_number_of_RepeatMasker = 100; #default 
# the number of small fasta files to split to from a large file for Blast_Reference_Genome
my $file_number_of_Blast_Ref_Genome = 2; #default
# the number of small fasta files to split to from a large file for Blast_N
my $file_number_of_Blast_N = 2; #default
# the number of small fasta files to split to from a large file for Blast_X
#my $file_number_of_Blast_X = 200; #default

#store job files here

my $HOME1=$log_dir;

if (! -d $HOME1) {
    `mkdir $HOME1`;
}
#store job files here
if (! -d $HOME1."/tmp") {
    `mkdir $HOME1"/tmp"`;
}
my $job_files_dir = $HOME1."/tmp";

#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR") {
    `mkdir $HOME1"/LSF_DIR"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR";

# obtain script path
my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

my $hold_RM_job = "norm"; 
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";

#directory suffix constants
my $REPEAT_MASKER_DIR_SUFFIX = "fa.cdhit_out_RepeatMasker";
my $BLAST_RefG_DIR_SUFFIX = "fa.cdhit_out.masked.goodSeq_RefGblast";
my $BLAST_NT_DIR_SUFFIX = "RefGfiltered_BLASTN";
my $BLASTX_NR_DIR_SUFFIX = "BNFiltered_BLASTX_NR";

# get sample list in the run, name should not contain "."
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
# --- NEW: chaining helpers (per-sample LSF dependencies) ---
my $CHAIN_MODE = ($step_number == 0 || $step_number >= 22) ? 1 : 0;  # chain in "run-all" modes
my $last_job_id = undef;  # reset per sample before submitting the first job

sub submit_with_dep_cmd {
    my ($cmd) = @_;
    if ($CHAIN_MODE && defined $last_job_id) {
        my $dep = qq{ -w "done($last_job_id)"};
        if ($cmd =~ / -o /) {
            $cmd =~ s/ -o /$dep -o /;
        } else {
            $cmd .= $dep;
        }
    }
    print $cmd, "\n";
    my $out = `$cmd`;
    my ($jobid) = ($out =~ /Job\s+<(\d+)>/);
    if ($CHAIN_MODE && defined $jobid) {
        $last_job_id = $jobid;
    }
    return $jobid;
}

close DH;

# check to make sure the input directory has correct structure
# &check_input_dir($run_dir);

# start data processsing
if ($step_number < 14 || $step_number>=22) {
	#begin to process each sample
	for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_dir_list[$i];
        print "Processing sample: ", $sample_name, "\n";
    	if (!($sample_name =~ /\./)) {
			$sample_full_path = $run_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
                # NEW: reset chain per sample
                $last_job_id = undef;
				$current_job_file="";
				if ($step_number == 0 || $step_number>=22) {#run the whole pipeline
					######################################################################
					#cd-hit
                    if($step_number==0) 
					{ &bsub_bwa();}

					######################################################################
					#RepeatMasker
					#split file for RepeatMasker
					#my $f_fa=$sample_full_path.".fa"; 
					#if(! -f $f_fa) { next; }
	
					if($step_number<=22) 
					{
					 &split_for_RepeatMasker(); }
                                        
					#submit RepeatMasker job array
					if($step_number<=23) 
					{
					&submit_job_array_RM();
					$hold_RM_job=$current_job_file; # to limit number repeatmasker jobs run in the cluster at the same time. Can be removed if the cluster is able to handle the volumn of data input/output. 
                                        }
					######################################################################
					#Sequence Quality Control
					if($step_number<=24)
					{ &seq_QC();}

					######################################################################
					#BLASTn against Reference Genome
					if($step_number<=25) 
					{
					&split_for_blast_RefG();}

					#submit Blast RefG job array
					if($step_number<=26) 
					{
					&submit_job_array_blast_RefG();}

					if($step_number<=27)
					{
					#parser Blast RefG file
					&parse_blast_RefG();}


					######################################################################
					#BLASTn against nt
					#pool and split files for BLASTn
					if($step_number<=28)
					{
					&pool_split_for_blast_N();}

					#submit BLASTn job array
					if($step_number<=29)
					{
					&submit_job_array_blast_N();}

					#parser BLASTn output file
					if($step_number<=30)
					{
					 &parse_blast_N();}

                    if($step_number<=31)
                    {
                     &blast_S();}

					if($step_number<=32){
					&report_for_each_sample();}

					#Assignment summary for each sample
					if($step_number<=33) {
					&summary_for_each_sample();}
				
				######################################################################
				#run the pipeline step by step
				}elsif ($step_number == 1) {
					&bsub_bwa();
				}elsif ($step_number == 2) {
					&split_for_RepeatMasker(1);
				}elsif ($step_number == 3) {
					&submit_job_array_RM(1); 
					$hold_RM_job=$current_job_file; # to limit number of repeatmasker jobs
				}elsif ($step_number == 4) {
					&seq_QC(1);
				}elsif ($step_number == 5) {
					&split_for_blast_RefG(1);
				}elsif ($step_number == 6) {
					&submit_job_array_blast_RefG(1);
				}elsif ($step_number == 7) {
					&parse_blast_RefG(1);
				}elsif ($step_number == 8) {
					&pool_split_for_blast_N(1);
				}elsif ($step_number == 9) {
					&submit_job_array_blast_N(1);
				}elsif ($step_number == 10) {
					&parse_blast_N(1);
				}elsif ($step_number == 11) {
                                        &blast_S(1);
				}elsif ($step_number == 12) {
					 &report_for_each_sample(1);
				}elsif ($step_number == 13) {
					 &summary_for_each_sample(1);
				}
			}
		}
	}
}

##########################################################################################
# generate report for the run
if (($step_number == 0) || ($step_number == 14) || ($step_number>=22)) {

	# -----------------------------
# Submit job to generate report (Docker + bsub < script.sh so #BSUB lines still apply)
# -----------------------------

print $yellow, "Submitting jobs for generating the report for the run ....", $normal, "\n";

$hold_job_file    = $current_job_file;
$current_job_file = "Run_report_".$$.".sh";

my $sh_file  = "$job_files_dir/$current_job_file";
my $lsf_out  = "$lsf_file_dir/$current_job_file.out";
my $lsf_err  = "$lsf_file_dir/$current_job_file.err";

open(REPRUN, ">$sh_file") or die $!;
print REPRUN "#!/bin/bash\n";
print REPRUN "#BSUB -n 1\n";
print REPRUN "#BSUB -R \"rusage[mem=40000]\"\n";
print REPRUN "#BSUB -M 40000000\n";
# print REPRUN "#BSUB -q ding-lab\n";
print REPRUN "#BSUB -o $lsf_out\n";
print REPRUN "#BSUB -e $lsf_err\n";
print REPRUN "#BSUB -J $current_job_file\n";
print REPRUN "#BSUB -w \"$hold_job_file\"\n";

# Only once (you had it duplicated)
print REPRUN "BAD_SEQ=fa.cdhit_out.masked.badSeq\n"; # output of RepeatMasker
print REPRUN "OUTPUT=".$run_dir."/Analysis_Report_gi_".$working_name."\n";

print REPRUN "if [ -f \"\$OUTPUT\" ]\n";            # file exists?
print REPRUN "then\n";
print REPRUN "  grep \"# Finished\" \"\${OUTPUT}\"\n";
print REPRUN "  CHECK=\$?\n";
print REPRUN "  while [ \${CHECK} -eq 1 ]\n";      # grep unsuccessful, file not finished
print REPRUN "  do\n";
print REPRUN "    ".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version."\n";
print REPRUN "    grep \"# Finished\" \"\${OUTPUT}\"\n";
print REPRUN "    CHECK=\$?\n";
print REPRUN "  done\n";
print REPRUN "else\n";                              # file does not exist
print REPRUN "  ".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version."\n";
print REPRUN "  grep \"# Finished\" \"\${OUTPUT}\"\n";
print REPRUN "  CHECK=\$?\n";
print REPRUN "  while [ \${CHECK} -eq 1 ]\n";
print REPRUN "  do\n";
print REPRUN "    ".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version."\n";
print REPRUN "    grep \"# Finished\" \"\${OUTPUT}\"\n";
print REPRUN "    CHECK=\$?\n";
print REPRUN "  done\n";
print REPRUN "fi\n";

close REPRUN;

# Docker-enabled submission while keeping #BSUB directives in the script
$bsub_com =
  "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active " .
  "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" " .
  "bsub -a 'docker(scao/virusscan:0.0.2)' < $sh_file\n";

print $bsub_com;
submit_with_dep_cmd($bsub_com);

}


#######################################################################
if ($step_number == 0) {
	print $green, "All jobs are submitted! You will get email notification when this run is completed.\n",$normal;
}

exit;


########################################################################
# subroutines 

sub check_input_dir {
	my ($input_dir) = @_;
	my $have_input_sample = 0;
	
	# get sample list in the run, name should not contain "."
	opendir(DH, $input_dir) or die "Cannot open dir $input_dir: $!\n";
	my @sample_list = readdir DH;
	close DH;
	
	for (my $i=0;$i<@sample_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_list[$i];
		if (!($sample_name =~ /\./)&&!($sample_name =~/Analysis_/)) {
			$have_input_sample = 1;
			$sample_full_path = $input_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				my $input_file = $input_dir."/".$sample_name."/".$sample_name.".bam";
				if (!(-e $input_file)) { # input file does not exist
					print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
					die;
				}
			}
			else { # input sample directory does not exist
				print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
				die;
			}
		}
	}

	if (!($have_input_sample)) { # does not have any input sample directory
		print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
		die;
	}

}

########################################################################
sub bsub_bwa{

    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

    $current_job_file = "j1_bwa_".$sample_name.".sh";

    my $IN_bam = $sample_full_path."/".$sample_name.".bam";

    if (! -e $IN_bam) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam is empty!", $normal, "\n\n";
    }

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out;`;
    `rm $lsf_err;`;

    open(my $BWA, ">$job_files_dir/$current_job_file") or die $!;

    print $BWA "#!/bin/bash\n";
	print $BWA "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";
    print $BWA "BWA_IN=".$sample_full_path."/".$sample_name.".bam\n";
    print $BWA "BWA_fq=".$sample_full_path."/".$sample_name.".fq\n";
    print $BWA "BWA_sai=".$sample_full_path."/".$sample_name.".sai\n";
    #print BWA "BWA_sam=".$sample_full_path."/".$sample_name.".sam\n";
    #print BWA "BWA_bam=".$sample_full_path."/".$sample_name.".realign.bam\n";
    #print BWA "BWA_mapped_bam=".$sample_full_path."/".$sample_name.".mapped.bam\n";
    print $BWA "BWA_mapped=".$sample_full_path."/".$sample_name.".mapped.reads\n";
    print $BWA "BWA_fa=".$sample_full_path."/".$sample_name.".fa\n";
	#print BWA 
	print $BWA 'if [ ! -s $BWA_mapped ]',"\n";
    print $BWA "    then\n";
	print $BWA "rm \${BWA_sai}","\n";
	print $BWA "rm \${BWA_fq}","\n";
	#print BWA "mkfifo \${BWA_sai}","\n";
	print $BWA "mkfifo \${BWA_fq}","\n";
	#0x100: secondary alignment
	#0x800: supplementary alignment
    #H: Hard clipping
	#S: Soft clipping
	print $BWA "samtools view -h \${BWA_IN} | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^gi/))) { print \$line;}\' | samtools view -Sb - | bamtools convert -format fastq > \${BWA_fq} \&","\n";
    #print BWA "bwa aln $bwa_ref -b0 \${BWA_IN} > \${BWA_sai} \&","\n";	
    print $BWA "bwa aln $bwa_ref \${BWA_fq} > \${BWA_sai}","\n";
    print $BWA 'rm ${BWA_fq}',"\n";
	print $BWA "mkfifo \${BWA_fq}","\n";
    print $BWA "samtools view -h \${BWA_IN} | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^gi/))) { print \$line;}\' | samtools view -Sb - | bamtools convert -format fastq > \${BWA_fq} \&","\n";
	#print BWA "samtools view -h \${BWA_IN} | gawk \'{if (substr(\$1,1,1)==\"\@\" || (and(\$2,0x4) || and(\$2,0x8) )) print}\' | samtools view -Sb - | bamtools convert -format fastq > \${BWA_fq} \&","\n";
	print $BWA "bwa samse $bwa_ref \${BWA_sai} \${BWA_fq} | grep -v \@SQ | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); if(\$ss[2]=~/^gi/) { print \$line; }\' > \${BWA_mapped}","\n";
	print $BWA "     ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
    print $BWA " 	".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
	print $BWA 'rm ${BWA_sai}',"\n";
	print $BWA 'rm ${BWA_fq}',"\n";
	print $BWA "else\n";
    print $BWA "     ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
    print $BWA "     ".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
    print $BWA "   fi\n";
    close $BWA;

    my $sh_file = $job_files_dir."/".$current_job_file;

    # Submit with docker + mounts (following seq_QC style)
    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active " .
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" " .
        "bsub -g /$compute_username/$group_name -q $q_name " .
        "-J \"$current_job_file\" " .
        "-n 1 " .
        "-R \"span[hosts=1] select[mem>30000] rusage[mem=30000]\" " .
        "-M 30000000 " .
        "-a 'docker(scao/virusscan:0.0.2)' " .
        "-o $lsf_out -e $lsf_err " .
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);

}


#####################################################################################

sub split_for_RepeatMasker {
    # split file for RepeatMasker

    $current_job_file = "j2_".$sample_name."_RM_split".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";

    `rm -f $lsf_out`;
    `rm -f $lsf_err`;

    open(my $RMSPLIT, ">", "$job_files_dir/$current_job_file") or die $!;

    print $RMSPLIT "#!/bin/bash\n";
    #print $RMSPLIT "set -euo pipefail\n\n";
    print $RMSPLIT "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";

    # inputs/dirs
    print $RMSPLIT "RMSPLIT_IN=".$sample_full_path."/".$sample_name.".fa\n";
    print $RMSPLIT "RM_DIR=".$sample_full_path."/".$sample_name.".$REPEAT_MASKER_DIR_SUFFIX\n";
    print $RMSPLIT "SAMPLE_DIR=".$sample_full_path."\n\n";

    # create dir if missing
    print $RMSPLIT 'if [ ! -d "${RM_DIR}" ]',"\n";
    print $RMSPLIT "then\n";
    print $RMSPLIT '  mkdir -p "${RM_DIR}"',"\n";
    print $RMSPLIT "fi\n\n";

    # run split + check, retry if needed
    print $RMSPLIT $run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.
                   ".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".
                   $sample_name.".fa.cdhit_out_file\n";
    print $RMSPLIT $run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
    print $RMSPLIT "CHECK=\$?\n";
    print $RMSPLIT "while [ \${CHECK} -eq 10 ]\n";
    print $RMSPLIT "do\n";
    print $RMSPLIT "  ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.
                   ".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".
                   $sample_name.".fa.cdhit_out_file\n";
    print $RMSPLIT "  ".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
    print $RMSPLIT "  CHECK=\$?\n";
    print $RMSPLIT "done\n";

    close($RMSPLIT);

    my $sh_file = $job_files_dir."/".$current_job_file;

    # Submit with docker and mounts (following seq_QC style)
    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active " .
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" " .
        "bsub -g /$compute_username/$group_name -q $q_name " .
        "-J \"$current_job_file\" " .
        "-n 1 " .
        "-R \"span[hosts=1] select[mem>30000] rusage[mem=30000]\" " .
        "-M 30000000 " .
        "-a 'docker(scao/virusscan:0.0.2)' " .
        "-o $lsf_out -e $lsf_err " .
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}


sub submit_job_array_RM {

    $current_job_file = "j3_".$sample_name."_RM".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".%I.out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".%I.err";

    `rm -f $lsf_out`;
    `rm -f $lsf_err`;

    open(my $RM, ">", "$job_files_dir/$current_job_file") or die $!;

    print $RM "#!/bin/bash\n";
    #print $RM "set -euo pipefail\n\n";

    # Make conda env tools visible by name (bwa, RepeatMasker, python3, etc.)
    print $RM "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";

    # robust array index
    print $RM "JOBIDX=\${LSB_JOBINDEX:-1}\n\n";

    # inputs/outputs
    print $RM "RM_dir=".$sample_full_path."/".$sample_name.".$REPEAT_MASKER_DIR_SUFFIX\n";
    print $RM "RMIN=".'${RM_dir}'."/".$sample_name.".fa.cdhit_out_file".'${JOBIDX}'.".fa\n";
    print $RM "RMOUT=".'${RMIN}'.".masked\n";
    print $RM "FINAL_OUT=".'${RM_dir}'."/".$sample_name.".fa.cdhit_out_file".'${JOBIDX}'.".fa.masked\n\n";

    # print $RM 'if [ -s "$RMIN" ]',"\n";
    # print $RM "then\n";
    # print $RM '  if [ ! -s "$FINAL_OUT" ]',"\n";
    # print $RM "  then\n";
    # # Now RepeatMasker can be called without full path
    # print $RM '     RepeatMasker -pa 4 -species "homo sapiens" "$RMIN"',"\n";
    # print $RM "  fi\n\n";

	print $RM 'if [ -s "$RMIN" ]',"\n";
	print $RM "then\n";
	print $RM '  attempt=1',"\n";
	print $RM '  max_attempts=100',"\n";
	print $RM '  while [ $attempt -le $max_attempts ] && [ ! -s "$FINAL_OUT" ]',"\n";
	print $RM "  do\n";
	print $RM '    echo "Running RepeatMasker attempt $attempt on $RMIN" >&2',"\n";
	print $RM '    RepeatMasker -pa 4 -species "homo sapiens" "$RMIN"',"\n";
	print $RM '    attempt=$((attempt+1))',"\n";
	print $RM "  done\n\n";
	print $RM '  if [ ! -s "$FINAL_OUT" ]',"\n";
	print $RM "  then\n";
	print $RM '    echo "ERROR: RepeatMasker failed after $max_attempts attempts; $FINAL_OUT missing" >&2',"\n";
	print $RM "    exit 1\n";
	print $RM "  fi\n";
	print $RM "fi\n";

    # copy actual RepeatMasker masked output to your expected name
    print $RM '  if [ -s "$RMOUT" ]',"\n";
    print $RM "  then\n";
    print $RM '    cp -f "$RMOUT" "$FINAL_OUT"',"\n";
    print $RM "  else\n";
    print $RM '    echo "WARNING: $RMOUT not produced; copying input to $FINAL_OUT" >&2',"\n";
    print $RM '    cp -f "$RMIN" "$FINAL_OUT"',"\n";
    print $RM "  fi\n";
    print $RM "fi\n";

    close($RM);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    # Mount the WHOLE famdb dir so both dfam39_full.0.h5 and dfam39_full.7.h5 are visible
    my $famdb_host = "/storage1/fs1/dinglab/Active/Projects/scao/database/dfam39_famdb";
    my $famdb_cont = "/opt/conda/envs/virusscan/share/RepeatMasker/Libraries/famdb";

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active /storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active $famdb_host:$famdb_cont:ro\" " .
        "bsub -g /$compute_username/$group_name -q $q_name " .
        "-J \"$current_job_file\[1-$file_number_of_RepeatMasker\]\" " .
        "$dep" .
        "-n 1 " .
        "-R \"span[hosts=1] select[mem>10000] rusage[mem=10000]\" " .
        "-M 10000000 " .
        "-a 'docker(scao/virusscan:0.0.2)' " .
        "-o $lsf_out -e $lsf_err " .
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub seq_QC {
    my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    } else {
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j4_".$sample_name."_QC".".sh";

  	my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";

    `rm -f $lsf_out`;
    `rm -f $lsf_err`;

    open(my $QC, ">", "$job_files_dir/$current_job_file") or die $!;

    print $QC "#!/bin/bash\n";
   # print $QC "set -euo pipefail\n\n";

    print $QC "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";

    print $QC "SAMPLE_DIR=".$sample_full_path."\n";
    print $QC "QC_OUT=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq\n";
    print $QC "f_fa=".$sample_full_path."/".$sample_name.".fa\n\n";

    print $QC 'if [ ! -f "$QC_OUT" ] && [ -s "$f_fa" ]',"\n";
    print $QC "then\n";
    print $QC "    ".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
    print $QC "    ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
    print $QC "    CHECK=\$?\n";
    print $QC "    while [ \${CHECK} -eq 10 ]\n";
    print $QC "    do\n";
    print $QC "        ".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
    print $QC "        ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
    print $QC "        CHECK=\$?\n";
    print $QC "    done\n";
    print $QC "else\n";
    print $QC "    ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
    print $QC "    CHECK=\$?\n";
    print $QC "    while [ \${CHECK} -eq 10 ]\n";
    print $QC "    do\n";
    print $QC "        ".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
    print $QC "        ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
    print $QC "        CHECK=\$?\n";
    print $QC "    done\n";
    print $QC "fi\n";

    close($QC);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active " .
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" " .
        "bsub -g /$compute_username/$group_name -q $q_name " .
        "-J \"$current_job_file\" " .
        "$dep" .
        "-n 1 " .
        "-R \"span[hosts=1] select[mem>10000] rusage[mem=10000]\" " .
        "-M 10000000 " .
        "-a 'docker(scao/virusscan:0.0.2)' " .
        "-o $lsf_out -e $lsf_err " .
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);

}

#####################################################################################
sub split_for_blast_RefG {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j5_".$sample_name."_RefG_split".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $RefGS, ">", "$job_files_dir/$current_job_file") or die $!;
    print $RefGS "#!/bin/bash\n";
    #print $RefGS "set -euo pipefail\n\n";
    print $RefGS "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";

    print $RefGS "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
    print $RefGS "SAMPLE_DIR=".$sample_full_path."\n\n";

    print $RefGS 'if [ ! -d "$RefG_DIR" ]',"\n";
    print $RefGS "then\n";
    print $RefGS "  mkdir -p \"\${RefG_DIR}\"\n";
    print $RefGS "fi\n\n";

    print $RefGS $run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.
                 ".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".
                 $sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
    print $RefGS $run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
    print $RefGS "CHECK=\$?\n";
    print $RefGS "while [ \${CHECK} -eq 10 ]\n";
    print $RefGS "do\n";
    print $RefGS "  ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.
                 ".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".
                 $sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
    print $RefGS "  ".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
    print $RefGS "  CHECK=\$?\n";
    print $RefGS "done\n";
    close($RefGS);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>10000] rusage[mem=10000]\" ".
        "-M 10000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub submit_job_array_blast_RefG {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j6_".$sample_name."_BRefG".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".%I.out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".%I.err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $RefG, ">", "$job_files_dir/$current_job_file") or die $!;
    print $RefG "#!/bin/bash\n";
    #print $RefG "set -euo pipefail\n\n";
    print $RefG "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";
    print $RefG "JOBIDX=\${LSB_JOBINDEX:-1}\n\n";

    print $RefG "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
    print $RefG "BlastRefGOUT=".'${RefG_DIR}'."/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${JOBIDX}'.".RefGblast.out\n";
    print $RefG "QUERY=".'${RefG_DIR}'."/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${JOBIDX}'.".fa\n\n";

    print $RefG 'if [ -s "$QUERY" ]',"\n";
    print $RefG "then\n";
    print $RefG '  CHECK=1',"\n";
    print $RefG '  while [ $CHECK -ne 0 ]',"\n";
    print $RefG "  do\n";
    print $RefG "    blastn -evalue 1e-9 -show_gis -num_threads 4 -num_descriptions 2 -num_alignments 2 -query \"\${QUERY}\" -out \"\${BlastRefGOUT}\" -db $reference_genome\n";
    print $RefG '    tail -10 "${BlastRefGOUT}" | grep -q Matrix',"\n";
    print $RefG '    CHECK=$?',"\n";
    print $RefG "  done\n";
    print $RefG "fi\n";
    close($RefG);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\[1-$file_number_of_Blast_Ref_Genome\]\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>20000] rusage[mem=20000]\" ".
        "-M 20000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub parse_blast_RefG {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j7_".$sample_name."_PRefG".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".%I.out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".%I.err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $PRefG, ">", "$job_files_dir/$current_job_file") or die $!;
    print $PRefG "#!/bin/bash\n";
    #print $PRefG "set -euo pipefail\n\n";
    print $PRefG "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";
    print $PRefG "JOBIDX=\${LSB_JOBINDEX:-1}\n\n";

    print $PRefG "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
    print $PRefG "BlastRefGOUT=".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${JOBIDX}'.".RefGblast.out\n";
    print $PRefG "BlastRefGIN=".'${RefG_DIR}'."/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${JOBIDX}'.".fa\n";
    print $PRefG "PARSED=".'${RefG_DIR}'."/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${JOBIDX}'.".RefGblast.parsed\n\n";

    print $PRefG 'if [ -s "$BlastRefGIN" ]',"\n";
    print $PRefG "then\n";
    print $PRefG '  CHECK=1',"\n";
    print $PRefG '  while [ $CHECK -ne 0 ]',"\n";
    print $PRefG "  do\n";
    print $PRefG "    ".$run_script_path."BLASTn_RefGenome_parser.pl \"\${RefG_DIR}\" \"\${BlastRefGOUT}\" $refrence_genome_taxonomy\n";
    print $PRefG '    tail -5 "$PARSED" | grep -q Summary',"\n";
    print $PRefG '    CHECK=$?',"\n";
    print $PRefG "  done\n";
    print $PRefG "fi\n";
    close($PRefG);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\[1-$file_number_of_Blast_Ref_Genome\]\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>10000] rusage[mem=10000]\" ".
        "-M 10000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub pool_split_for_blast_N {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j8_".$sample_name."_BN_split".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $BNS, ">", "$job_files_dir/$current_job_file") or die $!;
    print $BNS "#!/bin/bash\n";
    #print $BNS "set -euo pipefail\n\n";
    print $BNS "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";

    print $BNS "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
    print $BNS "SAMPLE_DIR=".$sample_full_path."\n";
    print $BNS "RefGFiltered_fa=".$sample_full_path."/".$sample_name.".RefGfiltered.fa\n";
    print $BNS "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n\n";

    print $BNS 'mkdir -p "$BN_DIR"',"\n";
    print $BNS 'rm -f "$RefGFiltered_fa" || true',"\n";
    print $BNS 'cat ${RefG_DIR}/*.RefGfiltered.fa >> "$RefGFiltered_fa"',"\n";

    print $BNS $run_script_path."check_split_BN.pl \${SAMPLE_DIR}\n";
    print $BNS "CHECK=\$?\n";
    print $BNS "echo \"[BN_split] initial CHECK=\${CHECK}\" \n";   # <-- add this
    print $BNS "while [ \${CHECK} -eq 10 ]\n";
    print $BNS "do\n";
    print $BNS "  ".$run_script_path."split_fasta.pl -i \"\${RefGFiltered_fa}\" -o \"\${BN_DIR}\" -n $file_number_of_Blast_N -p ".$sample_name.".RefGfiltered.fa_file\n";
    print $BNS "  ".$run_script_path."check_split_BN.pl \${SAMPLE_DIR}\n";
    print $BNS "  CHECK=\$?\n";
    print $BNS "done\n";
    close($BNS);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>10000] rusage[mem=10000]\" ".
        "-M 10000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub submit_job_array_blast_N {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j9_".$sample_name."_BN".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".%I.out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".%I.err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $BN, ">", "$job_files_dir/$current_job_file") or die $!;
    print $BN "#!/bin/bash\n";
   # print $BN "set -euo pipefail\n\n";
    print $BN "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";
    print $BN "JOBIDX=\${LSB_JOBINDEX:-1}\n\n";

    print $BN "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
    print $BN "BlastNOUT=".'${BN_DIR}'."/".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".blastn.out\n";
    print $BN "QUERY=".'${BN_DIR}'."/".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".fa\n\n";

    print $BN 'if [ -s "$QUERY" ]',"\n";
    print $BN "then\n";
    print $BN '  CHECK1=1',"\n";
    print $BN '  CHECK2=0',"\n";
    print $BN '  while [ ${CHECK1} -ne 0 ] || [ ${CHECK2} -eq 0 ]',"\n";
    print $BN "  do\n";
    print $BN "    blastn -evalue 1e-9 -show_gis -num_threads 4 -query \"\${QUERY}\" -out \"\${BlastNOUT}\" -db $db_BN\n";
    print $BN '    tail -5 "${BlastNOUT}" | grep -q Matrix',"\n";
    print $BN '    CHECK1=$?',"\n";
    print $BN '    grep -q "no longer exists in database" "${BlastNOUT}"',"\n";
    print $BN '    if [ $? -eq 0 ]; then CHECK2=0; else CHECK2=1; fi',"\n";
    print $BN "  done\n";
    print $BN "fi\n";
    close($BN);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\[1-$file_number_of_Blast_N\]\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>40000] rusage[mem=40000]\" ".
        "-M 40000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub parse_blast_N {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j10_".$sample_name."_PBN".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".%I.out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".%I.err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $PBN, ">", "$job_files_dir/$current_job_file") or die $!;
    print $PBN "#!/bin/bash\n";
    #print $PBN "set -euo pipefail\n\n";
    print $PBN "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";
    print $PBN "JOBIDX=\${LSB_JOBINDEX:-1}\n\n";

    print $PBN "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
    print $PBN "BlastNOUT=".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".blastn.out\n";
    print $PBN "BlastNIN=".'${BN_DIR}'."/".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".fa\n";
    print $PBN "PARSED=".'${BN_DIR}'."/".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".blastn.parsed\n\n";

    print $PBN 'if [ -s "$BlastNIN" ]',"\n";
    print $PBN "then\n";
    print $PBN '  CHECK=10',"\n";
    print $PBN '  while [ ${CHECK} -eq 10 ]',"\n";
    print $PBN "  do\n";
    print $PBN "    ".$run_script_path."BLASTn_NT_parser.v2.pl ".$log_dir." ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \"\${BlastNOUT}\"\n";
    print $PBN "    ".$run_script_path."check_Blast_parsed_file.pl \"\${PARSED}\"\n";
    print $PBN "    CHECK=\$?\n";
    print $PBN "  done\n";
    print $PBN "fi\n";
    close($PBN);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\[1-$file_number_of_Blast_N\]\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>10000] rusage[mem=10000]\" ".
        "-M 10000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub blast_S {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j11_".$sample_name."_blastS".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".%I.out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".%I.err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $PS, ">", "$job_files_dir/$current_job_file") or die $!;
    print $PS "#!/bin/bash\n";
    #print $PS "set -euo pipefail\n\n";
    print $PS "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";
    print $PS "JOBIDX=\${LSB_JOBINDEX:-1}\n\n";

    print $PS "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
    print $PS "BlastNparsed=".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".blastn.parsed\n";
    print $PS "BlastNIN=".'${BN_DIR}'."/".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".fa\n";
    print $PS "OUTPUT=".'${BN_DIR}'."/".$sample_name.".RefGfiltered.fa_file".'${JOBIDX}'.".blastn.summary\n\n";

    print $PS 'if [ -s "$BlastNIN" ]',"\n";
    print $PS "then\n";
    print $PS '  CHECK=1',"\n";
    print $PS '  while [ $CHECK -ne 0 ]',"\n";
    print $PS "  do\n";
    print $PS "    ".$run_script_path."blast_summary.v2.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \"\${BlastNparsed}\"\n";
    print $PS '    grep -q "Finished summary" "${OUTPUT}"',"\n";
    print $PS '    CHECK=$?',"\n";
    print $PS "  done\n";
    print $PS "fi\n";
    close($PS);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\[1-$file_number_of_Blast_N\]\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>10000] rusage[mem=10000]\" ".
        "-M 10000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}


#####################################################################################

sub report_for_each_sample {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j12_".$sample_name."_Rep".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $REP, ">", "$job_files_dir/$current_job_file") or die $!;
    print $REP "#!/bin/bash\n";
    #print $REP "set -euo pipefail\n\n";
    print $REP "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";

    print $REP "INPUT=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq\n";
    print $REP "REPORT=".$sample_full_path."/".$sample_name.".gi.AssignmentReport\n\n";

    print $REP 'CHECK=1',"\n";
    print $REP 'while [ ${CHECK} -ne 0 ]',"\n";
    print $REP "do\n";
    print $REP "  ".$run_script_path."assignment_report_virus_gi.v2.pl ".$sample_full_path." \"\${INPUT}\" $refrence_genome_taxonomy\n";
    print $REP '  grep -q "# Finished Assignment Report" "${REPORT}"',"\n";
    print $REP '  CHECK=$?',"\n";
    print $REP "done\n";
    close($REP);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>40000] rusage[mem=40000]\" ".
        "-M 40000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}

#####################################################################################

sub summary_for_each_sample {
    my ($step_by_step) = @_;
    if ($step_by_step) { $hold_job_file = ""; }
    else { $hold_job_file = $current_job_file; }

    $current_job_file = "j13_".$sample_name."_Sum".".sh";

    my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";
    `rm -f $lsf_out`; `rm -f $lsf_err`;

    open(my $SUM, ">", "$job_files_dir/$current_job_file") or die $!;
    print $SUM "#!/bin/bash\n";
    #print $SUM "set -euo pipefail\n\n";
    print $SUM "export PATH=/opt/conda/envs/virusscan/bin:\$PATH\n\n";

    print $SUM "OUTPUT=".$sample_full_path."/".$sample_name.".gi.AssignmentSummary\n";
    print $SUM "BAD_SEQ=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.badSeq\n\n";

    print $SUM 'CHECK=1',"\n";
    print $SUM 'while [ ${CHECK} -ne 0 ]',"\n";
    print $SUM "do\n";
    print $SUM "  ".$run_script_path."assignment_summary_gi.pl ".$sample_full_path." \"\${BAD_SEQ}\"\n";
    print $SUM '  grep -q "# Finished Assignment Summary" "${OUTPUT}"',"\n";
    print $SUM '  CHECK=$?',"\n";
    print $SUM "done\n";
    close($SUM);

    my $sh_file = "$job_files_dir/$current_job_file";

    my $dep = "";
    if (defined $hold_job_file && $hold_job_file ne "") {
        $dep = "-w \"$hold_job_file\" ";
    }

    $bsub_com =
        "LSF_DOCKER_VOLUMES=\"/storage1/fs1/dinglab/Active:/storage1/fs1/dinglab/Active ".
        "/storage1/fs1/songcao/Active:/storage1/fs1/songcao/Active\" ".
        "bsub -g /$compute_username/$group_name -q $q_name ".
        "-J \"$current_job_file\" $dep ".
        "-n 1 ".
        "-R \"span[hosts=1] select[mem>40000] rusage[mem=40000]\" ".
        "-M 40000000 ".
        "-a 'docker(scao/virusscan:0.0.2)' ".
        "-o $lsf_out -e $lsf_err ".
        "bash $sh_file\n";

    print $bsub_com;
    submit_with_dep_cmd($bsub_com);
}
