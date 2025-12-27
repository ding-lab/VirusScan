#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#color code
my $red    = "\e[31m";
my $gray   = "\e[37m";
my $yellow = "\e[33m";
my $green  = "\e[32m";
my $purple = "\e[35m";
my $cyan   = "\e[36m";
my $normal = "\e[0m";
my $version = "1.1";

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

# args
my $step_number = -1;
my $status_rerun = 0;
my $help = 0;
my $q_name = "";
my $run_dir = "";
my $log_dir = "";
my $compute_username = "";
my $group_name = "";

my $status = &GetOptions (
  "step=i"      => \$step_number,
  "sre=i"       => \$status_rerun,
  "groupname=s" => \$group_name,
  "users=s"     => \$compute_username,
  "rdir=s"      => \$run_dir,
  "log=s"       => \$log_dir,
  "q=s"         => \$q_name,
  "help"        => \$help
);

if ($help || $run_dir eq "" || $log_dir eq "" || $group_name eq "" || $compute_username eq "" || $step_number < 0) {
  print "wrong option\n";
  print $usage;
  exit;
}

#####################################################################################
# values need to be modified to adapt to local environment
my $email = "scao\@wustl\.edu";

# path and name of databases
my $db_BN  = "/storage1/fs1/songcao/Active/Database/nt";
my $db_BX  = "/storage1/fs1/songcao/Active/Database/nr";
my $bwa_ref = "/storage1/fs1/songcao/Active/Database/nt012414_RE_Split/nt012414_virus_abbr_cdhit98.fa";

# reference genome taxonomy classification and database location.
my $refrence_genome_taxonomy = "Homo";
my $reference_genome = "/storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa";

#####################################################################################
# docker image used by ALL steps (match step 1/2)
my $DOCKER_IMAGE = "scao/virusscan:0.0.2";

#####################################################################################
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name = (split(/\//, $run_dir))[-2];

my $file_number_of_RepeatMasker = 100;
my $file_number_of_Blast_Ref_Genome = 100;
my $file_number_of_Blast_N = 100;

#store job files here
my $HOME1 = $log_dir;

if (! -d $HOME1) {
  `mkdir $HOME1`;
}
if (! -d $HOME1."/tmp") {
  `mkdir $HOME1."/tmp"`;
}
my $job_files_dir = $HOME1."/tmp";

if (! -d $HOME1."/LSF_DIR") {
  `mkdir $HOME1."/LSF_DIR"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR";

# obtain script path
my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

my $hold_RM_job = "norm";
my $current_job_file = "";
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";

# directory suffix constants
my $REPEAT_MASKER_DIR_SUFFIX = "fa.cdhit_out_RepeatMasker";
my $BLAST_RefG_DIR_SUFFIX    = "fa.cdhit_out.masked.goodSeq_RefGblast";
my $BLAST_NT_DIR_SUFFIX      = "RefGfiltered_BLASTN";
my $BLASTX_NR_DIR_SUFFIX     = "BNFiltered_BLASTX_NR";

########################################################################
# Standardized LSF submission (make ALL steps look like step 1/2)
# - Job scripts are plain bash
# - All resource / queue / group / docker / out/err are set on the bsub command line
sub bsub_docker {
  my (%p) = @_;

  my $sh_file  = $p{sh_file}  or die "bsub_docker: missing sh_file\n";
  my $job_name = $p{job_name} || "job_$$";
  my $queue    = $p{queue}    || $q_name;
  my $group    = $p{group}    || "/$compute_username/$group_name";
  my $cpus     = $p{cpus}     || 1;
  my $mem_mb   = $p{mem_mb}   || 10000; # MB
  my $out      = $p{out}      || "$lsf_file_dir/$job_name.out";
  my $err      = $p{err}      || "$lsf_file_dir/$job_name.err";
  my $hold     = defined($p{hold}) ? $p{hold} : "";
  my $array_n  = $p{array_n}  || 0;
  my $extra_R  = $p{extra_R}  || "";

  my $R = "select[mem>$mem_mb] rusage[mem=$mem_mb]";
  $R .= " $extra_R" if $extra_R ne "";

  my $J = ($array_n && $array_n > 0) ? "$job_name\\[1-$array_n\\]" : $job_name;

  # LSF -M is KB (keep your existing convention: MB * 1000)
  my $M_kb = $mem_mb * 1000;

  my $cmd = "bsub -g $group -q $queue -n $cpus -R \"$R\" -M $M_kb"
          . " -a 'docker($DOCKER_IMAGE)' -o $out -e $err -J $J";
  $cmd .= " -w \"$hold\"" if ($hold ne "");
  $cmd .= " bash $sh_file\n";

  print $cmd;
  system($cmd) == 0 or die "bsub failed for $job_name\n";
}

########################################################################
# get sample list in the run, name should not contain "."
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check input structure
&check_input_dir($run_dir);

# start data processing
if ($step_number < 14 || $step_number >= 22) {
  for (my $i = 0; $i < @sample_dir_list; $i++) {
    $sample_name = $sample_dir_list[$i];
    if (!($sample_name =~ /\./)) {
      $sample_full_path = $run_dir."/".$sample_name;
      if (-d $sample_full_path) {
        print $yellow, "\nSubmitting jobs for the sample ", $sample_name, "...", $normal, "\n";
        $current_job_file = "";

        if ($step_number == 0 || $step_number >= 22) {
          if ($step_number == 0) { &bsub_bwa(); }

          if ($step_number <= 22) { &split_for_RepeatMasker(); }
          if ($step_number <= 23) { &submit_job_array_RM(); $hold_RM_job = $current_job_file; }

          if ($step_number <= 24) { &seq_QC(); }

          if ($step_number <= 25) { &split_for_blast_RefG(); }
          if ($step_number <= 26) { &submit_job_array_blast_RefG(); }
          if ($step_number <= 27) { &parse_blast_RefG(); }

          if ($step_number <= 28) { &pool_split_for_blast_N(); }
          if ($step_number <= 29) { &submit_job_array_blast_N(); }
          if ($step_number <= 30) { &parse_blast_N(); }
          if ($step_number <= 31) { &blast_S(); }

          if ($step_number <= 32) { &report_for_each_sample(); }
          if ($step_number <= 33) { &summary_for_each_sample(); }

        } elsif ($step_number == 1) {
          &bsub_bwa();
        } elsif ($step_number == 2) {
          &split_for_RepeatMasker(1);
        } elsif ($step_number == 3) {
          &submit_job_array_RM(1); $hold_RM_job = $current_job_file;
        } elsif ($step_number == 4) {
          &seq_QC(1);
        } elsif ($step_number == 5) {
          &split_for_blast_RefG(1);
        } elsif ($step_number == 6) {
          &submit_job_array_blast_RefG(1);
        } elsif ($step_number == 7) {
          &parse_blast_RefG(1);
        } elsif ($step_number == 8) {
          &pool_split_for_blast_N(1);
        } elsif ($step_number == 9) {
          &submit_job_array_blast_N(1);
        } elsif ($step_number == 10) {
          &parse_blast_N(1);
        } elsif ($step_number == 11) {
          &blast_S(1);
        } elsif ($step_number == 12) {
          &report_for_each_sample(1);
        } elsif ($step_number == 13) {
          &summary_for_each_sample(1);
        }
      }
    }
  }
}

##########################################################################################
# generate report for the run
if (($step_number == 0) || ($step_number == 14) || ($step_number >= 22)) {
  print $yellow, "Submitting jobs for generating the report for the run ....", $normal, "\n";
  $hold_job_file = $current_job_file;
  $current_job_file = "Run_report_".$$.".sh";

  open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
  print REPRUN "#!/bin/bash\n";
  print REPRUN "BAD_SEQ=fa.cdhit_out.masked.badSeq\n";
  print REPRUN "OUTPUT=".$run_dir."/Analysis_Report_gi_".$working_name."\n";
  print REPRUN 'if [ -f $OUTPUT ]',"\n";
  print REPRUN "then\n";
  print REPRUN '  grep "# Finished" ${OUTPUT}',"\n";
  print REPRUN '  CHECK=$?',"\n";
  print REPRUN '  while [ ${CHECK} -eq 1 ]',"\n";
  print REPRUN "  do\n";
  print REPRUN "    ".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version,"\n";
  print REPRUN '    grep "# Finished" ${OUTPUT}',"\n";
  print REPRUN '    CHECK=$?',"\n";
  print REPRUN "  done\n";
  print REPRUN "else\n";
  print REPRUN "  ".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version,"\n";
  print REPRUN '  grep "# Finished" ${OUTPUT}',"\n";
  print REPRUN '  CHECK=$?',"\n";
  print REPRUN '  while [ ${CHECK} -eq 1 ]',"\n";
  print REPRUN "  do\n";
  print REPRUN "    ".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version,"\n";
  print REPRUN '    grep "# Finished" ${OUTPUT}',"\n";
  print REPRUN '    CHECK=$?',"\n";
  print REPRUN "  done\n";
  print REPRUN "fi\n";
  close REPRUN;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 40000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file
  );
}

#######################################################################
# send email to notify the finish of the analysis
if (($step_number == 0) || ($step_number == 15) || ($step_number >= 22)) {
  print $yellow, "Submitting the job for sending an email when the run finishes ", $sample_name, "...", $normal, "\n";
  $hold_job_file = $current_job_file;
  $current_job_file = "Email_run_".$$.".sh";

  open(EMAIL, ">$job_files_dir/$current_job_file") or die $!;
  print EMAIL "#!/bin/bash\n";
  print EMAIL $run_script_path."send_email.pl ".$run_dir." ".$email."\n";
  close EMAIL;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 1000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file
  );
}

if ($step_number == 0) {
  print $green, "All jobs are submitted! You will get email notification when this run is completed.\n", $normal;
}
exit;

########################################################################
# subroutines

sub check_input_dir {
  my ($input_dir) = @_;
  my $have_input_sample = 0;

  opendir(DH, $input_dir) or die "Cannot open dir $input_dir: $!\n";
  my @sample_list = readdir DH;
  close DH;

  for (my $i = 0; $i < @sample_list; $i++) {
    $sample_name = $sample_list[$i];
    if (!($sample_name =~ /\./) && !($sample_name =~ /Analysis_/)) {
      $have_input_sample = 1;
      $sample_full_path = $input_dir."/".$sample_name;
      if (-d $sample_full_path) {
        my $input_file = $input_dir."/".$sample_name."/".$sample_name.".bam";
        if (!(-e $input_file)) {
          print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
          die;
        }
      } else {
        print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
        die;
      }
    }
  }

  if (!($have_input_sample)) {
    print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
    die;
  }
}

########################################################################
sub bsub_bwa {
  $current_job_file = "j1_bwa_".$sample_name.".".$$.".sh";

  my $IN_bam = $sample_full_path."/".$sample_name.".bam";
  if (! -e $IN_bam) {
    print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
    print "Warning: Died because there is no input bam file for bwa:\n";
    print "File $IN_bam does not exist!\n";
    die "Please check command line argument!", $normal, "\n\n";
  }
  if (! -s $IN_bam) {
    print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
    die "Warning: Died because $IN_bam is empty!", $normal, "\n\n";
  }

  my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
  my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";
  `rm $lsf_out`;
  `rm $lsf_err`;

  open(BWA, ">$job_files_dir/$current_job_file") or die $!;
  print BWA "#!/bin/bash\n";
  print BWA "BWA_IN=".$sample_full_path."/".$sample_name.".bam\n";
  print BWA "BWA_fq=".$sample_full_path."/".$sample_name.".fq\n";
  print BWA "BWA_sai=".$sample_full_path."/".$sample_name.".sai\n";
  print BWA "BWA_mapped=".$sample_full_path."/".$sample_name.".mapped.reads\n";
  print BWA "BWA_fa=".$sample_full_path."/".$sample_name.".fa\n";
  print BWA 'if [ ! -s $BWA_mapped ]',"\n";
  print BWA "then\n";
  print BWA "  rm -f \${BWA_sai}\n";
  print BWA "  rm -f \${BWA_fq}\n";
  print BWA "  mkfifo \${BWA_fq}\n";
  print BWA "  samtools view -h \${BWA_IN} | perl -ne '\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^gi/))) { print \$line;}' | samtools view -Sb - | bamtools convert -format fastq > \${BWA_fq} &\n";
  print BWA "  bwa aln $bwa_ref \${BWA_fq} > \${BWA_sai}\n";
  print BWA "  rm -f \${BWA_fq}\n";
  print BWA "  mkfifo \${BWA_fq}\n";
  print BWA "  samtools view -h \${BWA_IN} | perl -ne '\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^gi/))) { print \$line;}' | samtools view -Sb - | bamtools convert -format fastq > \${BWA_fq} &\n";
  print BWA "  bwa samse $bwa_ref \${BWA_sai} \${BWA_fq} | grep -v \\@SQ | perl -ne '\$line=\$_; \@ss=split(\"\\t\",\$line); if(\$ss[2]=~/^gi/) { print \$line; }' > \${BWA_mapped}\n";
  print BWA "  ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
  print BWA "  ".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
  print BWA "  rm -f \${BWA_sai}\n";
  print BWA "  rm -f \${BWA_fq}\n";
  print BWA "else\n";
  print BWA "  ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
  print BWA "  ".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
  print BWA "fi\n";
  close BWA;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 30000,
    out      => $lsf_out,
    err      => $lsf_err
  );
}

#####################################################################################
sub split_for_RepeatMasker {
  $current_job_file = "j2_".$sample_name."_RM_split_".$$.".sh";

  my $lsf_out = $lsf_file_dir."/".$current_job_file.".out";
  my $lsf_err = $lsf_file_dir."/".$current_job_file.".err";
  `rm $lsf_out`;
  `rm $lsf_err`;

  open(RMSPLIT, ">$job_files_dir/$current_job_file") or die $!;
  print RMSPLIT "#!/bin/bash\n";
  print RMSPLIT "RMSPLIT_IN=".$sample_full_path."/".$sample_name.".fa\n";
  print RMSPLIT "RM_DIR=".$sample_full_path."/".$sample_name.".$REPEAT_MASKER_DIR_SUFFIX\n";
  print RMSPLIT "SAMPLE_DIR=".$sample_full_path."\n\n";
  print RMSPLIT "if [ ! -d \${RM_DIR} ]\n";
  print RMSPLIT "then\n";
  print RMSPLIT "  mkdir \${RM_DIR}\n";
  print RMSPLIT "  ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".$sample_name.".fa.cdhit_out_file\n";
  print RMSPLIT "  ".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
  print RMSPLIT "  CHECK=$?\n";
  print RMSPLIT "  while [ \${CHECK} -eq 10 ]\n";
  print RMSPLIT "  do\n";
  print RMSPLIT "    ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".$sample_name.".fa.cdhit_out_file\n";
  print RMSPLIT "    ".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
  print RMSPLIT "    CHECK=$?\n";
  print RMSPLIT "  done\n";
  print RMSPLIT "else\n";
  print RMSPLIT "  ".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
  print RMSPLIT "  CHECK=$?\n";
  print RMSPLIT "  while [ \${CHECK} -eq 10 ]\n";
  print RMSPLIT "  do\n";
  print RMSPLIT "    ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".$sample_name.".fa.cdhit_out_file\n";
  print RMSPLIT "    ".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
  print RMSPLIT "    CHECK=$?\n";
  print RMSPLIT "  done\n";
  print RMSPLIT "fi\n";
  close RMSPLIT;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 30000,
    out      => $lsf_out,
    err      => $lsf_err
  );
}

#####################################################################################
sub submit_job_array_RM {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j3_".$sample_name."_RM_".$$.".sh";

  open(RM, ">$job_files_dir/$current_job_file") or die $!;
  print RM "#!/bin/bash\n";
  print RM "RM_IN=".$sample_full_path."/".$sample_name.".fa\n";
  print RM "RM_dir=".$sample_full_path."/".$sample_name.".$REPEAT_MASKER_DIR_SUFFIX\n";
  print RM "RMOUT=\${RM_dir}/".$sample_name.".fa.cdhit_out_file\${LSB_JOBINDEX}.fa.masked\n";
  print RM "RMIN=\${RM_dir}/".$sample_name.".fa.cdhit_out_file\${LSB_JOBINDEX}.fa\n";
  print RM "RMOTHER=\${RM_dir}/".$sample_name.".fa.cdhit_out_file\${LSB_JOBINDEX}.fa.out\n\n";
  print RM 'if [ -f $RMIN ]',"\n";
  print RM "then\n";
  print RM '  if [ ! -s $RMOUT ]',"\n";
  print RM "  then\n";
  print RM "    RepeatMasker -pa 4 \$RMIN\n";
  print RM "  fi\n\n";
  print RM '  if [ ! -f $RMOTHER ]',"\n";
  print RM "  then\n";
  print RM '    while [ ! -f $RMOTHER ]',"\n";
  print RM "    do\n";
  print RM "      RepeatMasker -pa 4 \$RMIN\n";
  print RM "    done\n";
  print RM "  fi\n\n";
  print RM '  if [ ! -f $RMOUT ]',"\n";
  print RM "  then\n";
  print RM "    cp \${RMIN} \${RMOUT}\n";
  print RM "  fi\n";
  print RM "fi\n";
  close RM;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 30000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file,
    array_n  => $file_number_of_RepeatMasker,
    extra_R  => "span[hosts=1]"
  );
}

#####################################################################################
sub seq_QC {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j4_".$sample_name."_QC_".$$.".sh";
  open(QC, ">$job_files_dir/$current_job_file") or die $!;
  print QC "#!/bin/bash\n";
  print QC "SAMPLE_DIR=".$sample_full_path."\n";
  print QC "QC_OUT=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq\n";
  print QC "f_fa=".$sample_full_path."/".$sample_name.".fa\n\n";
  print QC 'if [ ! -f $QC_OUT ] && [ -s $f_fa ]',"\n";
  print QC "then\n";
  print QC "  ".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
  print QC "  ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
  print QC "  CHECK=$?\n";
  print QC "  while [ \${CHECK} -eq 10 ]\n";
  print QC "  do\n";
  print QC "    ".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
  print QC "    ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
  print QC "    CHECK=$?\n";
  print QC "  done\n";
  print QC "else\n";
  print QC "  ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
  print QC "  CHECK=$?\n";
  print QC "  while [ \${CHECK} -eq 10 ]\n";
  print QC "  do\n";
  print QC "    ".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
  print QC "    ".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
  print QC "    CHECK=$?\n";
  print QC "  done\n";
  print QC "fi\n";
  close QC;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 10000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file
  );
}

#####################################################################################
sub split_for_blast_RefG {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j5_".$sample_name."_RefG_split_".$$.".sh";
  open(RefGS, ">$job_files_dir/$current_job_file") or die $!;
  print RefGS "#!/bin/bash\n";
  print RefGS "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
  print RefGS "SAMPLE_DIR=".$sample_full_path."\n\n";
  print RefGS 'if [ ! -d $RefG_DIR ]',"\n";
  print RefGS "then\n";
  print RefGS "  mkdir \${RefG_DIR}\n";
  print RefGS "  ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".$sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
  print RefGS "  ".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
  print RefGS "  CHECK=$?\n";
  print RefGS "  while [ \${CHECK} -eq 10 ]\n";
  print RefGS "  do\n";
  print RefGS "    ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".$sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
  print RefGS "    ".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
  print RefGS "    CHECK=$?\n";
  print RefGS "  done\n";
  print RefGS "else\n";
  print RefGS "  ".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
  print RefGS "  CHECK=$?\n";
  print RefGS "  while [ \${CHECK} -eq 10 ]\n";
  print RefGS "  do\n";
  print RefGS "    ".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".$sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
  print RefGS "    ".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
  print RefGS "    CHECK=$?\n";
  print RefGS "  done\n";
  print RefGS "fi\n";
  close RefGS;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 10000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file
  );
}

#####################################################################################
sub submit_job_array_blast_RefG {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j6_".$sample_name."_BRefG_".$$.".sh";
  open(RefG, ">$job_files_dir/$current_job_file") or die $!;
  print RefG "#!/bin/bash\n";
  print RefG "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
  print RefG "BlastRefGOUT=\${RefG_DIR}/".$sample_name.".fa.cdhit_out.masked.goodSeq_file\${LSB_JOBINDEX}.RefGblast.out\n";
  print RefG "QUERY=\${RefG_DIR}/".$sample_name.".fa.cdhit_out.masked.goodSeq_file\${LSB_JOBINDEX}.fa\n\n";
  print RefG 'if [ -s $QUERY ]',"\n";
  print RefG "then\n";
  print RefG '  if [ ! -f $BlastRefGOUT ]',"\n";
  print RefG "  then\n";
  print RefG "    blastn -evalue 1e-9 -show_gis -num_threads 4 -num_descriptions 2 -num_alignments 2 -query \${QUERY} -out \${BlastRefGOUT} -db $reference_genome\n";
  print RefG '    tail -10 ${BlastRefGOUT} | grep Matrix',"\n";
  print RefG "    CHECK=$?\n";
  print RefG '    while [ ${CHECK} -eq 1 ]',"\n";
  print RefG "    do\n";
  print RefG "      blastn -evalue 1e-9 -show_gis -num_threads 4 -num_descriptions 2 -num_alignments 2 -query \${QUERY} -out \${BlastRefGOUT} -db $reference_genome\n";
  print RefG '      tail -10 ${BlastRefGOUT} | grep Matrix',"\n";
  print RefG "      CHECK=$?\n";
  print RefG "    done\n";
  print RefG "  else\n";
  print RefG '    tail -10 ${BlastRefGOUT} | grep Matrix',"\n";
  print RefG "    CHECK=$?\n";
  print RefG '    while [ ${CHECK} -eq 1 ]',"\n";
  print RefG "    do\n";
  print RefG "      blastn -evalue 1e-9 -show_gis -num_threads 4 -num_descriptions 2 -num_alignments 2 -query \${QUERY} -out \${BlastRefGOUT} -db $reference_genome\n";
  print RefG '      tail -10 ${BlastRefGOUT} | grep Matrix',"\n";
  print RefG "      CHECK=$?\n";
  print RefG "    done\n";
  print RefG "  fi\n";
  print RefG "fi\n";
  close RefG;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 20000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file,
    array_n  => $file_number_of_Blast_Ref_Genome,
    extra_R  => "span[hosts=1]"
  );
}

#####################################################################################
sub parse_blast_RefG {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j7_".$sample_name."_PRefG_".$$.".sh";
  open(PRefG, ">$job_files_dir/$current_job_file") or die $!;
  print PRefG "#!/bin/bash\n";
  print PRefG "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
  print PRefG "BlastRefGOUT=".$sample_name.".fa.cdhit_out.masked.goodSeq_file\${LSB_JOBINDEX}.RefGblast.out\n";
  print PRefG "BlastRefGIN=\${RefG_DIR}/".$sample_name.".fa.cdhit_out.masked.goodSeq_file\${LSB_JOBINDEX}.fa\n";
  print PRefG "PARSED=\${RefG_DIR}/".$sample_name.".fa.cdhit_out.masked.goodSeq_file\${LSB_JOBINDEX}.RefGblast.parsed\n\n";
  print PRefG 'if [ -s $BlastRefGIN ]',"\n";
  print PRefG "then\n";
  print PRefG '  if [ ! -f $PARSED ]',"\n";
  print PRefG "  then\n";
  print PRefG "    ".$run_script_path."BLASTn_RefGenome_parser.pl \${RefG_DIR} \${BlastRefGOUT} $refrence_genome_taxonomy\n";
  print PRefG '    tail -5 ${PARSED} | grep Summary',"\n";
  print PRefG "    CHECK=$?\n";
  print PRefG '    while [ ${CHECK} -eq 1 ]',"\n";
  print PRefG "    do\n";
  print PRefG "      ".$run_script_path."BLASTn_RefGenome_parser.pl \${RefG_DIR} \${BlastRefGOUT} $refrence_genome_taxonomy\n";
  print PRefG '      tail -5 ${PARSED} | grep Summary',"\n";
  print PRefG "      CHECK=$?\n";
  print PRefG "    done\n";
  print PRefG "  else\n";
  print PRefG '    tail -5 ${PARSED} | grep Summary',"\n";
  print PRefG "    CHECK=$?\n";
  print PRefG '    while [ ${CHECK} -eq 1 ]',"\n";
  print PRefG "    do\n";
  print PRefG "      ".$run_script_path."BLASTn_RefGenome_parser.pl \${RefG_DIR} \${BlastRefGOUT} $refrence_genome_taxonomy\n";
  print PRefG '      tail -5 ${PARSED} | grep Summary',"\n";
  print PRefG "      CHECK=$?\n";
  print PRefG "    done\n";
  print PRefG "  fi\n";
  print PRefG "fi\n";
  close PRefG;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 10000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file,
    array_n  => $file_number_of_Blast_Ref_Genome
  );
}

#####################################################################################
sub pool_split_for_blast_N {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j8_".$sample_name."_BN_split_".$$.".sh";
  open(BNS, ">$job_files_dir/$current_job_file") or die $!;
  print BNS "#!/bin/bash\n";
  print BNS "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
  print BNS "SAMPLE_DIR=".$sample_full_path."\n";
  print BNS "RefGFiltered_fa=".$sample_full_path."/".$sample_name.".RefGfiltered.fa\n";
  print BNS "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n\n";
  print BNS 'if [ ! -d $BN_DIR ]',"\n";
  print BNS "then\n";
  print BNS "  mkdir \${BN_DIR}\n";
  print BNS "fi\n";
  print BNS 'if [ -f $RefGFiltered_fa ]',"\n";
  print BNS "then\n";
  print BNS "  rm \${RefGFiltered_fa}\n";
  print BNS "fi\n";
  print BNS "cat \${RefG_DIR}/*.RefGfiltered.fa >> \${RefGFiltered_fa}\n";
  print BNS "".$run_script_path."check_split_BN.pl \${SAMPLE_DIR}\n";
  print BNS "CHECK=$?\n";
  print BNS 'while [ ${CHECK} -eq 10 ]',"\n";
  print BNS "do\n";
  print BNS "  ".$run_script_path."split_fasta.pl -i \${RefGFiltered_fa} -o \${BN_DIR} -n $file_number_of_Blast_N -p ".$sample_name.".RefGfiltered.fa_file\n";
  print BNS "  ".$run_script_path."check_split_BN.pl \${SAMPLE_DIR}\n";
  print BNS "  CHECK=$?\n";
  print BNS "done\n";
  close BNS;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 10000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file
  );
}

#####################################################################################
sub submit_job_array_blast_N {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j9_".$sample_name."_BN_".$$.".sh";
  open(BN, ">$job_files_dir/$current_job_file") or die $!;
  print BN "#!/bin/bash\n";
  print BN "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
  print BN "BlastNOUT=\${BN_DIR}/".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.blastn.out\n";
  print BN "QUERY=\${BN_DIR}/".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.fa\n\n";
  print BN 'if [ -s $QUERY ]',"\n";
  print BN "then\n";
  print BN '  if [ ! -f $BlastNOUT ]',"\n";
  print BN "  then\n";
  print BN "    blastn -evalue 1e-9 -show_gis -num_threads 4 -query \${QUERY} -out \${BlastNOUT} -db $db_BN\n";
  print BN '    tail -5 ${BlastNOUT} | grep Matrix',"\n";
  print BN "    CHECK1=$?\n";
  print BN '    grep "no longer exists in database" ${BlastNOUT}',"\n";
  print BN "    CHECK2=$?\n";
  print BN '    while [ ${CHECK1} -eq 1 ] || [ ${CHECK2} -eq 0 ]',"\n";
  print BN "    do\n";
  print BN "      blastn -evalue 1e-9 -show_gis -num_threads 4 -query \${QUERY} -out \${BlastNOUT} -db $db_BN\n";
  print BN '      tail -5 ${BlastNOUT} | grep Matrix',"\n";
  print BN "      CHECK1=$?\n";
  print BN '      grep "no longer exists in database" ${BlastNOUT}',"\n";
  print BN "      CHECK2=$?\n";
  print BN "    done\n";
  print BN "  else\n";
  print BN '    tail -5 ${BlastNOUT} | grep Matrix',"\n";
  print BN "    CHECK1=$?\n";
  print BN '    grep "no longer exists in database" ${BlastNOUT}',"\n";
  print BN "    CHECK2=$?\n";
  print BN '    while [ ${CHECK1} -eq 1 ] || [ ${CHECK2} -eq 0 ]',"\n";
  print BN "    do\n";
  print BN "      blastn -evalue 1e-9 -show_gis -num_threads 4 -query \${QUERY} -out \${BlastNOUT} -db $db_BN\n";
  print BN '      tail -5 ${BlastNOUT} | grep Matrix',"\n";
  print BN "      CHECK1=$?\n";
  print BN '      grep "no longer exists in database" ${BlastNOUT}',"\n";
  print BN "      CHECK2=$?\n";
  print BN "    done\n";
  print BN "  fi\n";
  print BN "fi\n";
  close BN;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 40000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file,
    array_n  => $file_number_of_Blast_N,
    extra_R  => "span[hosts=1]"
  );
}

#####################################################################################
sub parse_blast_N {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j10_".$sample_name."_PBN_".$$.".sh";
  open(PBN, ">$job_files_dir/$current_job_file") or die $!;
  print PBN "#!/bin/bash\n";
  print PBN "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
  print PBN "BlastNOUT=".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.blastn.out\n";
  print PBN "BlastNIN=\${BN_DIR}/".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.fa\n";
  print PBN "PARSED=\${BN_DIR}/".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.blastn.parsed\n\n";
  print PBN 'if [ -s $BlastNIN ]',"\n";
  print PBN "then\n";
  print PBN '  if [ ! -f $PARSED ]',"\n";
  print PBN "  then\n";
  print PBN "    ".$run_script_path."BLASTn_NT_parser.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNOUT}\n";
  print PBN "    ".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
  print PBN "    CHECK=$?\n";
  print PBN '    while [ ${CHECK} -eq 10 ]',"\n";
  print PBN "    do\n";
  print PBN "      ".$run_script_path."BLASTn_NT_parser.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNOUT}\n";
  print PBN "      ".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
  print PBN "      CHECK=$?\n";
  print PBN "    done\n";
  print PBN "  else\n";
  print PBN "    ".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
  print PBN "    CHECK=$?\n";
  print PBN '    while [ ${CHECK} -eq 10 ]',"\n";
  print PBN "    do\n";
  print PBN "      ".$run_script_path."BLASTn_NT_parser.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNOUT}\n";
  print PBN "      ".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
  print PBN "      CHECK=$?\n";
  print PBN "    done\n";
  print PBN "  fi\n";
  print PBN "fi\n";
  close PBN;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 10000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file,
    array_n  => $file_number_of_Blast_N
  );
}

#####################################################################################
sub blast_S {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j11_".$sample_name."_blastS_".$$.".sh";
  open(PS, ">$job_files_dir/$current_job_file") or die $!;
  print PS "#!/bin/bash\n";
  print PS "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
  print PS "BlastNparsed=".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.blastn.parsed\n";
  print PS "BlastNIN=\${BN_DIR}/".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.fa\n";
  print PS "OUTPUT=\${BN_DIR}/".$sample_name.".RefGfiltered.fa_file\${LSB_JOBINDEX}.blastn.summary\n\n";
  print PS 'if [ -s $BlastNIN ]',"\n";
  print PS "then\n";
  print PS '  if [ ! -f $OUTPUT ]',"\n";
  print PS "  then\n";
  print PS "    ".$run_script_path."blast_summary.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNparsed}\n";
  print PS '    grep "Finished summary" ${OUTPUT}',"\n";
  print PS "    CHECK=$?\n";
  print PS '    while [ ${CHECK} -eq 1 ]',"\n";
  print PS "    do\n";
  print PS "      ".$run_script_path."blast_summary.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNparsed}\n";
  print PS '      grep "Finished summary" ${OUTPUT}',"\n";
  print PS "      CHECK=$?\n";
  print PS "    done\n";
  print PS "  else\n";
  print PS '    grep "Finished summary" ${OUTPUT}',"\n";
  print PS "    CHECK=$?\n";
  print PS '    while [ ${CHECK} -eq 1 ]',"\n";
  print PS "    do\n";
  print PS "      ".$run_script_path."blast_summary.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNparsed}\n";
  print PS '      grep "Finished summary" ${OUTPUT}',"\n";
  print PS "      CHECK=$?\n";
  print PS "    done\n";
  print PS "  fi\n";
  print PS "fi\n";
  close PS;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 10000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file,
    array_n  => $file_number_of_Blast_N
  );
}

#####################################################################################
sub report_for_each_sample {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j12_".$sample_name."_Rep_".$$.".sh";
  open(REP, ">$job_files_dir/$current_job_file") or die $!;
  print REP "#!/bin/bash\n";
  print REP "INPUT=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq\n";
  print REP "REPORT=".$sample_full_path."/".$sample_name.".gi.AssignmentReport\n";
  print REP 'if [ -f $REPORT ]',"\n";
  print REP "then\n";
  print REP '  grep "# Finished Assignment Report" ${REPORT}',"\n";
  print REP "  CHECK=$?\n";
  print REP '  while [ ${CHECK} -eq 1 ]',"\n";
  print REP "  do\n";
  print REP "    ".$run_script_path."assignment_report_virus_gi.pl ".$sample_full_path." \${INPUT} $refrence_genome_taxonomy\n";
  print REP '    grep "# Finished Assignment Report" ${REPORT}',"\n";
  print REP "    CHECK=$?\n";
  print REP "  done\n";
  print REP "else\n";
  print REP "  ".$run_script_path."assignment_report_virus_gi.pl ".$sample_full_path." \${INPUT} $refrence_genome_taxonomy\n";
  print REP '  grep "# Finished Assignment Report" ${REPORT}',"\n";
  print REP "  CHECK=$?\n";
  print REP '  while [ ${CHECK} -eq 1 ]',"\n";
  print REP "  do\n";
  print REP "    ".$run_script_path."assignment_report_virus_gi.pl ".$sample_full_path." \${INPUT} $refrence_genome_taxonomy\n";
  print REP '    grep "# Finished Assignment Report" ${REPORT}',"\n";
  print REP "    CHECK=$?\n";
  print REP "  done\n";
  print REP "fi\n";
  close REP;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 40000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file
  );
}

#####################################################################################
sub summary_for_each_sample {
  my ($step_by_step) = @_;
  $hold_job_file = $step_by_step ? "" : $current_job_file;

  $current_job_file = "j13_".$sample_name."_Sum_".$$.".sh";
  open(SUM, ">$job_files_dir/$current_job_file") or die $!;
  print SUM "#!/bin/bash\n";
  print SUM "OUTPUT=".$sample_full_path."/".$sample_name.".gi.AssignmentSummary\n";
  print SUM "BAD_SEQ=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.badSeq\n\n";
  print SUM 'if [ -f $OUTPUT ]',"\n";
  print SUM "then\n";
  print SUM '  grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
  print SUM "  CHECK=$?\n";
  print SUM '  while [ ${CHECK} -eq 1 ]',"\n";
  print SUM "  do\n";
  print SUM "    ".$run_script_path."assignment_summary_gi.pl ".$sample_full_path." \${BAD_SEQ}\n";
  print SUM '    grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
  print SUM "    CHECK=$?\n";
  print SUM "  done\n";
  print SUM "else\n";
  print SUM "  ".$run_script_path."assignment_summary_gi.pl ".$sample_full_path." \${BAD_SEQ}\n";
  print SUM '  grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
  print SUM "  CHECK=$?\n";
  print SUM '  while [ ${CHECK} -eq 1 ]',"\n";
  print SUM "  do\n";
  print SUM "    ".$run_script_path."assignment_summary_gi.pl ".$sample_full_path." \${BAD_SEQ}\n";
  print SUM '    grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
  print SUM "    CHECK=$?\n";
  print SUM "  done\n";
  print SUM "fi\n";
  close SUM;

  bsub_docker(
    sh_file  => "$job_files_dir/$current_job_file",
    job_name => $current_job_file,
    cpus     => 1,
    mem_mb   => 40000,
    out      => "$lsf_file_dir/$current_job_file.out",
    err      => "$lsf_file_dir/$current_job_file.err",
    hold     => $hold_job_file
  );
}
