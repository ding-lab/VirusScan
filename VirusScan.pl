#########Song Cao###########

## a simplified version of VirusScan pipeline ##
#
#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;

#use POSIX;

my $version = "katmai_v1.1";

#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m"; 

#usage information
(my $usage = <<OUT) =~ s/\t+//g;
This script will run the virus discovery pipeline by using nohup:
 
Pipeline version: $version

$yellow Usage: perl $0  --rdir --log --step $normal

<rdir> = full path of the folder holding files for this sequence run

<step> run this pipeline step by step. (running the whole pipeline if step number is 0)

<tf> file for third-party tools, including bwa, samtools, bamtools.  

$green		[1]  Run bwa for unmapped reads

$normal
OUT

my $run_dir="";
my $log_dir="";
my $help = 0;
my $step_number = -1;

my $status = &GetOptions (
      "step=i" => \$step_number,
      "rdir=s" => \$run_dir,
      "log=s"  => \$log_dir,
      "help" => \$help,
        );

if ($help || $run_dir eq "" || $log_dir eq "" || $step_number<0) 
{
      print $usage;
      exit;
}

print "run dir=",$run_dir,"\n";
print "log dir=",$log_dir,"\n";
print "step num=",$step_number,"\n";

## read tools ##

#samtools        /diskmnt/Software/samtools-1.2/samtools
#bwa     /diskmnt/Software/bwa-0.7.17/bwa
#bamtools        /home/scao/tools/anaconda2/bin/bamtools

my $tf="./source/path_tools.tsv";

open(IN_TF,"<$tf"); 
my $samtools;
my $bwa; 
my $bamtools;

while(<IN_TF>)
{

my $ltr=$_; 
chomp($ltr);
 
my @temp=split("\t",$ltr); 
if($temp[0] eq "samtools") { $samtools=$temp[1]; }
if($temp[0] eq "bwa") { $bwa=$temp[1]; }
if($temp[0] eq "bamtools") { $bamtools=$temp[1]; }

}

close IN_TF;

# software path
#my $cd_hit = "/gscuser/mboolcha/software/cdhit/cd-hit-est";

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# path and name of databases
#my $db_BN = "/gscuser/scao/gc3027/nt/nt";
#my $db_BX = "/gscuser/scao/gc3027/nr/nr";
#my $bwa_ref = "/gscuser/scao/gc3027/fasta/virus/virusdb_082414.fa";

## virus reference ##
my $bwa_ref="./source/virusref.fa";
 
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];

# To run jobs faster, split large fasta files to small ones. Split to specific number of 
# files instead of specific sequences in each small file, because the number of job array 
# cannot be determined if spliting to specific number of sequences in each file. Job 
# number is required by qsub ${SGE_TASK_ID}. The minimum size of each file is 4kb. 
# The number of files should be determined accourding to CPUs available in the computer
# cluster.


my $HOME1=$log_dir;

if(! -d $HOME1) { `mkdir $HOME1`; }

#store job files here
if (! -d $HOME1."/tmpvirus") {
    `mkdir $HOME1"/tmpvirus"`;
}
my $job_files_dir = $HOME1."/tmpvirus";

#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR_VIRUS") {
    `mkdir $HOME1"/LSF_DIR_VIRUS"`;
}

my $lsf_file_dir = $HOME1."/LSF_DIR_VIRUS";

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


# get sample list in the run, name should not contain "."
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure
&check_input_dir($run_dir);

# start data processsing
if ($step_number<2) {
	#begin to process each sample
	for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_dir_list[$i];
		if (!($sample_name =~ /\./)) {
			$sample_full_path = $run_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
				$current_job_file="";
				if ($step_number == 0) {#run the whole pipeline
					######################################################################
					#cd-hit
					{ 
					&bsub_bwa();
					}
				######################################################################
				#run the pipeline step by step
				}elsif ($step_number == 1) {
					&bsub_bwa();
				}
			}
		}
	}
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

    `rm $lsf_out`;
    `rm $lsf_err`;

    open(BWA, ">$job_files_dir/$current_job_file") or die $!;
    print BWA "#!/bin/bash\n";
   # print BWA "#BSUB -n 1\n";
   # print BWA "#BSUB -R \"rusage[mem=20000]\"","\n";
   # print BWA "#BSUB -M 20000000\n";
   # print BWA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
   # print BWA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
   # print BWA "#BSUB -J $current_job_file\n";
    print BWA "BWA_IN=".$sample_full_path."/".$sample_name.".bam\n";
    print BWA "BWA_fq=".$sample_full_path."/".$sample_name.".fq\n";
    print BWA "BWA_sai=".$sample_full_path."/".$sample_name.".sai\n";
    print BWA "BWA_sam=".$sample_full_path."/".$sample_name.".sam\n";
    print BWA "BWA_bam=".$sample_full_path."/".$sample_name.".remapped.bam\n"; 
    #print BWA "BWA_bam=".$sample_full_path."/".$sample_name.".realign.bam\n";
    #print BWA "BWA_mapped_bam=".$sample_full_path."/".$sample_name.".mapped.bam\n";
 
    print BWA "BWA_mapped=".$sample_full_path."/".$sample_name.".mapped.reads\n";
    print BWA "BWA_fa=".$sample_full_path."/".$sample_name.".fa\n";
	#print BWA 
	print BWA 'if [ ! -s $BWA_mapped ]',"\n";
    print BWA "    then\n";
	print BWA "rm \${BWA_sai}","\n";
	print BWA "rm \${BWA_fq}","\n";
	#print BWA "mkfifo \${BWA_sai}","\n";
	print BWA "mkfifo \${BWA_fq}","\n";
	#0x100: secondary alignment
	#0x800: supplementary alignment
    #H: Hard clipping
	#S: Soft clipping
     print BWA "$samtools view -h \${BWA_IN} | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^NC/))) { print \$line;}\' | $samtools view -Sb - | $bamtools convert -format fastq > \${BWA_fq} \&","\n";

     #print BWA "$samtools view -f 4 \${BWA_IN} | $samtools view -Sb - | $bamtools convert -format fastq > \${BWA_fq} \&","\n";	
    #print BWA "bwa aln $bwa_ref -b0 \${BWA_IN} > \${BWA_sai} \&","\n";	
     print BWA "$bwa aln $bwa_ref \${BWA_fq} > \${BWA_sai}","\n";
     print BWA 'rm ${BWA_fq}',"\n";
     print BWA "mkfifo \${BWA_fq}","\n";
     print BWA "$samtools view -h \${BWA_IN} | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^NC/))) { print \$line;}\' | $samtools view -Sb - | $bamtools convert -format fastq > \${BWA_fq} \&","\n";
	#print BWA "samtools view -h \${BWA_IN} | gawk \'{if (substr(\$1,1,1)==\"\@\" || (and(\$2,0x4) || and(\$2,0x8) )) print}\' | samtools view -Sb - | bamtools convert -format fastq > \${BWA_fq} \&","\n";
    # print BWA "$samtools view -f 4 \${BWA_IN} | $samtools view -Sb - | $bamtools convert -format fastq > \${BWA_fq} \&","\n";
     print BWA "$bwa samse $bwa_ref \${BWA_sai} \${BWA_fq} > \${BWA_sam}","\n";  	
     print BWA "grep -v \@SQ \${BWA_sam} | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); if(\$ss[2]=~/^gi/) { print \$line; }\' > \${BWA_mapped}","\n";
     print BWA "$samtools view -bT $bwa_ref \${BWA_sam} > \${BWA_bam}","\n"; 
	#print BWA "     ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
    #print BWA " 	".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
     print BWA 'rm ${BWA_sam}',"\n";
     print BWA 'rm ${BWA_sai}',"\n";
	#print BWA 'rm ${BWA_fq}',"\n";
	#print BWA "else\n";
    #print BWA "     ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
    #print BWA "     ".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
    print BWA "   fi\n";
    close BWA;

    #my $sh_file=$job_files_dir."/".$current_job_file;
    #$bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -J $current_job_file -o $lsf_out -e $lsf_err sh $sh_file\n";
    #system ( $bsub_com );
 
        my $sh_file=$job_files_dir."/".$current_job_file;

        $bsub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
        print $bsub_com;
        system ($bsub_com);


   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ( $bsub_com );
}

