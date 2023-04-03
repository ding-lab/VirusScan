
#!/usr/bin/perl
use strict;

my $usage = "
This script will read corresponding files in the given director and 
generate a report which contains SampleDescription, SequenceReport,
AssignmentSummary, InterestingReads.

perl $0 <run folder> <program version>
<run folder> = full path of the folder holding files for this sequence run

";

die $usage unless scalar @ARGV == 2;

my ( $dir, $version ) = @ARGV;

my @temp = split("/", $dir);
my $run_name = pop @temp;
my $outFile = $dir."/Analysis_Report_".$run_name.".txt";

my $virus_ref="./source/pathogenpdx.031823.fa"; 

open(IN,"<$virus_ref"); 

my %gi2name=();

while(<IN>)
{
my $l=$_; 
chomp($l); 
if($l=~/^>gi/) 
{
@temp=split(" ",$l); 
my $gi=$temp[0];
$gi=~s/\>//g; 
my $name=$temp[1]; 
#print $gi,"\n";
for(my $i=2;$i<@temp;$i++)
{
$name.=" ".$temp[$i];
}
$gi2name{$gi}=$name; 
}
}

open (OUT, ">$outFile") or die "can not open file $outFile!\n";

my ($wkday,$month,$day,$time,$year) = split(/\s+/, localtime);
print OUT "VirusScan ${version}; Processing date: $day-$month-$year\n";

my $c = "**************************************************************************\n";
my $c2 = "#########################################################################\n";
print OUT $c;

print OUT "Summary:\n\n";

print OUT "Sample","\t","totalreads","\t","Virus/Bacteria","\t","supportingreads","\t","RPHM","\n"; 

&generate_AssignmentSummary( $dir );

print OUT "\nEnd of Summary\n";
#print OUT $c ;

#print OUT "\n\nSequence Report\n\n";
#&generate_SequenceReport( $dir );
#print OUT "End of Sequence Report\n\n";
#print OUT $c ;

#print OUT "\n\nTaxonomy Assignment:\n\n";
#&generate_AssignmentSummary( $dir );
#print OUT "End of Assignment\n\n";
#print OUT $c ;

#print OUT "\n\nInteresting Reads\n\n";
#&generate_InterestingReads( $dir );
#print OUT "End of Interesting Reads\n\n";
#print "\n";

#print OUT "# Finished\n";

exit; 
#####################################################################
# Assignment Summary
sub generate_AssignmentSummary {
	my ( $dir ) = @_;
	my $n; 
	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $full_path = $dir."/".$name;
#		print $full_path,"\n";
		if (!($name =~ /\./)) {
			if (-d $full_path) { # is a directory
				my $Summary_file = $full_path."/".$name.".mapped.reads";
				my $readcount_file = $full_path."/".$name.".reads.count.tsv";
				if (-e $Summary_file) {
					my %vreads=();
					my $totreads;
					open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
					open (INC,$readcount_file) or die "can not open file $Summary_file!\n";

					while (<IN>) {
					my $l=$_; 
					chomp($l); 
					@temp=split("\t",$l); 
#					print $temp[2],"\n"; 
					#<STDIN>;
					if($temp[4]>0)
					{
					$n=$gi2name{$temp[2]}; 
#					print $n,"\n";
#					<STDIN>;
					$vreads{$n}{$temp[0]}++;							
					}
					}

					while (<INC>)
					{
					my $l=$_;
                                        chomp($l);
					$totreads=$l; 
					}

				if(keys %vreads)
				{	
				foreach $n (sort keys %vreads) 
				{
				my $count=0; 
				foreach $c (sort keys %{$vreads{$n}})
				{
					$count++; 
				}
				print OUT $name, "\t",$totreads, "\t", $n, "\t", $count,"\t",sprintf("%.4f",$count/$totreads*(100000000)),"\n"; } 
				} 	
				else { print OUT $name, "\t",$totreads,"\t","Virus/Bacteria","\t","0","\t","0","\n"; }
			
	#			print OUT $c2 ;
				close IN; 
			}
			}
		}
	}
}

