#### VirusScan version 1.1 ####

Author: Song Cao

Contact: scao@wustl.edu

Released on Apr 3 2019

Please cite the following paper for VirusScan pipeline:

Song Cao, Michael C. Wendl, Matthew A. Wyczalkowski, Kristine Wylie, Kai Ye, Reyka Jayasinghe, Mingchao Xie, Song Wu, Beifang Niu, Robert Grubb III, Kimberly J. Johnson, Hiram Gay, Ken Chen, Janet S. Rader,  John F. Dipersio, Feng Chen, and Li Ding, Divergent viral presentation among human tumors and adjacent normal tissues, Scientific Reports, 2016, 6:28294. 

A simplified version of VirusScan, no dependency on taxonomy database. A customized virus database was created by using viruses found in above TCGA panvirus paper.       


###Usage:###

git clone https://github.com/ding-lab/VirusScan.git
git checkout simplified

change ./source/path_tools.tsv for the software directory (samtools, bwa, bamtools) in your system. 

##Install third-party tools (samtools, bwa, bamtools)##

1. Install anaconda  

wget https://repo.anaconda.com/archive/Anaconda2-5.3.1-Linux-x86_64.sh

2. Install samtools, bwa, bamtools

conda install -c bioconda samtools

conda install -c bioconda bwa

conda install -c bioconda bamtools

## submit jobs in the directory that contains VirusScan scripts ##

Usage: perl VirusScan.pl  --rdir --log --step 

<rdir> = full path of the folder holding files for this sequence run

<log> = directory for log files 

<step> run this pipeline step by step. (running the whole pipeline if step number is 0)


[1]  Run bwa for unmapped reads againt virus reference

[2]  Generate summary for virus discovery


rdir: A folder contains different bam files for different samples: 

For example: 

work/sample1/sample1.bam 

work/sample2/sample2.bam

Warning: The prefix of the name of the bam file should be the same as the sample directory.



