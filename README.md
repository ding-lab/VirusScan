#### VirusScan ####
### released on Apr 25, 2016####

VirusScan pipeline is a fully automated and modular software package designed for the fast 
and accurate detection of known viruses from NGS data [1]. It works on LSF job scheduler. 

It was developed from VirusHunter pipeline [2], which focuses on identification of novel viruses for 454 reads. 
Compared to VirusHunter pipeline, VirusScan can work on Illlumina WGS, WES and RNA-Seq data and can return 
the discovery of known viruses from sequencing data very fast.  

Dependencies:

1. RepeatMasker: Download and install RepeatMasker from http://www.repeatmasker.org/RMDownload.html

2. BLAST Module: Download and install BLAST from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

3. MySQL DBI: See http://search.cpan.org/dist/DBI/. DBI may be included with your Linux distribution by default. This is a Perl module that allows Perl to interact directly with MySQL database.

4. BioPerl: See http://www.bioperl.org/wiki/Getting_BioPerl. BioPerl is used for parsing BLAST output files, and to construct taxonomy lineage tree from a taxonomy ID.  

5. NCBI nt database: Download NT database from ftp://ftp.ncbi.nih.gov/blast/db/

6. NCBI taxonomy database: Download NCBI taxonomy database from ftp://ftp.ncbi.nih.gov/pub/taxonomy/; See http://pathology.wustl.edu/VirusHunter/ for how to create taxonomy database. 

Usage: 

git clone https://github.com/ding-lab/VirusScan.git

perl VirusScan.pl < run_folder > < step_number >

run_folder: A folder contains different bam files for different samples: 
For example: 
work/sample1/sample1.bam 
work/sample2/sample2.bam

Warning: The prefix of the name of the bam file should be the same as the sample directory.

step_number: Integer between 1 and 33 which represents the following step: 

[1] Extract unmapped no-human reads from aligned bam file and map extracted reads to the viral database

[2] Split files for RepeatMasker

[3] Submit RepeatMasker job array

[4] Sequence Quality Control

[5] Split files for Blast Human Genome

[6] Submit Blast Human Genome job array

[7] Parse Human Genome Blast result

[8] Pool and split files for BlastN

[9] Submit BlastN job array

[10] Parse BlastN result

[11] Generate summary result for blastn output

[12] Assignment report for each sample

[13] Assignment summary for each sample

[14] Generate report for the run

[22] Run steps from 2 to 14

[23] Run steps from 3 to 14

[24] Run steps from 4 to 14

[25] Run steps from 5 to 14

[26] Run steps from 6 to 14

[27] Run steps from 7 to 14

[28] Run steps from 8 to 14

[29] Run steps from 9 to 14

[30] Run steps from 10 to 14

[31] Run steps from 11 to 14

[32] Run steps from 12 to 14

[33] Run steps from 13 to 14 

1. Song Cao, Michael C. Wendl, Matthew A. Wyczalkowski, Kristine Wylie, Kai Ye, Reyka Jayasinghe, Mingchao Xie, Song Wu, Beifang Niu, Robert Grubb III, Kimberly J. Johnson, Hiram Gay, Ken Chen, Janet S. Rader,  John F. Dipersio, Feng Chen, and Li Ding, Divergent viral presentation among human tumors and adjacent normal tissues, submitted. 

2. Identification of Novel Viruses Using VirusHunter -- An Automated Data Analysis Pipeline. Guoyan Zhao, Siddharth Krishnamurthy, Zhengqiu Cai, Vsevolod L. 
Popov, Amelia P. Travassos da Rosa, Hilda Guzman, Song Cao, Herbert W. Virgin, Robert B. Tesh and David Wang, PLoS One. 2013. 8(10):e78470. 
