#### VirusScan version 1.1 ####

Author: Song Cao

Released on Apr 25, 2016

VirusScan pipeline is a fully automated and modular software package designed for the fast 
and accurate detection of known viruses from NGS data [1]. It works on LSF job scheduler. 

It was developed from VirusHunter pipeline, which focuses on identification of novel viruses for 454 reads. 
Compared to VirusHunter pipeline, VirusScan can work on Illlumina WGS, WES and RNA-Seq data and can fastly return 
the discovery result of known viruses.  

Dependencies:

1. RepeatMasker: Download and install RepeatMasker from http://www.repeatmasker.org/RMDownload.html

2. BLAST Module: Download and install BLAST from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

3. MySQL DBI: See http://search.cpan.org/dist/DBI/. DBI may be included with your Linux distribution by default. This is a Perl module that allows Perl to interact directly with MySQL database.

4. BioPerl: See http://bioperl.org/. BioPerl is used for parsing BLAST output files, and to construct taxonomy lineage tree from a taxonomy ID.  

5. NCBI nt database: Download NT database from ftp://ftp.ncbi.nih.gov/blast/db/

6. NCBI taxonomy database: Download NCBI taxonomy database from ftp://ftp.ncbi.nih.gov/pub/taxonomy/. 


6.1. Create a directory to hold taxonomy file, e.g. taxdump_2016_06_20
Download taxdump.tar.gz file to the directory and Try "unzip taxdump.tar.gz" to unzip the file. 

6.2. Create MySQL Database for the taxonomy information

Ask your MySQL database administrator to create a MySQL database for taxonomy information, and grant privileges on this database to a suitable username. 
For example, ask your MySQL database administrator to use following commands to create a database named "vs_taxondb" and grant all privileges to the user "vs_taxonUser" with the password "vs_password".
$ mysql -u root -p
CREATE DATABASE test_taxondb;
GRANT ALL ON vs_taxondb.* TO 'vs_taxonUser'@'localhost' IDENTIFIED BY 'vs_password';
GRANT ALL ON vs_taxondb.* TO 'vs_taxonUser'@'%' IDENTIFIED BY 'vs_password';
QUIT;

6.3. Load gi-taxid into database for nucleotide sequences:
download gi_taxid_nucl.dmp.gz to the directory
unzip the file

Modify script "import_gi_taxid_nucl.sql " to replace the full path to the gi_taxid_nucl.dmp file with the actual full path in your local system at line " LOAD DATA LOCAL INFILE" in the script. 
The LOAD DATA INFILE statement reads rows from a text file into a table at a very high speed. The file name must be given as a literal string.

Load the gi_taxid_nucl.dmp content to a MySQL database using script " import_gi_taxid_nucl.sql" with the following command:
Cat import_gi_taxid_nucl.sql | mysql -h hostname --user=username databaseName --pass=password &

Warning: This can take a very long time. It is better to run it as a background task. 

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

[2] Split files for running RepeatMasker

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

Please cite the following paper for VirusScan pipeline:

1. Song Cao, Michael C. Wendl, Matthew A. Wyczalkowski, Kristine Wylie, Kai Ye, Reyka Jayasinghe, Mingchao Xie, Song Wu, Beifang Niu, Robert Grubb III, Kimberly J. Johnson, Hiram Gay, Ken Chen, Janet S. Rader,  John F. Dipersio, Feng Chen, and Li Ding, Divergent viral presentation among human tumors and adjacent normal tissues (2016), 6:28294. 

