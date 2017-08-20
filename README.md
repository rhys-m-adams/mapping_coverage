# mapping_coverage
#This repository details a basic pipeline for mapping RNA-seq data to a genomic reference. 
#  
#The following files were created:  
#analyze_coverage.py - outputs figures based on coverage and GC scores  
#correlation.pdf - sample output for test HIV data sets, with QC score threshold 38  
#coverage.pdf - sample output for test HIV data sets comparing GC to coverage, with QC filter of 38  
#install.md - files needed to be installed for coverage to work  
#mapping_report.bib, mapping_report.pdf, mapping_report.tex - tex file analyzing results   
#STAR_out - folder containg mapping figures with quality score thresholds 0, 10, 20, 30  
#coverage - main pipeline.  Usage is  
#coverage -1 forward_reads.fastq -2 reverse_reads.fastq -x reference_genome.fasta -Q 30 [optional: ] -o out_name -N 3 -t 4  
#-1  forward fastq reads  
#-2  referse fastq reads  
#-x  reference genome  
#-Q  reads with average quality scores less than this option are filtered out (e.g. less than 30)  
#-o  output name prefix  
#-N  parameter used by STAR. According to the manual, it(i.e. genomeSAindexNbases) must be scaled down to min(14, log2(GenomeLength)/2 - 2)  
#-t  number of threads used by STAR to map reads to genome.  
#  
#To complete the data set I downloaded the following files. 
#Using the sratoolkit from ncbi, I downloaded the files for the MiSeq data set:  
#http://www.ncbi.nlm.nih.gov/sra/?term=SRR961514  
#by running the command  
./fastq-dump --split-3 SRR961514  
  
#I downloaded the reference HIV genome in fasta format from:  
#http://www.ncbi.nlm.nih.gov/nuccore/K03455.1  
#The Miseq data sets were put into the ./reads directory, and the reference was placed into the ./reference directory  
#  
#To test my pipeline, I ran it filtering quality scores of 0, 10, 20, and 30.  
#I tested my pipeline with the following commands  

for QC in {0,10,20,30}  
do  
./coverage -1 ./reads/SRR961514_1.fastq -2 ./reads/SRR961514_2.fastq  -x ./reference/K03455.1.fasta -Q $QC -o ./STAR_out/pipeline_test_QC$QC  
done  

#The major results are output to   
#./STAR_out/pipline_test_QC<QC filter>coverage.pdf  
#./STAR_out/pipline_test_QC<QC filter>coverage.tsv  
#./STAR_out/pipline_test_QC<QC filter>correlation.pdf  

#I then tested if no output name yields coverage.pdf, coverage.tsv, and correlation.pdf files  
./coverage -1 ./reads/SRR961514_1.fastq -2 ./reads/SRR961514_2.fastq  -x ./reference/K03455.1.fasta -Q 38  

#The major results were output to   
#./coverage.pdf  
#./coverage.tsv  
#./correlation.pdf  
#satisfying the requirements of the pipeline.  
