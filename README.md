# mapping_coverage

#Using the sratoolkit from ncbi, I downloaded the files for the MiSeq data set:  
#http://www.ncbi.nlm.nih.gov/sra/?term=SRR961514  
#and ran the command  
#./fastq-dump --split-3 SRR961514  
#I downloaded the reference HIV genome in fasta format:  
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
