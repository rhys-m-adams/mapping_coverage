conda install -c bioconda bgv
conda install -c bioconda STAR
conda install -c bioconda samtools

#according to the manual, genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 2)
#5.623 = log2(9719)/2-1

STAR  --runMode genomeGenerate --genomeSAindexNbases 3 --runThreadN 4 --genomeDir ./index/ --genomeFastaFiles ./reference/K03455.1.fasta 

./fastq-dump --split-3 SRR961514

STAR --genomeDir ./index/ --runThreadN 4 --readFilesIn ./reads/SRR961514_1.fastq ./reads/SRR961514_2.fastq --outFileNamePrefix map1

samtools view -bq 2 file.bam > filtered.bam #quality score filter
