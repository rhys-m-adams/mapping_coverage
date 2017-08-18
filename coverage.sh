mkdir $3_index
STAR  --runMode genomeGenerate --genomeSAindexNbases 3 --runThreadN 4 --genomeDir $3_index --genomeFastaFiles $3
STAR --genomeDir $3_index runThreadN 4 --readFilesIn $1 $2 --outFileNamePrefix $4

samtools view -bS $4Aligned.out.sam > $4.bam
#samtools view -bq 2 file.bam > filtered.bam #quality score filter
samtools sort $4.bam -o $4.sort.bam
samtools depth $4.sort.bam > $4.coverage
samtools faidx $3
samtools mpileup -g -f $3 $4.sort.bam >$4.bcf
bcftools view -v snps $4.bcf > $4.vcf
samtools index $4.sort.bam

#printf "genome $3 \n load $4.sort.bam \n load $4.vcf \n sort \n collapse \n snapshot $4.png \n exit \n" > igv_script.sh
#igv -b igv_script.sh
