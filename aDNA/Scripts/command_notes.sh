mkdir 35_read_length 45_read_length 55_read_length 65_read_length 
mkdir 1X_coverage 2X_coverage 5X_coverage 10X_coverage 20X_coverage 50X_coverage full_coverage test
mkdir 1X_coverage 2X_coverage 5X_coverage 10X_coverage 20X_coverage 50X_coverage full_coverage


while read sample_id; do
  mv /home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/${sample_id}.bam_reads*.fastq.gz /home/projects/ht3_bmem/people/nikthu/data/
done </home/projects/ht3_bmem/people/nikthu/scripts/config.yaml



samtools view -h mapped.bam | grep -e '^@' -e 'readName' | samtools stats | grep '^SN' | cut -f 2-


SN      average first fragment length: 


samtools stats HG00096.bam | grep "average first fragment length:" | awk '{print $NF}'



C=$(echo "$A / $sum" | bc -l)


$original_coverage * $cut_bam_read_length / $original_average_read_length 

echo "$a $b" | awk '{print $1 - $2}'