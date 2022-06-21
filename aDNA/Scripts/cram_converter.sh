#!/bin/sh
### Account information
#PBS -W group_list=ht3_bmem -A ht3_bmem
### Mailing
#PBS -m e -M s152993@student.dtu.dk
### Resources X core, X node(s)
#PBS -l nodes=1:ppn=40:thinnode
### Requested memory (GB RAM)
#PBS -l mem=120GB
### Requesting time
#PBS -l walltime=06:00:00
### Output files
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

#qsub -F "<input_1.fq.gz> <path/to/output_sample_name>" <path/to/this_shellscript.sh>
#qsub -F "/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/cram_downloads/HG01341.cram /home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/HG01341.bam  " /home/projects/cu_10148/people/nikthu/scripts/cram_converter.sh
#iqsub example
# IN=/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/cram_downloads/HG00096.cram
# OUT=/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/HG00096.bam

cd /home/projects/cu_10148/people/nikthu/tmp

module load samtools/1.9
module load bedtools/2.28.0

IN=$1
OUT=$2
N_CORE=40

tempdir=$(mktemp -d)

#chr6 coordinates taken from the HLA-LA article (widest region, I could find to include as much as needed)
#reference genome taken from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

#Extract all aligned reads from cram file, sort and fix mates.
samtools view -b -@ $N_CORE -T /home/projects/cu_10148/people/nikthu/data/ref_genome_data/GRCh38_full_analysis_set_plus_decoy_hla.fa -o - "${IN}" | samtools sort --threads $N_CORE -n - | samtools fixmate --threads $N_CORE -m - $tempdir/fullbam_sorted.bam

#Extract all reads mapping to the HLA region in chr6 as well as noted HLA contigs and their mates.
bedtools pairtobed -type either -abam $tempdir/fullbam_sorted.bam -b /home/projects/cu_10148/people/nikthu/data/ref_genome_data/hla_headers.bed > $tempdir/tmp1.bam

#Extract unaligned reads, where mate is also unaligned
samtools view -b -u -o $tempdir/tmp2.bam -f 12 -T /home/projects/cu_10148/people/nikthu/data/ref_genome_data/GRCh38_full_analysis_set_plus_decoy_hla.fa -@ $N_CORE "${IN}" 

#Merge all bam files
samtools merge -@ $N_CORE -u - $tempdir/tmp[12].bam | samtools sort --threads $N_CORE -o "${OUT}"

sync

#Convert it all to paired fastqfiles (-1, -2), singleton reads (-s) as well as reads, where READ1 and READ2 FLAG bits set are either both set or both unset (-0)
samtools collate --threads $N_CORE "${OUT}" -O | samtools fastq -@ $N_CORE -0 "${OUT}_reads_null.fastq.gz" -1 "${OUT}_reads1.fastq.gz" -2 "${OUT}_reads2.fastq.gz" -s "${OUT}_reads_singletons.fastq.gz" -

rm -r $tempdir