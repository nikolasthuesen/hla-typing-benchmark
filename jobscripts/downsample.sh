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

module load samtools/1.9

IN=$1
OUT=$2
FRACTION=$3

N_CORE=40

#Downsample BAM
samtools view -@ $N_CORE -bs 42."${FRACTION}" "${IN}" > "${OUT}"


#Make fastqfiles from new bamfile
samtools collate --threads $N_CORE "${OUT}" -O | samtools fastq -@ $N_CORE -0 "${OUT}_reads_null.fastq.gz" -1 "${OUT}_reads1.fastq.gz" -2 "${OUT}_reads2.fastq.gz" -s "${OUT}_reads_singletons.fastq.gz" -