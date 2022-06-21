#!/bin/sh
### Account information
#PBS -W group_list=ht3_bmem -A ht3_bmem
### Mailing
#PBS -m e -M s152993@student.dtu.dk
### Resources X core, X node(s)
#PBS -l nodes=1:ppn=8:thinnode
### Requested memory (GB RAM)
#PBS -l mem=30GB
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

module load bwa/0.7.17
module load samtools/1.9

IN_FASTQ1=$1
IN_FASTQ2=$2
READ_LENGTH=$3
OUT_BAM=$4

N_CORE=8

tempdir=$(mktemp -d)

#Cut read pairs 1 and 2 - save to tmp files
/home/projects/ht3_bmem/people/nikthu/gargammel/src/fragSim --fq -l ${READ_LENGTH} ${IN_FASTQ1} > $tempdir/tmp_reads1.fq.gz

/home/projects/ht3_bmem/people/nikthu/gargammel/src/fragSim --fq -l ${READ_LENGTH} ${IN_FASTQ2} > $tempdir/tmp_reads2.fq.gz

#Map read pairs to bam file for later downsampling
bwa mem /home/projects/cu_10148/people/nikthu/data/ref_genome_data/GRCh38_full_analysis_set_plus_decoy_hla.fa -t ${N_CORE} $tempdir/tmp_reads1.fq.gz $tempdir/tmp_reads2.fq.gz | samtools view -b -o - | samtools sort --threads $N_CORE -o "${OUT_BAM}"

#remove tmp dir
rm -r $tempdir