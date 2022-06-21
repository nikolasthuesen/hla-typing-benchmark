#!/bin/sh
###Note: No commands may be executed until after the #PBS lines
### Notice how the options are used on the #PBS lines
### Account information
#PBS -W group_list=ht3_bmem -A ht3_bmem
###
### Mailing
#PBS -m e -M s152993@student.dtu.dk
###
### Resources X core, X node(s)
#PBS -l nodes=1:ppn=40:thinnode
###
### Requested memory (GB RAM)
#PBS -l mem=60GB
###
### Requesting time
#PBS -l walltime=06:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

#qsub -F "<input_1.fq.gz> <path/to/output_sample_name>" <path/to/this_shellscript.sh>
#qsub -F "/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/HG01357.cram" /home/projects/cu_10148/people/nikthu/scripts/index_cram.sh
#iqsub example
# IN=/home/projects/cu_10148/people/nikthu/data/1000G_first_benchmark/NA11992.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram
# OUT=/home/projects/cu_10148/people/nikthu/tmp/NA11992.bam

module load samtools/1.9

IN=$1
N_CORE=40

samtools index -@N_CORE "${IN}"