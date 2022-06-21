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
#PBS -l mem=120GB
###
### Requesting time
#PBS -l walltime=01:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

#qsub -F "/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/HG01357.bam" /home/projects/cu_10148/people/nikthu/scripts/index_bam.sh


module load samtools/1.9

IN=$1
N_CORE=40

samtools index -@N_CORE "${IN}"