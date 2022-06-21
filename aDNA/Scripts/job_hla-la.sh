#!/bin/sh
### Account information
#PBS -W group_list=ht3_bmem -A ht3_bmem
###
### Mailing
#PBS -m e -M s152993@student.dtu.dk
###
### Resources X core, X node(s)
#PBS -l nodes=1:ppn=10:thinnode
###
### Requested memory (GB RAM)
#PBS -l mem=180GB
###
### Requesting time
#PBS -l walltime=48:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

# Load all required modules for the job

#Use:
#qsub -F "<path/to/input.bam> <path/to/outputfolder> <specific_output_directory_name>" <path/to/this_shellscript.sh> 
#Example
#qsub -F "/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/HG00365.bam /home/projects/cu_10148/people/nikthu/output/hla-la HG01049" /home/projects/cu_10148/people/nikthu/scripts/tool_jobscripts/job_hla-la.sh
#qsub -F "/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/cram_downloads/HG00365.cram /home/projects/cu_10148/people/nikthu/output/hla-la HG00365" /home/projects/cu_10148/people/nikthu/scripts/tool_jobscripts/job_hla-la.sh
#qsub -F "/home/projects/cu_10148/people/nikthu/tmp2/NA18504/NA18504_end_sorted.bam /home/projects/cu_10148/people/nikthu/tmp2 NA18504_HLA_LA_test" /home/projects/cu_10148/people/nikthu/scripts/tool_jobscripts/job_hla-la.sh


#For iqsub:
# IN1=/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/NA19909.cram
# IN2=/home/projects/cu_10148/people/nikthu/output/hla-la
# IN3=NA19909

module load tools
module load anaconda3/4.4.0

source activate hla-la

module load samtools/1.9

cd /home/projects/cu_10148/people/nikthu/conda/envs/hla-la/opt/hla-la/working

N_CORE=10

echo $1 $2 $3

HLA-LA.pl --BAM $1 --graph PRG_MHC_GRCh38_withIMGT --sampleID $3 --maxThreads $N_CORE --workingDir $2

#HLA-LA.pl --BAM ${IN1} --graph PRG_MHC_GRCh38_withIMGT --sampleID ${IN3} --samtools_T /home/projects/cu_10148/people/nikthu/data/ref_genome_data/GRCh38_full_analysis_set_plus_decoy_hla.fa --maxThreads 10 --workingDir ${IN2}


#HLA-LA.pl --BAM $1 --graph PRG_MHC_GRCh38_withIMGT --sampleID $3 --maxThreads 40 --workingDir $2

#iqsub example:
#HLA-LA.pl --BAM /home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/cram_downloads/NA18504.cram --graph PRG_MHC_GRCh38_withIMGT --sampleID NA18504 --samtools_T /home/projects/cu_10148/people/nikthu/data/ref_genome_data/GRCh38_full_analysis_set_plus_decoy_hla.fa --maxThreads 10 --workingDir /home/projects/cu_10148/people/nikthu/tmp2

#HLA-LA.pl --BAM /home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/NA19909.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID NA19909 --maxThreads $N_CORE --workingDir /home/projects/cu_10148/people/nikthu/tmp

#Log memory and time use.

qstat -f $PBS_JOBID | grep "Job Id" > "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/stc-seq/hla-la/${3}_performance.log"
qstat -f $PBS_JOBID | grep "resources_used" >> "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/hla-la/${3}_performance.log"


#Move the important result to the 1000 genomes database
cp "/home/projects/cu_10148/people/nikthu/output/hla-la/${3}/hla/R1_bestguess_G.txt" "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/hla-la/${3}.txt"

sync

rm -r "/home/projects/cu_10148/people/nikthu/output/hla-la/${3}"
