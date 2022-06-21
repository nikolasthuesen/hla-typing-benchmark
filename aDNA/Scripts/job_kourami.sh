#!/bin/sh
### Account information
#PBS -W group_list=ht3_bmem -A ht3_bmem
###
### Mailing
#PBS -m e -M s152993@student.dtu.dk
###
### Resources X core, X node(s)
### allignandextract script has 8 cores as standard. I just used those
#PBS -l nodes=1:ppn=10:thinnode
###
### Requested memory (GB RAM)
#PBS -l mem=50GB
###
### Requesting time
#PBS -l walltime=03:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

#Use:
#qsub -F "<distinguishable_sample_name.bam> </path/to/input> </path/to/outputfolder>" <paht/to/this_shellscript.sh> 
#example: 
#qsub -F "/home/projects/cu_10148/people/nikthu/data/5X/NA20753.bam NA20753 " /home/projects/cu_10148/people/nikthu/scripts/tool_jobscripts/job_kourami.sh


IN_1=$1
IN_2=$2

source activate kourami

# Load all required modules for the job
module load tools
# module load anaconda3/4.4.0
module load java/1.8.0
module load apache-maven/3.6.0 
module load bwa/0.7.15
module load bamtools/2.5.1 
module load samtools/1.9

export PATH=/home/projects/cu_10148/people/nikthu/tools/kourami/bamUtil:$PATH

#Convert fastq files to BAM files by aligning to the reference given by Kourami:
#Read, extract and realign to Kourami:

cd /home/projects/cu_10148/people/nikthu/tools/kourami/kourami-0.9.6/bam_aligned_to_kourami

/home/projects/cu_10148/people/nikthu/tools/kourami/kourami-0.9.6/scripts/alignAndExtract_hs38DH.sh -d /home/projects/cu_10148/people/nikthu/tools/kourami/kourami-0.9.6/db -r /home/projects/cu_10148/people/nikthu/tools/kourami/kourami-0.9.6/resources/hs38NoAltDH.fa $IN_2 $IN_1

#sleep 10s

#Perform analysis

cd /home/projects/cu_10148/people/nikthu/output/kourami

java -jar /home/projects/cu_10148/people/nikthu/tools/kourami/kourami-0.9.6/target/Kourami.jar -d /home/projects/cu_10148/people/nikthu/tools/kourami/kourami-0.9.6/db -o $IN_2 "/home/projects/cu_10148/people/nikthu/tools/kourami/kourami-0.9.6/bam_aligned_to_kourami/${IN_2}_on_KouramiPanel.bam"


#Move the important result to the 1000 genomes database
cp "/home/projects/cu_10148/people/nikthu/output/kourami/${IN_2}.result" "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/kourami/${IN_2}.txt"

#Log memory and time use.
qstat -f $PBS_JOBID | grep "Job Id" > "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/kourami/${IN_2}_performance.log"

qstat -f $PBS_JOBID | grep "resources_used" >> "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/kourami/${IN_2}_performance.log"

