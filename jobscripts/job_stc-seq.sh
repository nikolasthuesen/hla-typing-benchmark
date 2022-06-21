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
#PBS -l mem=60GB
###
### Requesting time
#PBS -l walltime=06:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

#Use: IMPORTANT: Input should be unzipped and should be provided without .fq ending
#qsub -F "<prefix_name>" <paht/to/this_shellscript.sh>
#qsub -F "HG00339" /home/projects/cu_10148/people/nikthu/scripts/tool_jobscripts/job_stc-seq_single_core.sh

#bash
#bash /home/projects/cu_10148/people/nikthu/scripts/tool_jobscripts/job_stc-seq.sh /home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/HG00380.bam_reads1.fastq /home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/HG00380.bam_reads2.fastq HG00380

# Load all required modules for the job

cd /home/projects/cu_10148/people/nikthu/output/stc-seq

module load bowtie/1.2.1.1
module load bowtie2/2.3.4.1

module load jdk/14
module load bbmap/38.35

module load gcc
module load intel/perflibs/2019
module load R/3.6.1

#bbmap on computerome needs the .sh addition
alias bbmap="bbmap.sh"

#Add bin directory to path
export PATH=$PATH:/home/projects/cu_10148/people/nikthu/tools/stc-seq/bin

###First time running this program - generating all possible artificial reads (70bp)
###/home/projects/cu_10148/people/nikthu/tools/stc-seq/bin/getArtifical_reads /home/projects/cu_10148/people/nikthu/tools/stc-seq/data/hla.exon.fast70bp.fastq

#make sure, that all the relevant files also exist in the stc-seq output folder (outcommented, so as to not try to copy every run)
#cp -r /home/projects/cu_10148/people/nikthu/tools/stc-seq/data/* ./

if [[ ! -f /home/projects/cu_10148/people/nikthu/output/stc-seq/cigar.txt ]] ; then
    cp -r /home/projects/cu_10148/people/nikthu/tools/stc-seq/data/* ./
fi


#For Iqsub:
#sample_id=NA18504

sample_id=$1

IN_1="/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/${sample_id}.bam_reads1.fastq.gz"
IN_2="/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/${sample_id}.bam_reads2.fastq.gz"
IN_3="/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/${sample_id}.bam_reads_null.fastq.gz"
IN_4="/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/${sample_id}.bam_reads_singletons.fastq.gz"
IN_5=/home/projects/cu_10148/people/nikthu/output/stc-seq
IN_6="${sample_id}"

mkdir "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}"


#cat $IN_1 $IN_2 | gunzip > "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}_combined.fq"

cat $IN_1 $IN_2 $IN_3 $IN_4 | gunzip > "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined.fq"

IN="/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined"



/home/projects/cu_10148/people/nikthu/tools/stc-seq/script/step_1.sh "${IN}"

/home/projects/cu_10148/people/nikthu/tools/stc-seq/script/step_2.sh "${IN}"

#Added Rscript in front, as R isn't executed from /usr/bin/R
Rscript /home/projects/cu_10148/people/nikthu/tools/stc-seq/bin/Initial_screening.R "${IN}_position-array.txt" 7 15 "${IN}_second_retain_allele.txt"    

/home/projects/cu_10148/people/nikthu/tools/stc-seq/script/step_3.sh "${IN}"

#Move the important result to the 1000 genomes database
cp "${IN}_report_null.txt" "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/stc-seq/${IN_6}.txt"

sync
#rm "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined.fq" "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined_first_retain_allele."* "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined_second_retain_allele."* "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined_real_allele_p"*  "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined_de_duplication."* "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}/${IN_6}_combined_false_allele_pool_"*

#mv "${IN}_"* /home/projects/cu_10148/people/nikthu/output/stc-seq

#rm cigar.txt  connect_exon_70bp.fasta  exon-70bp.fastq  hla.exon.fasta  nullAllele_list.txt


#Log memory and time use.
qstat -f $PBS_JOBID | grep "Job Id" > "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/stc-seq/${IN_6}_performance.log"
qstat -f $PBS_JOBID | grep "resources_used" >> "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/stc-seq/${IN_6}_performance.log"


rm -r "/home/projects/cu_10148/people/nikthu/output/stc-seq/${IN_6}"


