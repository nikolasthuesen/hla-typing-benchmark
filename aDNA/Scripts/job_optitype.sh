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
#PBS -l mem=30GB
###
### Requesting time
#PBS -l walltime=06:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

# Use:
# qsub -F "<sample_id> <downsample_folder_name> " <path/to/this_shellscript.sh> 
# example:
# qsub -F "NA11832 100X" /home/projects/cu_10148/people/nikthu/scripts/tool_jobscripts/job_optitype.sh

# iqsub

sample_id=$1
read_length=$2
coverage=$3


IN_1="/home/projects/ht3_bmem/people/nikthu/shortened_reads/${read_length}_read_length/${coverage}X_coverage/${sample_id}.bam_reads1.fastq.gz"
IN_2="/home/projects/ht3_bmem/people/nikthu/shortened_reads/${read_length}_read_length/${coverage}X_coverage/${sample_id}.bam_reads2.fastq.gz"
OUTPUT_FOLDER="/home/projects/ht3_bmem/people/nikthu/output/${read_length}_read_length/${coverage}X_coverage/"


N_CORE=10

# Load all required modules for the job

#source activate optitype

module load razers3/3.4.0
module load samtools/1.9
module load hdf5/1.8.21
module load cbc/2.9.5

###Include SAMtools and ILP solver in PATH environment variable (maybe not needed, as I have used module load)
# export PATH=/services/tools/samtools/1.9/bin:$PATH
export PATH=/services/tools/cbc/2.9.5/bin/cbc:$PATH

###Add hdf5 to library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/services/tools/hdf5/1.8.21/lib



###Create new environment varible containing the path to your HDF5 installation.
export HDF5_DIR=/services/tools/hdf5/1.8.21/

###Make sure python imports are installed:
# pip install numpy==1.15.4
# pip install pyomo
# pip install pysam
# pip install matplotlib

# pip install tables
# pip install pandas
# pip install future

#For testing:
#python /home/projects/cu_10148/people/nikthu/tools/optitype/OptiTypePipeline.py -i ${tempdir}/${IN_6}_combined.fastq.gz -d -v -o $IN_5 --prefix "combinedtest"

python /home/projects/cu_10148/people/nikthu/tools/optitype/OptiTypePipeline.py -i $IN_1 $IN_2 -d -v -o $OUTPUT_FOLDER --prefix $sample_id

rm "${OUTPUT_FOLDER}/${sample_id}_coverage_plot.pdf"

#Log memory and time use.
qstat -f $PBS_JOBID | grep "Job Id" > "${OUTPUT_FOLDER}/${sample_id}_performance.log"
qstat -f $PBS_JOBID | grep "resources_used" >> "${OUTPUT_FOLDER}/${sample_id}_performance.log"
