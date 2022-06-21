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
#PBS -l mem=50GB
###
### Requesting time
#PBS -l walltime=10:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

# Load all required modules for the job

# Use: IMPORTANT! for some reason, this script would only work, if I started in the data folder. All files 
# qsub -F "<sample_id> <downsample_folder_name>"" <path/to/this_shellscript.sh> 


sample_id=$1
coverage=$2

IN_1="/home/projects/cu_10148/people/nikthu/data/${coverage}/${sample_id}.bam_reads1.fastq.gz"
IN_2="/home/projects/cu_10148/people/nikthu/data/${coverage}/${sample_id}.bam_reads2.fastq.gz"
IN_3="/home/projects/cu_10148/people/nikthu/data/${coverage}/${sample_id}.bam_reads_null.fastq.gz"
IN_4="/home/projects/cu_10148/people/nikthu/data/${coverage}/${sample_id}.bam_reads_singletons.fastq.gz"
IN_5=/home/projects/cu_10148/people/nikthu/output/hisatgenotype
IN_6=$sample_id


module load tools
module load anaconda3/4.4.0
module load samtools/1.9

cd /home/projects/cu_10148/people/nikthu/tools/hisatgenotype/

export PATH=/home/projects/cu_10148/people/nikthu/tools/hisatgenotype/:/home/projects/cu_10148/people/nikthu/tools/hisatgenotype/hisat2:/home/projects/cu_10148/people/nikthu/tools/hisatgenotype/indicies:$PATH
export PYTHONPATH=/home/projects/cu_10148/people/nikthu/tools/hisatgenotype/hisatgenotype_modules:$PYTHONPATH

tempdir=$(mktemp -d)

#hisatgenotype --base hla --locus-list A,B,C,DRB1,DQB1 -p 10 --pp 2 --type-primary-exons --in-dir / --out-dir $IN_5 -1 $IN_1 -2 $IN_2
hisatgenotype --base hla --locus-list A,B,C,DRB1,DQB1 -p 10 --in-dir / --out-dir $IN_5 -1 $IN_1 -2 $IN_2


sync

# #Comparison between using all files and only paired:
# #Concat the unpaired reads:
# cat $IN_1 $IN_2 $IN_3 $IN_4 > "${tempdir}/${IN_6}_combined.fq.gz"
# #Using paired read info, but without unpaired reads, singletons etc.
# hisatgenotype --base hla --locus-list A,B,C,DRB1,DQB1 -p 10 --in-dir / --out-dir "${IN_5}/paired" -1 $IN_1 -2 $IN_2
# #Using all files as one fastqfile
# hisatgenotype --base hla --locus-list A,B,C,DRB1,DQB1 -p 10 --in-dir / --out-dir "${IN_5}/U" -U ${tempdir}/${IN_6}_combined.fq.gz 


#Bam file command
#hisatgenotype --base hla --locus-list A,B,C,DRB1,DQB1 -p 10 --bamfile /home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/NA11832.bam --in-dir / --out-dir /home/projects/cu_10148/people/nikthu/tmp

#Old command:
#hisatgenotype -x genotype_genome --base hla -z /home/projects/cu_10148/people/nikthu/tools/hisatgenotype/indicies -p 40 -1 $1 -2 $2 --in-dir / --out-dir $3

#iqsub example:
#hisatgenotype  --base hla -p 40 -1 /home/projects/cu_10148/people/nikthu/data/1000G_testdata/SRR359102_1.fastq.gz -2 /home/projects/cu_10148/people/nikthu/data/1000G_testdata/SRR359102_2.fastq.gz --in-dir / --out-dir /home/projects/cu_10148/people/nikthu/output/test
#hisatgenotype  --base hla -p 40 -1 /home/projects/cu_10148/people/nikthu/data/jens_testdata/chr6_1.fq.gz -2 /home/projects/cu_10148/people/nikthu/data/jens_testdata/chr6_2.fq.gz --in-dir / --out-dir /home/projects/cu_10148/people/nikthu/output/hisatgenotype
#hisatgenotype  --base hla -p 40 -1 /home/projects/cu_10148/people/nikthu/data/1000G_first_benchmark/NA11992.bam_reads1.fastq -2 /home/projects/cu_10148/people/nikthu/data/1000G_first_benchmark/HG00704.bam_reads2.fastq --in-dir / --out-dir /home/projects/cu_10148/people/nikthu/output/hisatgenotype
#hisatgenotype  --base hla -p 10 -1 /home/projects/cu_10148/people/nikthu/tools/hisatgenotype/ILMN/NA12892.extracted.1.fq.gz -2 /home/projects/cu_10148/people/nikthu/tools/hisatgenotype/ILMN/NA12892.extracted.2.fq.gz --in-dir / --out-dir /home/projects/cu_10148/people/nikthu/tmp


#hisatgenotype -x [GENOME] --base [GENE_GROUP] -z [INDEX_DIR] [OPTIONS] -1 [FASTQ_PAIR1] -2 [FASTQ_PAIR2]

#Log memory and time use.
qstat -f $PBS_JOBID | grep "Job Id" > "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/hisatgenotype/${IN_6}_performance.log"
qstat -f $PBS_JOBID | grep "resources_used" >> "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/hisatgenotype/${IN_6}_performance.log"


#Move the important result to the 1000 genomes database
cp "/home/projects/cu_10148/people/nikthu/output/hisatgenotype/assembly_graph-hla.${IN_6}_bam_reads1_fastq_gz-hla-extracted-"*"_fq.report" "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/hisatgenotype/${IN_6}.txt"

rm -r $tempdir
