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
#PBS -l mem=20GB
###
### Requesting time
#PBS -l walltime=01:00:00
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e /home/projects/cu_10148/people/nikthu/bin
#PBS -o /home/projects/cu_10148/people/nikthu/bin

#qsub -F " <path/to/input_1.cram>" <path/to/this_shellscript.sh>
#qsub -F "<prefix>" /home/projects/cu_10148/people/nikthu/scripts/measure_depth.sh
#iqsub example
# IN=/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/cram_downloads/NA19818.cram

module load htslib/1.9
module load mosdepth/0.2.6

IN="/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/cram_downloads/${1}.cram"

mosdepth --by /home/projects/cu_10148/people/nikthu/data/ref_genome_data/output_1000G_Exome.v1.bed -f /home/projects/cu_10148/people/nikthu/data/ref_genome_data/GRCh38_full_analysis_set_plus_decoy_hla.fa --threads 10 --no-per-base --thresholds 0,1,2,5 "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/depth/${1}.depth" "${IN}"

#keep HG00380.depth.mosdepth.region.dist.txt and HG00380.depth.mosdepth.global.dist.txt and HG00380.depth.mosdepth.summary.txt
rm "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/depth/${1}.depth.regions.bed.gz" "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/depth/${1}.depth.regions.bed.gz.csi" "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/depth/${1}.depth.thresholds.bed.gz" "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/depth/${1}.depth.thresholds.bed.gz.csi"