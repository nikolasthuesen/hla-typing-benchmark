coverage="100X"

configfile: "config_" + coverage + ".yaml"

def mem(gb=0,mb=0,kb=0,b=0):
    return gb*10**9 + mb*10**6 + kb*10**3 + b

data_dir = "/home/projects/cu_10148/people/nikthu/data/" + coverage + "/"
script_dir = "/home/projects/cu_10148/people/nikthu/scripts/"
output_dir = "/home/projects/cu_10148/people/nikthu/output/1000_genomes_results/"

rule all:
    input:
        expand(output_dir + "kourami/{sample_id}.txt", sample_id = config["sample_downsample_fraction"]),
        expand(output_dir + "hisatgenotype/{sample_id}.txt", sample_id = config["sample_downsample_fraction"]),
        expand(output_dir + "stc-seq/{sample_id}.txt", sample_id = config["sample_downsample_fraction"]),
        expand(output_dir + "hla-la/{sample_id}.txt", sample_id = config["sample_downsample_fraction"]),
        expand(output_dir + "optitype/{sample_id}.txt", sample_id = config["sample_downsample_fraction"])


rule downsample:
    input:
        bam = "/home/projects/cu_10148/people/nikthu/data/1000G_ten_benchmark/{sample_id}.bam"
    output:
        bam = data_dir + "{sample_id}.bam",
        fastq1 = data_dir + "{sample_id}.bam_reads1.fastq.gz",
        fastq2 = data_dir + "{sample_id}.bam_reads2.fastq.gz",
        fastq0 = data_dir  + "{sample_id}.bam_reads_null.fastq.gz",
        singletons = data_dir  + "{sample_id}.bam_reads_singletons.fastq.gz"

    resources:
        walltime = 10000
    run:
        fraction = config["sample_downsample_fraction"][wildcards.sample_id]
        shell(script_dir + "downsample.sh {input.bam} {output.bam} {fraction}")

rule index_bam:
    input:
        bam = data_dir + "{sample_id}.bam"
    threads: 40
    resources: 
        mem = mem(gb = 120),
        walltime = 10000
    output:
        bai = data_dir + "{sample_id}.bam.bai"
    run:
        shell(script_dir + "index_bam.sh {input.bam}")

rule submit_kourami:
    input:
        bam = data_dir + "{sample_id}.bam",
        bai = data_dir + "{sample_id}.bam.bai"
    threads: 10
    resources: 
        mem = mem(gb = 50),
        walltime = 20000
    output:
        kourami = output_dir + "kourami/{sample_id}.txt",
    run:
        shell(script_dir + "tool_jobscripts/job_kourami.sh {input.bam} {wildcards.sample_id}")

rule submit_optitype:
    input:
        fastq1 = data_dir + "{sample_id}.bam_reads1.fastq.gz",
        fastq2 = data_dir + "{sample_id}.bam_reads2.fastq.gz"
    threads: 10
    resources: 
        mem = mem(gb = 50),
        walltime = 20000
    output:
        optitype = output_dir + "optitype/{sample_id}.txt"
    run:
        shell(script_dir + "tool_jobscripts/job_downsample_optitype.sh {wildcards.sample_id} " + coverage)


rule submit_hisatgenotype:
    input:
        fastq1 = data_dir + "{sample_id}.bam_reads1.fastq.gz",
        fastq2 = data_dir + "{sample_id}.bam_reads2.fastq.gz"
    threads: 10
    resources: 
        mem = mem(gb = 50),
        walltime = 20000
    output:
        hisatgenotype = output_dir + "hisatgenotype/{sample_id}.txt"
    run:
        shell(script_dir + "tool_jobscripts/job_downsample_hisatgenotype.sh {wildcards.sample_id} " + coverage)


rule submit_stc_seq:
    input:
        fastq1 = data_dir + "{sample_id}.bam_reads1.fastq.gz",
        fastq2 = data_dir + "{sample_id}.bam_reads2.fastq.gz",
        fastq0 = data_dir + "{sample_id}.bam_reads_null.fastq.gz",
        singletons = data_dir + "{sample_id}.bam_reads_singletons.fastq.gz"
    threads: 10
    resources:
        mem = mem(gb = 100),
        walltime = 20000
    output:
        stc_seq = output_dir + "stc-seq/{sample_id}.txt"
    run:
        shell(script_dir + "tool_jobscripts/job_downsample_stc-seq.sh {wildcards.sample_id} " + coverage)


rule submit_hla_la:
    input:
        bam = data_dir + "{sample_id}.bam",
        bai = data_dir + "{sample_id}.bam.bai"
    threads: 10
    resources: 
        mem = mem(gb = 185),
        walltime = 172800
    output:
        hla_la = output_dir + "hla-la/{sample_id}.txt"
    run:
        shell(script_dir + "tool_jobscripts/job_hla-la.sh {input.bam} /home/projects/cu_10148/people/nikthu/output/hla-la {wildcards.sample_id}")

