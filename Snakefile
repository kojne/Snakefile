configfile: "config.yaml"

DATA_URL = "https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz"
DATA_TAR = "snakemake-tutorial-data-v5.4.5.tar.gz"

rule all:
    input:
        "plots/quals.svg"


rule download_data:
    output:
        DATA_TAR
    shell:
        "wget {DATA_URL} -O {output}"

rule extract_data:
    input:
        DATA_TAR
    output:
        directory("data")
    shell:
        "mkdir -p {output} && tar -xzvf {input} --strip-components=1 -C {output}"

rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"



rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"

