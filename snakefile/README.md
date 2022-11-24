### Snakemake
Use python interpreter  
Workflows are defined in terms of rules that define how to create output files from input files. Dependencies between the rules are determined automatically, creating a DAG (directed acyclic graph) of jobs that can be automatically parallelized.  
Snakemake only re-runs jobs if one of the input files is newer than one of the output files or one of the input files will be updated by another job.  

On terminal
1. Download dataset  
`curl -L https://github.com/snakemake/snakemake-tutorial-data/archive/v5.24.1.tar.gz -o snakemake-tutorial-data.tar.gz`  
2. extract data, create folder data, and file environment, run  
`tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"` (Linux) OR  
`tar -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"` (MacOS)  
3. Create environment with required software  
`conda activate base`   
`mamba env create --name snakemake-tutorial --file environment.yaml`  
`conda install -n base -c conda-forge mamba`   
4. Activate env    
`conda activate snakemake-tutorial`  
***
example 
```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
 ```
 Here we name the rule "bwa_map", inputs are files that are expected to be used, and outputs are files to be created. Shell is followed by python strings with shell command to execute. Concatenate input if multiple files. Pipe output of first command to next command. 
 
 ```
snakemake -np mapped_reads/A.bam
 ```
 tell snakemake to generate mapped_reads/A.bam. `-n`only show execution plan, instead of performing the steps, `-p` print the resulting shell command for illustration
 ```
 snakemake --cores 1 mapped_reads/A.bam
 ```
 Specifying the number of cores to use
 
 The previous example specifies for 1 sample "A", to automate processing for all samples, use generalizing rules of named wildcards.
 Replace the A with wildcard {sample}, and yields
 ```
 rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
 ```
 use expand function when using multiple wildcards
 When using own custum script:
 ```
 rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
 ```
 Script paths are always relative to the referring Snakefile.  
 for Rscripts use `snakemake@input[["myfile"]].`
 
 Target rule:
 Rule names as targets, if requested rule does not have wildcards
 Snakemake will define the first rule of the Snakefile as the target. Hence, it is best practice to have a rule all at the top of the workflow which has all typically desired target files as input files.
```
rule all:
    input:
        "plots/quals.svg"
snakemake -n
```
Here, this means that we add a rule, to the top of our workflow. When executing Snakemake with the execution plan for creating the final file, summarization of all results

Advanced:
```
 rule:
    input:
    output:
    threads: 8
    shell:
snakemake --cores 10
```
Thread: multiple thread to speed up computation
cores: execution CPU cores, if no number is given, all available cores are used

samples to be considered can be provided by Python list in snakefile
To be customizable and adaptable to new data use config
Config can be writen in JSON or YAML and are used with the configfile directive 

To the top of snakefile add:
```
configfile: "config.yaml"

This file yields:
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
``` 
snakefile will load config file and store its content into globally available dictionary named config 
```
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
```
The expand functions in the list of input files of the rule bcftools_call are executed during the initialization phase. In this phase, we don’t know about jobs, wildcard values and rule dependencies. Hence, we cannot determine the FASTQ paths for rule bwa_map from the config file in this phase, because we don’t even know which jobs will be generated from that rule. Instead, we need to defer the determination of input files to the DAG phase. This can be achieved by specifying an input function instead of a string as inside of the input directive. For the rule bwa_map this works as follows:
```
def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```
Rule parameter:
Shell commands may need additional paprameters to be set depending on the wildcard values fo the job
For this, Snakemake allows to define arbitrary parameters for rules with the params directive
```
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    threads: 8
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"
```
Log:
Store log, organized when running parallel
```
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```
Marking output files temporally to save storage space:
```
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```
BAM is deleted once coresponding job has been executed 
If we wish to protect the final BAM from accidental del/mod
```
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```
        
 Automatic deployment of software dependencies:
 Specify isolated software env for a whole wokflow with `envs/samtools.yaml` that yields:
 ```
channels:
  - bioconda
  - conda-forge
dependencies:
  - samtools =1.9

execute snakemake with:  
snakemake --use-conda --cores 1
```
        



