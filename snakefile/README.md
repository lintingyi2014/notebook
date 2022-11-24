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
