#!/bin/bash
module load bcftools

# Define the input and output files
input_files="PROJ207_chr*.dose.vcf.gz"
output_file="PROJ207_MergedImputed.vcf.gz"

# Run bcftools concat with wildcard for input files
bcftools concat --allow-overlaps -o "$output_file" $input_files
