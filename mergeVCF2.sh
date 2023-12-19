#!/bin/bash

# Define the input and output files
input_files="Aging_chr*.dose.vcf.gz"
output_file="Aging_MergedImputed.vcf.gz"

# Run bcftools concat with wildcard for input files
bcftools concat --allow-overlaps -o "$output_file" $input_files
