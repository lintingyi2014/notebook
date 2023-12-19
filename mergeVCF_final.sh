#!/bin/bash
module load bcftools

# Define the input and output files
input_aging="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/20231016_imputation/AGING/chr*.dose.vcf.gz"
output_aging="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/20231016_imputation/Aging_MergedImputed.vcf.gz"

input_proj207="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/20231016_imputation/PROJ207/chr*.dose.vcf.gz"
output_proj207="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/20231016_imputation/Proj207_MergedImputed.vcf.gz"

input_proj236="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/20231016_imputation/PROJ236/chr*.dose.vcf.gz"
output_proj236="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/20231016_imputation/Proj236_MergedImputed.vcf.gz"

# Run bcftools concat with wildcard for input files
bcftools concat --allow-overlaps -o "$output_aging" $input_aging
bcftools concat --allow-overlaps -o "$output_proj207" $input_proj207
bcftools concat --allow-overlaps -o "$output_proj236" $input_proj236

# Define the paths to the input VCF files
vcf1=$(echo "$output_aging")
vcf2=$(echo "$output_proj207")
vcf3=$(echo "$output_proj236")

# Define the output file name for the merged VCF
output_final="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/20231016_imputation/MGS1234_MergedImputed.vcf.gz"

# Step 1: Create tabix indexes for individual input VCF files
tabix -p vcf "$vcf1"
tabix -p vcf "$vcf2"
tabix -p vcf "$vcf3"

echo "Current working directory: $(pwd)"
echo "VCF1 Index File: $vcf1.tbi"
echo "VCF2 Index File: $vcf2.tbi"
echo "VCF3 Index File: $vcf3.tbi"

# Step 2: Merge the three VCF files directly
bcftools merge "$vcf1" "$vcf2" "$vcf3" -o "$output_final"

# Step 3: Create tabix index for final merged VCF
tabix -p vcf "$output_final"

echo "Merged VCF file is saved as $output_final"
echo "Merged VCF Index File: $output_final.tbi"










