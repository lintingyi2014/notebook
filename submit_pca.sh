#!/bin/bash
input_bfile="/data/HumanRNAProject/Aging_samples/Genotype/Imputation/MGS1234/mgs1234"

python pca.py \
    --study_in "$input_bfile" \
    --format_in bfile \
    --output_file "${input_bfile}_pca" \
    --verbose \
    --min_maf 0.05 \
    --max_missingness 0.05 \
    --max_r2 0.1 \
    --window_size 200 \
    --step_size 100 \
    --plink1_path /data/HumanRNAProject/Aging_samples/Genotype/SHARDS-main/plink \
    --plink2_path /data/HumanRNAProject/Aging_samples/Genotype/SHARDS-main/plink2 \
    --underscore_sample_id
