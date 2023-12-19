#!/bin/bash

# Get user input
echo "Enter the VCF file path:"
read vcf_file

echo "Enter the project name:"
read project_name

# Module loading
module load bcftools
module load python

# Get the number of samples from the VCF file
n=$(bcftools query -l "$vcf_file" | wc -l)

# Print the number of samples
echo "Number of samples: $n"

# Calculate MONO_cutoff using the value of n
MONO_cutoff=$(echo "scale=6; 1 / (2.1 * $n)" | bc)

# Print the calculated MONO_cutoff
echo "MONO_cutoff: $MONO_cutoff"

# Step 1: Convert VCF to binary PLINK format
plink2 --vcf "$vcf_file" --double-id --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10000 --make-bed --out "$project_name"

# Step 2: Calculate the Missingness Rate
plink2 --bfile "$project_name" --missing sample-only --out "$project_name"
awk '$5 > 0.05' "$project_name".smiss > outputID.smiss.txt
cat  outputID.smiss.txt

# Step 3: Ancestry Check
python /data/HumanRNAProject/Aging_samples/Genotype/SHARDS-main/src/shards.py \
  --study_vcf "$vcf_file" \
  --reference_vcf /data/HumanRNAProject/Aging_samples/Genotype/1000G_meei/CCDG_14151_B01_GRM_WGS_2020-08-05_all_chr.filtered.shapeit2-duohmm-phased.biallelic.vcf.gz \
  --reference_populations /data/HumanRNAProject/Aging_samples/Genotype/SHARDS-main/data/1000_genomes_ancestry_hg38.tsv \

# Step 4: Sex Check
# input sex information, check .fam file for FID and IID format
plink --bfile "$project_name"_nomiss --update-sex sex_check.txt --make-bed -out "$project_name"_nomiss
plink --bfile "$project_name"_nomiss --check-sex --out "$project_name"_nomiss
awk '/PROBLEM/' "$project_name"_nomiss.sexcheck > outputID.sexcheck.txt
cat outputID.sexcheck.txt

# Step 5: Heterozygosity Test
plink --bfile "$project_name"_nomiss --chr 1-22 --het --homozyg-kb 1000 --out "$project_name"_nomiss
awk '{ print $0, (NR>1) ? $5-$4 : "subsHet" }' "$project_name"_nomiss.het > "$project_name"_nomiss.het.v2
awk '{ print $0, (NR>1) ? $7/$5 : "meanHet" }' "$project_name"_nomiss.het.v2 > "$project_name"_nomiss.het.v3
cat "$project_name"_nomiss.het.v3 | head -n 10
awk '{s+=$8; ss+=$8^2} END{print m=s/NR, sqrt(ss/NR-m^2)}' "$project_name"_nomiss.het.v3
cat  "$project_name"_nomiss.het.v3

# Step 6: IBD Calculation
plink --bfile "$project_name"_nomiss --indep-pairwise 200 100 0.1 --out "$project_name"_nomiss_LD_pruned
plink --bfile "$project_name"_nomiss --extract "$project_name"_nomiss_LD_pruned.prune.in --genome --out "$project_name"_nomiss_LD_pruned
awk 'NR==1 || $9 > 0.1875' "$project_name"_nomiss_LD_pruned.genome > outputID.ibdcheck.txt
awk 'NR==FNR {ids[$1]; next} $2 in ids' IBD_pairs_callrate.txt "$project_name".smiss > IBD_pairs_callrate_FMISS.txt

# Step 7: Hard Call Variant QC
# Split chrX into nonPAR and PAR regions
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT --split-par b38 --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR

# Remove variants with missingness rate greater than 2%
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR --geno 0.02 --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss

# Check allele frequency
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss --freq --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq

# Remove monomorphic sites 
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq --maf "$MONO_cutoff" --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq_mono

# Rename variant IDs
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq_mono --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 1222 --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq_mono_newid

# Remove indels >50 bps
awk 'NR>1 {print $2, ((length($5)>length($6)?length($5):length($6)))}' "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq_mono_newid.bim > indel_len.tsv
sort -k2,2 -n -r indel_len.tsv | awk '$2 > 50' > indel_output.tsv
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq_mono_newid --exclude indel_output.tsv --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq_mono_newid_indel

# Remove mitochondrial variants 
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove --not-chr MT --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT

# Split PAR+Autosomes and nonPAR for each subpopulation
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT -chr 1-22,PAR1,PAR2 --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_Autosomes_PAR
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT -chr X --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_X

# Assign heterozygous sites on the non-PAR region of chrX in males to missing.
# First of all, we need to extract non-PAR
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT --keep-males --chr X --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_Male_chrX
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_Male_chrX --set-hh-missing --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_Male_chrX_sethhmissing
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT --remove-males --chr X --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_Female_chrX

# Step 8: Merge male and female samples
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_splitPAR_miss_freq_mono_newid_indel_remove --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales --merge-list merge_list.txt --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR
awk 'NR>1 {print $1, $2}' "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR.missnp > missing_snps.txt
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR --exclude missing_snps.txt --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_missingSNPs

# Merge with PAR regions
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT --chr PAR1,PAR2 --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_PAR
plink --merge-list mergelist.txt --allow-extra-chr --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR

# Merge with nonPAR regions
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT -chr 1-22 --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_freq_mono_newid_indelremove_removeMT_Autosome
plink --merge-list mergelist2.txt --allow-extra-chr --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr
plink2 --bfile project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr --rm-dup force-first --make-bed --out project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup

# Step 9: Remove duplicate samples
awk '{ print $1, $2 }' "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_missingSNPs.fam | sort | uniq -d > duplicate_samples.txt
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_missingSNPs --remove duplicate_samples.txt --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup

# Step 10: Calculate Hardy-Weinberg Equilibrium (HWE)
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup --hardy --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE
awk '$9 < 0.00001' "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE.hwe > low_hwe_snps.txt
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup --exclude low_hwe_snps.txt --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE

# Step 13: Extract Autosomes, PAR, and ChrX
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_missingSNPs --chr 1-22 --hardy midp --out HWE_Auto
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_missingSNPs --chr PAR1,PAR2 --hardy midp --out HWE_PAR
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_missingSNPs --chr X --hardy midp --out HWE_X

# Check HWE p-values
awk 'NR>1 && $10 < 1e-6' HWE_Auto.hardy | wc -l
awk 'NR>1 && $10 < 1e-6' HWE_PAR.hardy | wc -l
awk 'NR>1 && $14 < 1e-6' HWE_X.hardy | wc -l

# Extract problematic SNPs for HWE
awk 'NR>1 && $10 < 1e-6' HWE_Auto.hardy > HWE_Auto_fail.txt
awk 'NR>1 && $10 < 1e-6' HWE_PAR.hardy > HWE_PAR_fail.txt
awk 'NR>1 && $14 < 1e-6' HWE_X.hardy > HWE_X_fail.txt

# Merge problematic SNPs
cat HWE_Auto_fail.txt > HWE_fail_merged.txt 
cat HWE_PAR_fail.txt >> HWE_fail_merged.txt
cat HWE_X_fail.txt >> HWE_fail_merged.txt

# Remove problematic SNPs
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_missingSNPs --exclude HWE_fail_merged.txt --make-bed --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_cleaned

# Convert cleaned dataset to VCF
plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_cleaned --recode vcf id-paste=iid --out "$project_name"_cleaned

# Sort the VCF file
bcftools sort "$project_name"_cleaned.vcf -Oz -o "$project_name"_cleaned_sorted.vcf.gz

# Index the cleaned and sorted VCF
tabix -p vcf "$project_name"_cleaned_sorted.vcf.gz

# Step 14: Remove missingness rate > 2%
plink2 --vcf "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge.vcf.gz --geno 0.02 --recode vcf --out "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_miss_vcf
bgzip -c project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_miss_vcf.vcf >project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_miss_vcf.vcf.gz
bcftools annotate --rename-chrs /data/HumanRNAProject/Aging_samples/Genotype/hg38_to_vcf_mapping.txt -o project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_miss_vcf_renamechr.vcf.gz -O z project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_miss_vcf.vcf.gz
tabix -p vcf project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_miss_vcf_renamechr.vcf.gz

bcftools index -s project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_miss_vcf_renamechr.vcf.gz | cut -f 1 | while read C; do bcftools view -O z -o project_207_hg38_SplitByChr.${C}.vcf.gz project_207_hg38_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_miss_vcf_renamechr.vcf.gz -r "${C}"; done

zless project_207_hg38_SplitByChr.chrX.vcf.gz | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | bgzip -c > project_207_hg38_SplitByChr.chrX_4.2.vcf.gz

# Generate a VCF for each chromosome for TOP
for i in {1..22}
do
  plink2 --bfile "$project_name"_nomiss_noIBD_splitPAR_miss_mono_newid_indelremove_removeMT_sethhmissing_mergeMalesFemales_mergePAR_mergeChr_rmdup_HWE_merge_missingSNPs --chr $i --recode vcf --out "$project_name"_TOP_chr"$i"
done
