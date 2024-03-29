### VCF formatting and merging (individual VCF -> combined VCF)

### Change linebreaks to next line
```
:%s/\r/\r/g 
```
### Replace 'Integer' with 'Float' in all files
```
grep -RiIl 'Integer' | xargs sed -i 's/Integer/Float/g'
```
### Remove all INFO fields except FORMAT/GT from VCF
```
bcftools annotate -x ^FORMAT/GT hg38_205271050012_R01C01.vcf 
```
### Loop through all VCF files, remove INFO fields except FORMAT/GT, and save as separate files
```
for i in *.vcf; do bcftools annotate -x ^FORMAT/GT $i -o ${i%.vcf}_GT.vcf; done 
```
### Parallel compression of VCF files using bgzip
```
parallel bgzip {}::: *_GT.vcf
```
### Parallel indexing of compressed VCF files using tabix
```
parallel tabix {} ::: *_GT.vcf.gz
```
### Check the number of variants in each VCF, ensuring the merged file has more variants
### (Assumption: #variants in merged file > #variants in each individual VCF)
### This issue has been solved

### Decompress the merged VCF file
```
bgzip -d Merged_184.vcf.gz 
```
### Filter out a specific individual (e.g., 185:205271050012_R01C01) from the merged VCF
```
vcftools --remove-indv 185:205271050012_R01C01 --vcf Merged_184.vcf --recode --out Merged_184_filtered.vcf
```
### Normalize the VCF file, splitting multiallelic variants
```
bcftools norm --multiallelics -any Merged_184_filtered.vcf.recode.vcf -o Merged_184_split_multiallelics.vcf
```
