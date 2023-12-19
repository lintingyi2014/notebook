### Generate statistics for the original VCF file
```
bcftools stats project_207.vcf
```
### Load Python module
```
module load python
```
### Fix ALT field in VCF file using a Python script
### Note: Ensure the reference file is indexed using 'samtools faidx reference.fa'
```
python fix_vcf4.py -f hg19.fa -o /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207.vcf /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/project_207.vcf
```
### Quick workaround for 'Contig '1' is not defined in the header' issue
```
bgzip -c project_207.vcf > project_207.vcf.gz
tabix -p vcf project_207.vcf.gz
```
### Fix ALT field in the compressed VCF file
```
python fix_vcf4.py -f hg19.fa -o /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207.vcf.gz /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/project_207.vcf.gz
```
### Handle KeyError caused by mismatch between '1' in fasta and 'chr1' in VCF
### Convert fasta 'chr1' to VCF '1'
### Download assembly report for chromosome renaming
```
report_dir='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405'
wget -N "${report_dir}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
```
### Extract useful columns from assembly report
```
for k in *assembly_report.txt
do
    out=$(echo $k | sed 's/.txt/.chrnames2/')
    grep -e '^[^#]' $k | awk '{ print $1, $11 }' > $out
done
```
### Rename chromosomes in VCF file using assembly report
```
bcftools annotate \
  --rename-chrs /data/HumanRNAProject/Aging_samples/Genotype/GCF_000001405.25_GRCh37.p13_assembly_report.chrnames2 \
  --threads 10 -Oz \
  -o project_207_config2.vcf.gz \
  project_207.vcf.gz 
```
### Index the annotated VCF file
```
tabix -p vcf project_207_config2.vcf.gz
```
### FixTypeError in Python script
```
python fix_vcf6.py -f hg19.fa -o /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207_config2.vcf.gz /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/project_207_config2.vcf.gz
```
### Validate variants using GATK
```
gatk ValidateVariants \
   -R /data/HumanRNAProject/Aging_samples/Genotype/hg19.fa \
   -V /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207_config2.vcf.gz \
   --dbsnp /data/HumanRNAProject/Aging_samples/Genotype/GRCh37.dbSNP153.vcf.gz
```
### Check file integrity of the compressed VCF
```
gunzip -t /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207_config2.vcf.gz
```
### Rename and compress the VCF file
```
mv fixALT_project207_config2.vcf.gz fixALT_project207_config2.vcf
bgzip /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207_config2.vcf
```
### Index the compressed and renamed VCF file
```
tabix -p vcf fixALT_project207_config2.vcf.gz
```
### Sort the VCF file
```
bcftools sort -o fixALT_project207_config2_sorted.vcf.gz -O z fixALT_project207_config2.vcf.gz
```
### Validate variants after sorting
```
gatk ValidateVariants \
   -R /data/HumanRNAProject/Aging_samples/Genotype/hg19.fa \
   -V /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207_config2_sorted.vcf.gz \
   --dbsnp /data/HumanRNAProject/Aging_samples/Genotype/GRCh37.dbSNP153.vcf.gz
```
### Resolve unsorted positions issue
```
bcftools view -r 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX' -O z -o fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz  fixALT_project207_config2_sorted.vcf.gz 
```
### Index the sorted VCF file
```
tabix -p vcf fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz
```
### Validate variants after sorting and filtering
```
gatk ValidateVariants \
   -R /data/HumanRNAProject/Aging_samples/Genotype/hg19.fa \
   -V /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz  \
   --dbsnp /data/HumanRNAProject/Aging_samples/Genotype/GRCh37.dbSNP153.vcf.gz
```
### Display lines containing chromosome codes in the filtered VCF
```
grep -n -E '^#|chr[0-9XYM]+$' fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz 
```
### Check file integrity of the gzipped VCF
```
gzip -d fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz
```
### Compress the VCF file
```
bgzip -c fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf > fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz
```
### Generate statistics for the final VCF
```
bcftools stats fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz
```
### Perform liftover to hg38
```
java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar LiftoverVcf \
     I=fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz\
     O=project_207_hg38.vcf.gz \
     CHAIN=/data/HumanRNAProject/Aging_samples/Genotype/hg19ToHg38.over.chain.gz \
     REJECT=rejected_variants.vcf \
     R=/data/HumanRNAProject/Aging_samples/Genotype/hg38.fa
```
### Handle exception caused by malformed VCF record with allele '-'
### Display the problematic VCF record
```
zcat  fixALT_project207_config2_sorted_ChrAuto_ChrX.vcf.gz| head -n 726 | tail -n 1
```
# Fix REF field containing '-'
```
python fix_vcf7.py
```
