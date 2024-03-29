### Request interactive session with 5 CPUs
```
sinteractive --cpus-per-task=5
```

### Perform liftover of VCF file from b37 to hg38
```
module load picard
java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar LiftoverVcf \
     I=input.vcf \
     O=lifted_over.vcf \
     CHAIN=b37tohg38.chain \
     REJECT=rejected_variants.vcf \
     R=reference_sequence.fasta
```
### Exit the interactive session
```
Exit
```
### Perform liftover of another VCF file from hg19 to hg38
```
java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar LiftoverVcf \
     I=project_207.vcf.gz\
     O=project_207_hg38.vcf.gz \
     CHAIN=/data/HumanRNAProject/Aging_samples/Genotype/hg19ToHg38.over.chain.gz \
     REJECT=rejected_variants.vcf \
     R=/data/HumanRNAProject/Aging_samples/Genotype/hg38.fa
```
### Copy necessary files to the Genotype directory
```
cp /fdb/genomebrowser/fasta/hg38.fa /data/HumanRNAProject/Aging_samples/Genotype/
cp /fdb/genomebrowser/fasta/hg38.dict /data/HumanRNAProject/Aging_samples/Genotype/
cp /fdb/genomebrowser/fasta/hg38.fa.fai /data/HumanRNAProject/Aging_samples/Genotype/
```
### Copy necessary files for hg19 to the Genotype directory
```
cp /fdb/genomebrowser/fasta/hg19.fa /data/HumanRNAProject/Aging_samples/Genotype/
cp /fdb/genomebrowser/fasta/hg19.dict /data/HumanRNAProject/Aging_samples/Genotype/
cp /fdb/genomebrowser/fasta/hg19.fa.fai /data/HumanRNAProject/Aging_samples/Genotype/
```
### Copy dbSNP files to the Genotype directory
```
/fdb/dbSNP/organisms/human_9606_b150_GRCh37p13/00-All.vcf.gz /data/HumanRNAProject/Aging_samples/Genotype/
/fdb/dbSNP/organisms/human_9606_b150_GRCh37p13/00-All.vcf.gz.tbi /data/HumanRNAProject/Aging_samples/Genotype/
```
### Copy additional directories
```
cp -r human_9606_b153_GRCh38p12 /data/HumanRNAProject/Aging_samples/Genotype/
cp -r human_9606_b153_GRCh37p13 /data/HumanRNAProject/Aging_samples/Genotype/
```
### Compress and index VCF file
```
bgzip -c project_207.vcf > project_207.vcf.gz
tabix -p vcf hg38.fa
```
### Annotate VCF file with chromosome renaming
```
bcftools annotate --rename-chrs /data/HumanRNAProject/Aging_samples/Genotype/hg38_to_vcf_mapping.txt -o project_207_config.vcf.gz -O z project_207.vcf.gz
```
### Download assembly reports for renaming chromosomes
```
report_dir='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405'
wget -N "${report_dir}/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt"
wget -N "${report_dir}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
```
### Extract useful columns from assembly reports
```
for k in *assembly_report.txt
do
    out=$(echo $k | sed 's/.txt/.chrnames/')
    grep -e '^[^#]' $k | awk '{ print $7, $11 }' > $out
done
```
### Extract alternative useful columns from assembly reports
```
for k in *assembly_report.txt
do
    out=$(echo $k | sed 's/.txt/.chrnames2/')
    grep -e '^[^#]' $k | awk '{ print $1, $11 }' > $out
done
```
### Annotate VCF file with renamed chromosomes using assembly reports
```
bcftools annotate \
  --rename-chrs /data/HumanRNAProject/Aging_samples/Genotype/GCF_000001405.25_GRCh37.p13_assembly_report.chrnames \
  --threads 10 -Oz \
  -o GRCh37.dbSNP153.vcf.gz \
  /fdb/dbSNP/organisms/human_9606_b153_GRCh37p13/VCF/GCF_000001405.25.gz

bcftools annotate \
  --rename-chrs /data/HumanRNAProject/Aging_samples/Genotype/GCF_000001405.39_GRCh38.p13_assembly_report.chrnames \
  --threads 10 -Oz \
  -o GRCh38.dbSNP153.vcf.gz \
  /fdb/dbSNP/organisms/human_9606_b153_GRCh38p12/VCF/GCF_000001405.38.gz
```
### Index annotated VCF files
```
gatk IndexFeatureFile -I GRCh37.dbSNP153.vcf.gz
gatk IndexFeatureFile -I GRCh38.dbSNP153.vcf.gz
```
### Annotate original VCF file with another set of renamed chromosomes
```
bcftools annotate \
  --rename-chrs /data/HumanRNAProject/Aging_samples/Genotype/GCF_000001405.25_GRCh37.p13_assembly_report.chrnames2 \
  --threads 10 -Oz \
  -o project_207_config2.vcf.gz \
  project_207.vcf.gz 
```
### Index annotated VCF file
```
gatk IndexFeatureFile -I project_207_config2.vcf.gz
tabix -p vcf project_207_config2.vcf.gz
```
### Validate variants against dbSNP for GRCh37
```
module load GATK
gatk ValidateVariants \
   -R /data/HumanRNAProject/Aging_samples/Genotype/hg19.fa \
   -V /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/fixed_project_207_config2.vcf.gz \
   --dbsnp /data/HumanRNAProject/Aging_samples/Genotype/GRCh37.dbSNP153.vcf.gz
```
### Validate variants against dbSNP for GRCh38
```
gatk ValidateVariants \
   -R /data/HumanRNAProject/Aging_samples/Genotype/hg19.fa \
   -V /data/HumanRNAProject/Aging_samples/Genotype/Rinki_genotype/project_207_config.vcf.gz\
   --dbsnp /fdb/dbSNP/organisms/human_9606_b153_GRCh38p12/VCF/GCF_000001405.38.gz
```
### Perform another liftover to hg38
```
java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar LiftoverVcf \
     I=fixed_project_207_config2.vcf.gz\
     O=project_207_hg38.vcf.gz \
     CHAIN=/data/HumanRNAProject/Aging_samples/Genotype/hg19ToHg38.over.chain.gz \
     REJECT=rejected_variants.vcf \
     R=/data/HumanRNAProject/Aging_samples/Genotype/hg38.fa
```
