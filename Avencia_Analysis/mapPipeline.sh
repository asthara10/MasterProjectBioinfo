#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=8000

#SBATCH -e stderr.err
#SBATCH -o stdout.out

module load Python/3.5.2-foss-2016b
module load BBMap
module load Bowtie2/2.3.2-foss-2016b
module load SAMtools/1.6-foss-2016b

### Get sequence names which first 24 nucleotides form part of the LTR

python3 retrievePositions.py -f "results/LV-11-4-2018_S4_L001_R1_001.fa" "results/LV-20-6-2018_S1_L001_R1_001.fa" "results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001.fa" "results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001.fa" "results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001.fa" "results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001.fa" -t "results/LV-11-4-2018_S4_L001_R1_001.txt" "results/LV-20-6-2018_S1_L001_R1_001.txt" "results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001.txt" "results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001.txt" "results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001.txt" "results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001.txt" -l "LTR.fa" -s "mappedHuman.sai" -o "filteredNames.txt"

### Filter fastq files including only the sequences that contain the LRT

filterbyname.sh in=results/LV-11-4-2018_S4_L001_R1_001.fastq out=results/LV-11-4-2018_S4_L001_R1_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/LV-11-4-2018_S4_L001_R2_001.fastq out=results/LV-11-4-2018_S4_L001_R2_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/LV-20-6-2018_S1_L001_R1_001.fastq out=results/LV-20-6-2018_S1_L001_R1_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/LV-20-6-2018_S1_L001_R2_001.fastq out=results/LV-20-6-2018_S1_L001_R2_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001.fastq out=results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001.fastq out=results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R2_001.fastq out=results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R2_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R2_001.fastq out=results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R2_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001.fastq out=results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R2_001.fastq out=results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R2_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001.fastq out=results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001_filtered.fastq names=filteredNames.txt include=t
filterbyname.sh in=results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R2_001.fastq out=results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R2_001_filtered.fastq names=filteredNames.txt include=t

### Trimming the LTR sequence from each read

python3 trimSequences.py -f "results/LV-11-4-2018_S4_L001_R1_001_filtered.fastq" "results/LV-20-6-2018_S1_L001_R1_001_filtered.fastq" "results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001_filtered.fastq" "results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001_filtered.fastq" "results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001_filtered.fastq" "results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001_filtered.fastq" -fw "56" -rv "62" -o "forward"
python3 trimSequences.py -f "results/LV-11-4-2018_S4_L001_R2_001_filtered.fastq" "results/LV-20-6-2018_S1_L001_R2_001_filtered.fastq" "results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R2_001_filtered.fastq" "results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R2_001_filtered.fastq" "results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R2_001_filtered.fastq" "results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R2_001_filtered.fastq" -fw "56" -rv "62" -o "reverse"

### Perform the mapping

bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/LV-11-4-2018_S4_L001_R1_001_filtered_trim.fastq -2 results/LV-11-4-2018_S4_L001_R2_001_filtered_trim.fastq -S LV-S4.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/LV-20-6-2018_S1_L001_R1_001_filtered_trim.fastq -2 results/LV-20-6-2018_S1_L001_R2_001_filtered_trim.fastq -S LV-S1.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001_filtered_trim.fastq -2 results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R2_001_filtered_trim.fastq -S NILV-S5.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001_filtered_trim.fastq -2 results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R2_001_filtered_trim.fastq -S NILV-S2.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001_filtered_trim.fastq -2 results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R2_001_filtered_trim.fastq -S NILV-S6.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001_filtered_trim.fastq -2 results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R2_001_filtered_trim.fastq -S NILV-S3.sai

### Running samtools to obatin bam files

samtools view -bS LV-S1.sai > LV-S1.bam
samtools sort -o LV-S1.sorted.bam LV-S1.bam
samtools index -b LV-S1.sorted.bam LV-S1.sorted.bam.bai
samtools view -bS NILV-S3.sai > NILV-S3.bam
samtools sort -o NILV-S3.sorted.bam NILV-S3.bam
samtools index -b NILV-S3.sorted.bam NILV-S3.sorted.bam.bai
samtools view -bS LV-S4.sai > LV-S4.bam
samtools sort -o LV-S4.sorted.bam LV-S4.bam
samtools index -b LV-S4.sorted.bam LV-S4.sorted.bam.bai
samtools view -bS NILV-S6.sai > NILV-S6.bam
samtools sort -o NILV-S6.sorted.bam NILV-S6.bam
samtools index -b NILV-S6.sorted.bam NILV-S6.sorted.bam.bai
samtools view -bS NILV-S5.sai > NILV-S5.bam
samtools sort -o NILV-S5.sorted.bam NILV-S5.bam
samtools index -b NILV-S5.sorted.bam NILV-S5.sorted.bam.bai
samtools view -bS NILV-S2.sai > NILV-S2.bam
samtools sort -o NILV-S2.sorted.bam NILV-S2.bam
samtools index -b NILV-S2.sorted.bam NILV-S2.sorted.bam.bai

