#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=8000

#SBATCH -e stderr.err
#SBATCH -o stdout.out

module load Bowtie2/2.3.2-foss-2016b

bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/5116-Control-of-mutated_S10_L001_R1_001.fastq -2 results/5116-Control-of-mutated_S10_L001_R2_001.fastq -S 5116-Control-S10.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/5116-Control-of-mutated_S7_L001_R1_001.fastq -2 results/5116-Control-of-mutated_S7_L001_R2_001.fastq -S 5116-Control-S7.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/5116-Mutated_S9_L001_R1_001.fastq -2 results/5116-Mutated_S9_L001_R2_001.fastq -S 5116-Mutated-S9.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/5116-Reverted_S8_L001_R1_001.fastq -2 results/5116-Reverted_S8_L001_R2_001.fastq -S 5116-Reverted-S8.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/LV-11-4-2018_S4_L001_R1_001.fastq -2 results/LV-11-4-2018_S4_L001_R2_001.fastq -S LV-S4.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/LV-20-6-2018_S1_L001_R1_001.fastq -2 results/LV-20-6-2018_S1_L001_R2_001.fastq -S LV-S1.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001.fastq -2 results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R2_001.fastq -S NILV-S5.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001.fastq -2 results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R2_001.fastq -S NILV-S2.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001.fastq -2 results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R2_001.fastq -S NILV-S6.sai
bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001.fastq -2 results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R2_001.fastq -S NILV-S3.sai
