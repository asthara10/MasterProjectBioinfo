#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=8000

#SBATCH -e stderr4.err
#SBATCH -o stdout4.out

module load Bowtie2/2.3.2-foss-2016b
module load SAMtools/1.6-foss-2016b

bowtie2 -x humanRefGen/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 results/LV-20-6-2018_S1_L001_R1_001.fastq -2 results/LV-20-6-2018_S1_L001_R2_001.fastq -S LV-S1.sai

samtools view -bS LV-S1.sai > LV-S1.bam
samtools sort -o LV-S1.sorted.bam LV-S1.bam
samtools index -b LV-S1.sorted.bam LV-S1.sorted.bam.bai
