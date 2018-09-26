#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=8000

#SBATCH -e stderr3.err
#SBATCH -o stdout3.out

module load SAMtools/1.6-foss-2016b

samtools view -bS 5116-Control-S10.sai > 5116-Control-S10.bam
samtools sort -o 5116-Control-S10.sorted.bam 5116-Control-S10.bam
samtools index -b 5116-Control-S10.sorted.bam 5116-Control-S10.sorted.bam.bai
samtools view -bS LV-S1.sai > LV-S1.bam
samtools sort -o LV-S1.sorted.bam LV-S1.bam
samtools index -b LV-S1.sorted.bam LV-S1.sorted.bam.bai
samtools view -bS NILV-S3.sai > NILV-S3.bam
samtools sort -o NILV-S3.sorted.bam NILV-S3.bam
samtools index -b NILV-S3.sorted.bam NILV-S3.sorted.bam.bai
samtools view -bS 5116-Control-S7.sai > 5116-Control-S7.bam
samtools sort -o 5116-Control-S7.sorted.bam 5116-Control-S7.bam
samtools index -b 5116-Control-S7.sorted.bam 5116-Control-S7.sorted.bam.bai
samtools view -bS LV-S4.sai > LV-S4.bam
samtools sort -o LV-S4.sorted.bam LV-S4.bam
samtools index -b LV-S4.sorted.bam LV-S4.sorted.bam.bai
samtools view -bS 5116-Mutated-S9.sai > 5116-Mutated-S9.bam
samtools sort -o 5116-Mutated-S9.sorted.bam 5116-Mutated-S9.bam
samtools index -b 5116-Mutated-S9.sorted.bam 5116-Mutated-S9.sorted.bam.bai
samtools view -bS 5116-Reverted-S8.sai > 5116-Reverted-S8.bam
samtools sort -o 5116-Reverted-S8.sorted.bam 5116-Reverted-S8.bam
samtools index -b 5116-Reverted-S8.sorted.bam 5116-Reverted-S8.sorted.bam.bai
samtools view -bS NILV-S6.sai > NILV-S6.bam
samtools sort -o NILV-S6.sorted.bam NILV-S6.bam
samtools index -b NILV-S6.sorted.bam NILV-S6.sorted.bam.bai
samtools view -bS NILV-S5.sai > NILV-S5.bam
samtools sort -o NILV-S5.sorted.bam NILV-S5.bam
samtools index -b NILV-S5.sorted.bam NILV-S5.sorted.bam.bai
samtools view -bS NILV-S2.sai > NILV-S2.bam
samtools sort -o NILV-S2.sorted.bam NILV-S2.bam
samtools index -b NILV-S2.sorted.bam NILV-S2.sorted.bam.bai
