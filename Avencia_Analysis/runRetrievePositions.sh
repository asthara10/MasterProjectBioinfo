#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=8000

#SBATCH -e stderr2.err
#SBATCH -o stdout2.out

module load Python/3.5.2-foss-2016b

python3 retrievePositions.py -f "results/5116-Control-of-mutated_S10_L001_R1_001.fa" "results/5116-Control-of-mutated_S7_L001_R1_001.fa" "results/5116-Mutated_S9_L001_R1_001.fa" "results/5116-Reverted_S8_L001_R1_001.fa" "results/LV-11-4-2018_S4_L001_R1_001.fa" "results/LV-20-6-2018_S1_L001_R1_001.fa" "results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001.fa" "results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001.fa" "results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001.fa" "results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001.fa" -t "results/5116-Control-of-mutated_S10_L001_R1_001.txt" "results/5116-Control-of-mutated_S7_L001_R1_001.txt" "results/5116-Mutated_S9_L001_R1_001.txt" "results/5116-Reverted_S8_L001_R1_001.txt" "results/LV-11-4-2018_S4_L001_R1_001.txt" "results/LV-20-6-2018_S1_L001_R1_001.txt" "results/NILV-SP-IN-AAVS1-11-4-2018_S5_L001_R1_001.txt" "results/NILV-SP-IN-AAVS1-20-6-2018_S2_L001_R1_001.txt" "results/NILV-SP-IN-CCR5-11-4-2018_S6_L001_R1_001.txt" "results/NILV-SP-IN-CCR5-20-6-2018_S3_L001_R1_001.txt" -l "LTR.fa" -s "mappedHuman.sai" -o "filteredNames.txt"
