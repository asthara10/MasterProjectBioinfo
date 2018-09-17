#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=100000

#SBATCH -e stderr2.err
#SBATCH -o stdout2.out

module load Python/3.5.2-foss-2016b

python3 SimulateSequencing.py -f "DigestedFragments.fasta" -as "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" -al "TGCAATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" -ac "ATCTCTCTCTTAAAAAAAAAAAAAAAAAAAAAAATTGAGAGAGAT" -l "GG" "CC" "TGCAC" "ACGTG" -r "CCTGCA" "GGACGT" "G" "C" -o "SimulatedReads.fasta"

