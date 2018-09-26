#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=8000

#SBATCH -e stderr.err
#SBATCH -o stdout.out

bash runRetrievePositions.sh
bash filterFastq.sh
bash runSamTools.sh
bash runBowtie2.sh
