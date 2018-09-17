#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=100000

#SBATCH -e stderr.err
#SBATCH -o stdout.out

module load R/3.4.2-foss-2016b
module load Python/3.5.2-foss-2016b

Rscript --vanilla DigestGenome.R

