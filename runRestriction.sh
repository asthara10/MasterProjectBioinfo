#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=10000

#SBATCH -e stderr1.err
#SBATCH -o stdout1.out

module load R/3.4.2-foss-2016b
module load Python/3.5.2-foss-2016b

Rscript --vanilla DigestGenome.R

