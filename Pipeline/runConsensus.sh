#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.mir01@estudiant.upf.edu
#SBATCH --mem=8000

#SBATCH -e stderr4.err
#SBATCH -o stdout4.out

module load R/3.4.2-foss-2016b
module load GSL/2.1-foss-2016b
module load MAFFT/7.305-foss-2016b-with-extensions

Rscript --vanilla CreateConsensus.R

