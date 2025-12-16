#!/bin/bash

#SBATCH --job-name=dada2_hickman
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH -t 7-00:00:00
#SBATCH --mem=350G
#SBATCH -n 3



source /home/rbuecki/.guix-profile/etc/profile
Rscript dada_2_hickman_2024.R 
# Rscript dada_2_hickman_2024_parallel.R 
