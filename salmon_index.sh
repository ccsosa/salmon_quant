#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load salmon
time salmon index -t  /users/ccsosaa/data/transcriptome/gentrome.fa -d /users/ccsosaa/data/transcriptome/decoys.txt -p 12 -i salmon_index --gencode
 
