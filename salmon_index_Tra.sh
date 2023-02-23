#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load salmon/1.9.0
time salmon index -i transcripts_index --kmerLen 17 --tmpdir /users/ccsosaa/workdir/tmp --keepFixedFasta --transcripts /users/ccsosaa/data/transcriptome/IRGSP-1.0_cds_2022-09-01.fasta 
