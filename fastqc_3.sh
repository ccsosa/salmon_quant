#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

##https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/06_qc_running_fastqc_sbatch.html
module load fastqc/0.11.9

fastqc -d /users/ccsosaa/workdir/tmp/ -o /users/ccsosaa/workdir/fastqc_ribodetector/ -t 6 *.fq.gz
