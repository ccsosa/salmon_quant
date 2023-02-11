#!/bin/bash
#SBATCH --mem=20G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load lang/python/3.9.7
multiqc /users/ccsosaa/workdir/ --ignore "/users/ccsosaa/workdir/tmp/" "/users/ccsosaa/workdir/trim_galore/" "/users/ccsosaa/workdir/fastqc/" --outdir /users/ccsosaa/workdir/multiqc/ 

