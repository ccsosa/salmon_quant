#!/bin/bash
#SBATCH --mem=40
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=4
#SBATCH --partition=FULL
#SBATCH --threads-per-core=1
##SBATCH --nodelist=node15

module load lang/python/3.9.7

ribodetector_cpu -t 20 \
  -l 40 \
  -i /users/ccsosaa/workdir/trim_galore/SRR4947480_1_trimmed.fq.gz \
  -e rrna \
  --chunk_size 256 \
  -o  /users/ccsosaa/workdir/ribodetector/SRR4947480_1_nonrrna.fq.gz