#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load salmon/1.9.0

for i in *_1_trimmed.fq.gz
do
   n=${i%%_1_trimmed.fq.gz} # strip part of file name

   salmon quant -i /users/ccsosaa/data/salmon_index/salmon_index \
   --libType A -o /users/ccsosaa/workdir/salmon/${n}\
   --seqBias --gcBias --validateMappings \
   -r ${n}_1_trimmed.fq.gz
done
