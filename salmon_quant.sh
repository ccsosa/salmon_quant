#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load salmon/1.9.0

for i in *_1_val_1_nonrrna.fq.gz
do
   n=${i%%_1_val_1_nonrrna.fq.gz} # strip part of file name

   salmon quant -i /users/ccsosaa/data/salmon_index/salmon_index \
   --libType A -o /users/ccsosaa/workdir/salmon/${n}\
   -1 ${n}_1_val_1_nonrrna.fq.gz -2 ${n}_2_val_2_nonrrna.fq.gz --softclip  --validateMappings --seqBias --gcBias --minScoreFraction 0.6
done
