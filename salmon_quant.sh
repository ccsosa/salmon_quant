#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load salmon/1.9.0

for i in *_1_val_1.fq.gz
do
   n=${i%%_1_val_1.fq.gz} # strip part of file name

   salmon quant -i /users/ccsosaa/data/salmon_index/salmon_index \
   --libType A -o /users/ccsosaa/workdir/salmon/${n}\
   --seqBias --gcBias --validateMappings \
   -1 ${n}_1_val_1.fq.gz -2 ${n}_2_val_2.fq.gz 
done
