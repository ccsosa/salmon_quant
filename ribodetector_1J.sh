#!/bin/bash
#SBATCH --mem=40
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=4
#SBATCH --partition=FULL
#SBATCH --threads-per-core=1
##SBATCH --nodelist=node15

module load lang/python/3.9.7

for f in *_1_val_1.fq.gz;# for each sample

do
	n=${f%%_1_val_1.fq.gz} # strip part of file name
	  ribodetector_cpu -t 4 \
		  -l 100 \
		  -i /users/ccsosaa/workdir/trim_galore/${n}_1_val_1.fq.gz /users/ccsosaa/workdir/trim_galore/${n}_2_val_2.fq.gz\
		  -e rrna \
		  --chunk_size 256 \
		  -o  /users/ccsosaa/workdir/ribodetector/${n}_1_val_1_nonrrna.fq.gz /users/ccsosaa/workdir/ribodetector/${n}_2_val_2_nonrrna.fq.gz
done