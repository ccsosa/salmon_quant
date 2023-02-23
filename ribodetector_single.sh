#!/bin/bash
#SBATCH --mem=40
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=4
#SBATCH --partition=FULL
#SBATCH --threads-per-core=1
##SBATCH --nodelist=node15

module load lang/python/3.9.7

for f in *_1_trimmed.fq.gz;# for each sample

do
	n=${f%%_1_trimmed.fq.gz} # strip part of file name
	  ribodetector_cpu -t 4 \
		  -l 100 \
		  -i /users/ccsosaa/workdir/trim_galore/${n}_1_trimmed.fq.gz \
		  -e rrna \
		  --chunk_size 256 \
		  -o  /users/ccsosaa/workdir/ribodetector/${n}_1_trimmed_nonrrna.fq.gz
done