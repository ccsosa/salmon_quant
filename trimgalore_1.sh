#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err
#SBATCH --cpus-per-task=2
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load TrimGalore/0.6.7

for f in *_1.fastq.gz;# for each sample

do
	n=${f%%_1.fastq.gz} # strip part of file name
    	trim_galore -q 30 --gzip -o /users/ccsosaa/workdir/trim_galore/  ${n}_1.fastq.gz
done
