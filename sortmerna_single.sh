#!/bin/bash
#SBATCH --mem=40
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --partition=FULL
#SBATCH --threads-per-core=1
##SBATCH --nodelist=node15

module load sortmerna/4.3.6

for f in *_1_trimmed.fq.gz;# for each sample
do
  n=${f%%_1_trimmed.fq.gz} # strip part of file name
  sortmerna --ref /users/ccsosaa/rRNA_databases_v4.3.4/smr_v4.3_default_db.fasta  --threads 4 --reads /users/ccsosaa/workdir/trim_galore/${n}_1_trimmed.fq.gz --workdir /users/ccsosaa/workdir/tmp --aligned rRNA_reads --fastx --zip-out YES --other ${n}_1_tr_non_rRNA_reads
  mv  /users/ccsosaa/workdir/trim_galore/${n}_1_tr_non_rRNA_reads.fq.gz /users/ccsosaa/workdir/sortmerna/${n}_1_tr_non_rRNA_reads.fq.gz
  mv  /users/ccsosaa/workdir/trim_galore/rRNA_reads.fq.gz /users/ccsosaa/workdir/sortmerna/${n}_1_tr_rRNA_reads.fq.gz
  mv  /users/ccsosaa/workdir/trim_galore/rRNA_reads.log /users/ccsosaa/workdir/sortmerna/${n}_1_rRNA_reads.log
  rm -rfv /users/ccsosaa/workdir/tmp
done