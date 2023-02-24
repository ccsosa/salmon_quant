#!/bin/bash
#SBATCH --mem=40
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=1
#SBATCH --ntasks=6
#SBATCH --partition=FULL
#SBATCH --threads-per-core=1
##SBATCH --nodelist=node15

module load sortmerna/4.3.6

for f in *_1_val_1.fq.gz;# f
do
  n=${f%%_1_val_1.fq.gz} # strip part of file name
  sortmerna --ref /users/ccsosaa/rRNA_databases_v4.3.4/smr_v4.3_default_db.fasta  --threads 6 --reads /users/ccsosaa/workdir/trim_galore/${n}_1_val_1.fq.gz --reads /users/ccsosaa/workdir/trim_galore/${n}_2_val_2.fq.gz --workdir /users/ccsosaa/workdir/tmp --aligned rRNA_reads --fastx --zip-out YES --paired_in --other non_rRNA_reads --out2
  
  mv /users/ccsosaa/workdir/trim_galore/non_rRNA_reads_fwd.fq.gz /users/ccsosaa/workdir/sortmerna/${n}_1_val_1_nonrrna.fq.gz
  mv /users/ccsosaa/workdir/trim_galore/non_rRNA_reads_rev.fq.gz /users/ccsosaa/workdir/sortmerna/${n}_2_val_2_nonrrna.fq.gz

  mv /users/ccsosaa/workdir/trim_galore/rRNA_reads_fwd.fq.gz /users/ccsosaa/workdir/sortmerna/${n}_1_val_1_rrna.fq.gz
  mv /users/ccsosaa/workdir/trim_galore/rRNA_reads_rev.fq.gz /users/ccsosaa/workdir/sortmerna/${n}_2_val_2_rrna.fq.gz
 
  mv /users/ccsosaa/workdir/trim_galore/rRNA_reads.log /users/ccsosaa/workdir/sortmerna/${n}_rRNA_reads.log

  rm -rfv /users/ccsosaa/workdir/tmp
done