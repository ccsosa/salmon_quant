#!/bin/bash
#SBATCH --mem=20G
#SBATCH -o %j.out
##SBATCH -e %j.err
#SBATCH --cpus-per-task=1
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load SRAToolkit/3.0.1
time fastq-dump -I --split-files --gzip -M 36 SRR4947485 SRR4947493 SRR4947501 SRR4947509 SRR4947489 SRR4947497 SRR4947505 SRR4947484 SRR4947492 SRR4947500 SRR4947508 SRR4947488 SRR4947496 SRR4947504 SRR4947481�SRR4947480