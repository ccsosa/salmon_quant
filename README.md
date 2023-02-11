# salmon_quant
These are bash files to run salmon quantification and multiqc


# Description:

Do the following steps:
- 0.) Preparing salmon index (salmon_index.sh)
- 1.) Download fastq files from NCBI (fastq_down.sh)
- 2.) Run fastqc to observe sequence quality (fastqc_1.sh or fastqc_1J.sh for single and paired reads respectively)
- 3.) Run Trimgalore to remove adapters and polyA tails (trimgalore_1.sh or trimgalore_1J.sh for single and paired reads respectively)
- 4.) Run fastqc again to observe sequence quality after trimming.(fastqc_2.sh)
- 5.) Run Salmon quantification step(salmon_quant_single.sh or salmon_quant.sh for single and paired reads respectively)
- 6.) Run multiqc to summarize all steps in one HTML report (multiqc.sh)


