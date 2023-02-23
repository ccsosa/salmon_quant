# salmon_quant
These are bash files to run salmon quantification and multiqc

##Steps to run decoy aware salmon index
```
#https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
grep "^>" < IRGSP-1.0_genome.fasta | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat IRGSP-1.0_cds_2022-09-01.fasta IRGSP-1.0_genome.fasta > gentrome.fa
```

# Description:

Do the following steps:
- 0.) Preparing salmon index (salmon_index.sh to use decoy-aware option or salmon_index_Tra.sh  to use only the transcriptome file)
- 1.) Download fastq files from NCBI (fastq_down.sh)
- 2.) Run fastqc to observe sequence quality (fastqc_1.sh or fastqc_1J.sh for single and paired reads respectively)
- 3.) Run Trimgalore to remove adapters and polyA tails (trimgalore_1.sh or trimgalore_1J.sh for single and paired reads respectively)
- 4.) Run fastqc again to observe sequence quality after trimming.(fastqc_2.sh)
- 4.1) optional: Run Ribodetector to reduce the rRNA contamination (ribodetector_single.sh or ribodetector_1J.sh for single and paired reads respectively)
- 4.2) optional: Run fastqc to see effects of ribodetector (fastqc_3.sh)
-4.1) optional: Run sortmerna to reduce the rRNA contamination (sortmerna_single.sh)
- 5.) Run Salmon quantification step(salmon_quant_single.sh or salmon_quant.sh for single and paired reads respectively)
- 6.) Run multiqc to summarize all steps in one HTML report (multiqc.sh)

#R codes:

- 1.) Run edgeR_func_QL_combined_chunks.R to do edgeR pairwise comparison with the salmon counts
- 2.) Run gene_intersects.R to get upset plots and see genes overlap among contrasts for Up and downregulated genes
- 3.) Run compile_topGO.R to summarize topGO results
- 4.) Run Join_files.R to obtain core and pan stress contrasts  genes and Log2 fold change and TPM z-normalized heatmaps for exploration. Also, this scripts subset files to use with vHRR 
- 5.) Run CO_NET2.R to obtain a co-expression network analysis, modules and the graph to be read in Cytoscape





