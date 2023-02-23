# Load EGAD and the data files 
library(EGAD);require(data.table);require(igraph);require(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#https://bioc.ism.ac.jp/packages/3.12/bioc/vignettes/EGAD/inst/doc/EGAD.html
data_dir <- "D:/TESIS_PHD/CHAPTER3/DEG"
core_genes <- paste0(data_dir,"/","atlas_al_pan.tsv")
core_files <- fread(core_genes,header = T)
row.names(core_files) <- core_files$`#`;
genes <- core_files$`#`  
core_files$`#`  <- NULL
#row.names(core_files) <- colnames(subset_core)
core_files <- t(core_files); colnames(core_files) <- genes

datExpr0 <- as.data.frame(core_files)
# datExpr0 <- datExpr0[1:12,]
# datExpr0 <- as.data.frame(log2(datExpr0+1))
 

#datExpr0 <- cor(datExpr0)
#write.table(datExpr02,"D:/TESIS_PHD/CHAPTER3/jenn_Samplescor.tsv",)
# dd <- data.frame(names = genes,sum=colSums(datExpr0))
# dd$test <- dd$sum>quantile(dd$sum,0.25)
# dd <- dd[which(dd$test==TRUE),]
# datExpr0 <- datExpr0[,dd$names]
# core_files <- core_files[rowSums(core_files)>0]
# core_files2 <- core_files
# core_files2 <- log2(core_files2)

 gsg = goodSamplesGenes(datExpr0, verbose = 3);
 gsg$allOK
 
 if (!gsg$allOK){
   # Optionally, print the gene and sample names that were removed:
   if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
   if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
   # Remove the offending genes and samples from the data:
   datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
 }
 par(mfrow = c(1,1));
 sampleTree = hclust(dist(datExpr0), method = "average");
 plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
      cex.axis = 1.5, cex.main = 2)
 
 
 
 enableWGCNAThreads(nThreads = 4)
 
 
 # Choose a set of soft-thresholding powers
 powers = c(c(1:10), seq(from = 12, to=120, by=2))
 # Call the network topology analysis function
 sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
 
 
 par(mfrow = c(1,2));
 cex1 = 0.9;
 # Scale-free topology fit index as a function of the soft-thresholding power
 plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
      main = paste("Scale independence"));
 text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      labels=powers,cex=cex1,col="red");
 # this line corresponds to using an R^2 cut-off of h
 abline(h=0.80,col="red")
 # Mean connectivity as a function of the soft-thresholding power
 plot(sft$fitIndices[,1], sft$fitIndices[,5],
      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
      main = paste("Mean connectivity"))
 text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

 
 collectGarbage();
 

 
 stress_net = blockwiseModules(datExpr0, power = 30, maxBlockSize = 15000,
                               TOMType = "unsigned", minModuleSize = 30,
                               reassignThreshold = 0, mergeCutHeight = 0.25,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "stressTOM",
                               verbose = 3)
 disableWGCNAThreads()
gc()
 stress_mergedColors = labels2colors(stress_net$colors)
 
 # Plot the dendrogram and the module colors underneath
 # plotDendroAndColors(control_net$dendrograms[[1]], control_mergedColors[control_net$blockGenes[[1]]],
 #                     "Module colors",
 #                     dendroLabels = FALSE, hang = 0.03,
 #                     addGuide = TRUE, guideHang = 0.05)
 
 plotDendroAndColors(stress_net$dendrograms[[1]], stress_mergedColors[stress_net$blockGenes[[1]]],
                     "Module colors",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05)
 
 
 stress_moduleLabels = stress_net$colors
 stress_moduleColors = labels2colors(stress_net$colors)
 stress_MEs = stress_net$MEs;
 stress_geneTree = stress_net$dendrograms[[1]];
 
 
 # Define numbers of genes and samples
 nGenes = ncol(datExpr0);
 nSamples = nrow(datExpr0);
 # Recalculate MEs with color labels
 # MEs0 = moduleEigengenes(datExpr, control_moduleColors)$eigengenes
 MEs0 = moduleEigengenes(datExpr0, stress_moduleColors)$eigengenes
 MEs = orderMEs(MEs0)

 disableWGCNAThreads()
 gc()
 
 
 # Calculate eigengenes
 MEList = moduleEigengenes(datExpr0, colors = stress_moduleColors)
 MEs = MEList$eigengenes
 # Calculate dissimilarity of module eigengenes
 MEDiss = 1-cor(MEs);
 
 # Cluster module eigengenes
 METree = hclust(as.dist(MEDiss), method = "average");
 # Plot the result
 sizeGrWindow(7, 6)
 plot(METree, main = "Clustering of module eigengenes",
      xlab = "", sub = "")
 
 MEDissThres = 0.25
 # Plot the cut line into the dendrogram
 abline(h=MEDissThres, col = "red")
 # Call an automatic merging function
 merge = mergeCloseModules(datExpr0, stress_mergedColors, cutHeight = MEDissThres, verbose = 3)
 # The merged module colors
 mergedColors = merge$colors;
 # Eigengenes of the new merged modules:
 mergedMEs = merge$newMEs;
 # Recalculate topological overlap if needed
 TOM = TOMsimilarityFromExpr(datExpr0, power = 30)
 
 
 # Select modules
 modules = unique(stress_moduleColors)#c("brown", "red");
 # Select module probes
 probes = names(datExpr0)
 inModule = is.finite(match(modules, modules));
 modProbes = probes[inModule];
 #modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
 # Select the corresponding Topological Overlap
 modTOM = TOM[inModule, inModule];
 
 
 setwd(data_dir)
 cyt = exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.02,
                                nodeNames = modProbes,
                                altNodeNames = modProbes,
                                nodeAttr = stress_moduleColors);
 

 
 