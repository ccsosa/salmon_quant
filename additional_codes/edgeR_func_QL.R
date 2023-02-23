#' @title Automatic run of edgeR to obtain differential expressed genes
#'
#' @name DEG_edgeR_func
#' @description compareGDEG_edgeR_func function provides a simple workflow to obtain differential gene expressed from nf-core/rnaseq outputs
#' This function needs run first nf-core/rnaseq pipeline first. The steps done by this code are the following:
#' 1.) Create output dir folder structure with two subfolders: csv and graphics
#' 2.) Read Salmon outcomes (salmon_tx2gene.tsv, and tximport)
#' 3.) Read the metadata file
#' 4.) Create  groups for the samples based on the samples files provided for nf-core/rnaseq using combinatorics or using a group_vect object
#' 5.) Save raw counts obtained from salmon files
#' 6.) Creates a sample count table and calculate normalizing factors
#' 7.) Creates a exploratory plot (plotMDS) and save it (MDS.pdf)
#' 8.) Creates the design object to be use for differential expressed genes step
#' 9.) Estimate the dispersion  and save a plot of it (plotBCV.pdf)
#' 10.) Creates the contrast object to be use for differential expressed genes step (list)
#' 11.) Run in parallel the next steps per contrast:
#'      - Use a Fit a quasi-likelihood negative binomial generalized log-linear model (glmQLFit) with the design and contrasts
#'      - Use glmQLFTest test to get the DEG outcome table
#'      - Run a Benjamini-Hochberg procedure to obtain false discovery rate values (FDR)
#'      - Split up and downregulated genes using the sign of the log2fold change
#'      - Save raw and filtered by a p value results
#'      - Save plots of glmQLFit and biological coefficient variance
#' 12.) Creates a summary file for all up and downregulated genes for all the contrasts provided
#' 13.) Return a list with two slots: contrasts and the results of the contrasts (DEG)
#' 
#' @param folder_name Folder name with the results obtained by nf-core/rnaseq results. 
#'  If you have a folder name "2" this name will be used to call the metadata table and to name the folder with the outputs
#' @param data_dir Folder name with the results obtained by nf-core/rnaseq for the Salmon counts
#' @param out_dir Folder where the outcomes will be written
#' @param met_dir  Folder where the metadata file will be loaded. Please name your metadata as folder name (e.g. run2.csv")
#' @param pval p-value used to filter the results of differential expressed genes (default value = 0.05)
#' @param numCores numeric, Number of cores to use for the process (default value numCores=2)
#' @param plot_MDS This is a boolean value to indicate if the exploratory PCA plot should be displayed in the R session
#' @param group_vect (default value numCores=2) This is a character object which represents the groups available in the metadata data.
#'  If the value is NULL the script automatically will obtain groups using the sample names. For instance if there are four samples named
#'  CONTROL_SRR1, CONTROL_SRR2,DROUGHT_SRR3, DROUGHT_SRR4, the groups will be control and drought respectively. If values are provided, they
#'  must respect the order used in the metadata file. For the last example the group_vect object will be:
#'  group_vect <- c("CONTROL","CONTROL","DROUGHT","DROUGHT") 
#' @return This function will return a list with two slots: contrasts and the results of the contrasts (DEG) (contrasts and contrasts_results respectively)
#' @examples
#'
#' ####Loading parameters
#' #Directories
#' 
#' #Folder name
#' folder_name <- "2" #(If there are several runs this helps to only get results from chunks)
#' #Folder to save outputs
#' #Salmon results folder
#' data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
#' out_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"
#' #Folder where the nf-core/rnaseq sample metadata is available
#' met_dir <- "/scratch/bis_klpoe/chsos/data/sample_files/DONE"
#' #Other parameters
#' #p value for filtering 
#' pval = 0.05 
#' #Number of CPU cores to use in parallel
#' numCores <- 4 
#' #If plot should be appears in the R session
#' plot_MDS <- TRUE 
#' 
#' #Running function
#' 
#' x <- DEG_edgeR_func(data_dir = data_dir,
#'                     folder_name = folder_name,
#'                     out_dir = out_dir,
#'                     met_dir = met_dir,
#'                     pval = pval,
#'                     plot_MDS = TRUE,
#'                     numCores = 4,
#'                     )
#'
#' @importFrom utils combn setTxtProgressBar txtProgressBar p.adjust
#' @importFrom parallel makeCluster parLapplyLB stopCluster detectCores
#' @importFrom edgeR DGEList calcNormFactors cpm estimateDisp glmQLFit plotQLDisp
#' @importFrom limma plotMDS makeContrasts plotMD
#' @importFrom ggplot2 ggplot
#' @export


DEG_edgeR_func <- function(folder_name, data_dir,pval=0.05,data_dir,out_dir,met_dir,plot_MDS=T,numCores=2,group_vect=NULL){
  
  
  message(paste("Processing folder:",folder_name))
  library(ggplot2);library(edgeR);
  library(tximport);library(parallel)
  ##############################################################################
  # Detecting number of cores for parallelizing
  
  x_det <- NULL
  x_det <- parallel::detectCores()
  message(paste("Detecting cores, total cores are: ",x_det))
  ##############################################################################
  if (numCores > x_det) {
    stop("Number of cores exceed the maximum allowed by the machine,
         use a coherent number of cores such as four")
  }
  
  ##############################################################################
  #Creating output dir folder
  
  out_dir_1 <- paste0(out_dir,"/",folder_name)
  if(!dir.exists(out_dir_1)){
    dir.create(paste0(out_dir,"/",folder_name))
  }
  ##############################################################################
  # Defining plot dirs for contrasts
  plt_dir <- paste0(out_dir_1,"/graphics")
  if(!dir.exists(plt_dir)){
    dir.create(plt_dir)
  }
  
  # Defining csv dir for contrasts
  csv_dir <- paste0(out_dir_1,"/csv")
  if(!dir.exists(csv_dir)){
    dir.create(csv_dir)
  }
  
  ##############################################################################
  #Reading salmon files
  message("Reading salmon files")
  
  tx2gene <- read.table(paste0(data_dir,"/","salmon_tx2gene.tsv"))
  #reading metadata
  #tx2gene
  meta <- read.csv(paste0(met_dir,"/","run",folder_name,".csv"), header = TRUE, sep=",")
  meta$trt <- sapply(strsplit(meta$sample,"_"),"[[",1)
  
  
  #reading quant files
  list1 = paste0(data_dir,"/", meta$sample,"/", "quant.sf")
  
  #read salmon tximport
  txi <- tximport::tximport(files=list1, type = "salmon", tx2gene = tx2gene)
  data <- as.data.frame(txi$counts)
  ab <-  as.data.frame(txi$abundance)
  df <- colSums(data)
  
  
  write.table(df, paste0(out_dir_1,"/","library_size.txt"), col.names = F, quote = F)
  
  ##############################################################################
  #Testing that groups have more than one sample. (this avoid NA dispersion values)
  trt_summary <- as.data.frame(tapply(meta$trt,meta$trt,length))
  trt_summary$group <- row.names(trt_summary)
  trt_summary <- trt_summary[,c(2,1)]
  colnames(trt_summary) <- c("group","count")
  
  
  #Using group file if it is available
  if(!is.null(group_vect)){
    meta$trt <-   group_vect
  } else {
    #If a group have only one sample the function stops!
    if(any(trt_summary$count==1)){
      stop("Each group must have at least two samples, please check and provide
           an object group_vect with the group and run again.")
    }
  }
  
  
  
  ##############################################################################
  message("Defining groups with the metadata provided and levels")
  
  # Defining treatments
  colnames(data) <- meta$sample
  colnames(ab) <- meta$sample
  trt <- factor(meta$trt)#,levels = levels)
  
  #Defining possible combinations
  cmb_un <- as.data.frame(t(combn(levels(trt),2)))
  colnames(cmb_un) <- c("SOURCE","TARGET")
  cmb_un$COMP <- paste0(cmb_un$TARGET, "-", cmb_un$SOURCE)
  
  message("Possible combinations available:",nrow(cmb_un))
  print(cmb_un)
  # Save normalized (but not filtered) CPM
  y <- edgeR::DGEList(counts = data,group = trt)
  y$samples$lib.size <- colSums(y$counts)
  y <- edgeR::calcNormFactors(y)
  
  
  message("Saving raw data obtained (tpm and cpm)")
  
  #saving cpm
  write.table(edgeR::cpm(y), paste0(out_dir_1,"/","CPM_normalized.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  #save tpm
  write.table(ab, paste0(out_dir_1,"/","abundance_drought_tpm_genelevel.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  
  ##############################################################################
  # Filtering
  message("Filtering using cpm 2")
  keep <- rowSums(edgeR::cpm(data)>1)>=2 
  data_kept <- data[keep,]
  ##############################################################################
  #normal factors
  message("Normalizing")
  # Normalization and MDS plot
  d <- edgeR::DGEList(data_kept,group =trt)
  d$samples$lib.size <- colSums(d$counts)
  d <- edgeR::calcNormFactors(d)
  
  ##############################################################################
  #Plotting MDS
  t <- limma::plotMDS(d,plot = F)
  MDS_df <- data.frame(t$x, t$y)
  
  #saving MDS with shapes (it plots numbers)
  
  MDS<- ggplot2::ggplot(data = MDS_df) +
    geom_point(aes(x = t.x, y =t.y,shape = meta$trt),size = 4)  +
    theme_bw() +
    #    aes(label = meta$sample)+ #allpoints) +
    theme(legend.box="horizontal") +
    xlab("Leading logFC dim1") +
    ylab("Leading logFC dim2") +
    ggtitle("") +
    scale_shape_manual(values = c(1:length(trt), letters, LETTERS, "0", "1"))+
    theme(panel.grid =element_blank())
  
  MDS <- MDS +
    guides(shape = guide_legend(title = "Groups"))
  
  ggsave(paste0(out_dir_1,"/","MDS.pdf"), MDS,height = 8, width = 8)
  
  if(plot_MDS==TRUE){
    plot(MDS)  
  }
  ##############################################################################
  #Refining design to create contrasts
  #design
  trt <- trt
  #using treatments and 0 to do easily
  design <- model.matrix( ~0 + trt )
  colnames( design ) <- levels( trt )
  
  ##############################################################################
  message("Estimate dispersion")
  #estimateDisp
  d <- edgeR::estimateDisp(y = d, design = design,robust = T)
  pdf(paste0(out_dir_1,"/","plotBCV.pdf"),height = 8, width = 8)
  plotBCV(d)
  dev.off() 
  
  ##############################################################################
  #contrasts
  message("Creating contrasts")
  
  message("Using automatical combinations, analyze carefully!")
  if(nrow(cmb_un)>1){
    contrasts <- list()
    for(i in 1:nrow(cmb_un)){
      contrasts[[i]] <- limma::makeContrasts(contrasts =  cmb_un$COMP[[i]],
                                      levels=levels(trt))
    };rm(i)
  } else {
    contrasts <- list(limma::makeContrasts( contrasts = cmb_un$COMP[[1]] ,
                                     levels=levels(trt)))
  }
  
  message(paste("Contrasts provided: ",length(contrasts)))
  ##############################################################################
  #obtaining DEG per contrast 
  message(paste0("Calculating DE genes in parallel, using ",numCores," cores"))
  if(length(contrasts)==1){
    warning("Only one contrast will be run")
  }
  if(length(contrasts)>0){
    #Using parallel approach
    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("contrasts","d","design",
                                          "pval","cmb_un","plt_dir",
                                          "csv_dir"),envir=environment())
    
    cont_results <- parallel::parLapplyLB(cl,
                                          X = seq_len(length(contrasts)),
                                          fun = function (i){
                                            #glmQLFit
                                            fit <- edgeR::glmQLFit(y = d, design = design,contrast=contrasts[[i]]) 
                                            pdf(paste0(plt_dir,"/",cmb_un$COMP[[i]],"_BCV.pdf"),height = 8, width = 8)
                                            edgeR::plotQLDisp(fit)
                                            dev.off()
                                            
                                            #testing glmQFTest
                                            qlf.2vs1 <- edgeR::glmQLFTest(glmfit = fit,contrast = contrasts[[i]])
                                            qlf.2vs1$table$FDR <- p.adjust(qlf.2vs1$table$PValue,method = "BH")
                                            qlf.2vs1$table$status <- as.factor(sign(qlf.2vs1$table$logFC))
                                            #defining names for values 1 and -1 
                                            qlf.2vs1$table$status_name <- NA
                                            qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="1")] <- "UP"
                                            qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="-1")] <- "DOWN"
                                            #qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="1")] <- cmb_un$TARGET[[i]]
                                            #qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="-1")] <- cmb_un$SOURCE[[i]]
                                            #subsetting and saving
                                            qlf.2vs1_final <- qlf.2vs1$table
                                            qlf.2vs1_final2 <- qlf.2vs1_final[which(qlf.2vs1_final$FDR < pval),]
                                            
                                            message(paste0("Genes kept after filtering: ",nrow(qlf.2vs1_final2)), "printing groups...")
                                            print(tapply(qlf.2vs1_final2$status_name,qlf.2vs1_final2$status_name,length))
                                            
                                            #saving in csv tables (full and filtered)
                                            write.csv(qlf.2vs1_final, file= paste0(csv_dir,"/","glmQLFTest_",
                                                                                   cmb_un$COMP[[i]],"_","pval_",
                                                                                   as.character(pval),"_","fulltable.csv"),
                                                      row.names = T,quote = F)
                                            write.csv(qlf.2vs1_final2, file= paste0(csv_dir,"/","glmQLFTest_",
                                                                                    cmb_un$COMP[[i]],"_","pval_",
                                                                                    as.character(pval),"_","filtered.csv"),
                                                      row.names = T,quote = F) 
                                            
                                            
                                            #plot significant genes
                                            
                                            pdf(paste0(plt_dir,"/","plotMD_glmQLFTest_",cmb_un$COMP[[i]],".pdf"),height = 8, width = 8)
                                            limma::plotMD(qlf.2vs1)
                                            #abline(h=c(-1, 1), col="blue")
                                            dev.off()
                                            
                                            
                                            q1 <- data.frame(count = tapply(qlf.2vs1_final2$status_name,qlf.2vs1_final2$status_name,length))
                                            q1$groups <- row.names(q1)
                                            q1 <- q1[,c(2,1)]
                                            print(q1)
                                            #saving in csv tables (full and filtered)
                                            write.csv(q1, file= paste0(out_dir_1,"/","glmQLFTest_",
                                                                       cmb_un$COMP[[i]],"_","pval_",
                                                                       as.character(pval),"_","summary.csv"),
                                                      row.names = F)
                                            
                                            #Do not move it is need to return qlf.2vs1_final
                                            qlf.2vs1_final2
                                          })
    
    parallel::stopCluster(cl)
  } else {
    stop("no contrasts to run.") #NEW
    cont_results <- NULL # NEW
  }
  
  final_list <- list(contrasts = contrasts,
                     contrasts_results=cont_results)
  message("DONE!")
  return(final_list)
  
}

##############################################################################
#parameters
folder_name <- "2" #folder name
pval = 0.05 #p value for filtering
numCores <- 4 #number of cores to use in parallel
plot_MDS <- TRUE #if plot should be appears in R session
#groups if it not easy to get from the headers
#group_vect <- c("DT6H","CK0H","ABA6H","CK0H","DT6H","ABA6H","CK0H","DT6H","ABA6H")
group_vect <- NULL
##############################################################################
#Directories
#where are the salmon files
data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
#defining output folder
out_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"
#path where the nf core metadata is available
met_dir <- "/scratch/bis_klpoe/chsos/data/sample_files/DONE"
##############################################################################
x <- DEG_edgeR_func(folder_name=folder_name,
                    data_dir=data_dir,
                    pval=pval,
                    out_dir=out_dir,
                    met_dir=met_dir,
                    plot_MDS=TRUE,
                    numCores=4,
                    group_vect=group_vect)
