
#doing for each combination
topGO_func_plaza_file <- function(indir,folder,ta_go_file,custom_format,contrasts,pval){
  
  ################################################################################  
  #calling libraries  
  require(data.table);require(topGO);require(dplyr);
  require(GO.db);require(parallel)
  #require(forcats);require(ggrepel) # for manhattan plot
  ################################################################################
  #files to used
  file <- paste0(indir,"/",folder,"/csv/glmQLFTest_",contrasts,"_pval_0.05_filtered.csv")
  if(any(!file.exists(file))){
    stop("Fix contrasts or folder!")
  }
  #background
  file_back <- paste0(indir,"/",folder,"/csv/glmQLFTest_",contrasts,"_pval_0.05_fulltable.csv")
  ################################################################################
  message("Reading GOA file")
  #reading goa file downloaded from plaza
  ta_go_file <- fread(ta_go_file)
  
  if(custom_format==F){
    #changing "#gene_id" colname to "gene_id"
    colnames(ta_go_file)[1] <- unlist(strsplit(colnames(ta_go_file)[1],"#"))[[2]]
    #subsetting to reduce size
    ta_go_file <- ta_go_file[,.(gene_id,go,description)]
  } else {
    colnames(ta_go_file)<- c("go","gene_id","na")
    #subsetting to reduce size
    ta_go_file <- ta_go_file[,.(gene_id,go)]  
    
  }
  
  #obtain children GO terms
   # GOBPCHILDREN_list <- as.list(GOBPCHILDREN)
   # GOBPCHILDREN_list <- GOBPCHILDREN_list[!is.na(GOBPCHILDREN_list)]
  # GOOBSOLETE_list <-  as.list(GOOBSOLETE)
  # GOOBSOLETE_list <- GOOBSOLETE_list[!is.na(GOOBSOLETE_list)]
  
  ################################################################################
  
  message("######################################")
  message("Running functional enrichment analyses")
  message("######################################")
  
  ################################################################################
  
  x <- lapply(seq_len(length(file)),function(i){
    #i <- 1
    message("   ")
    message(paste("Running", contrasts[[i]]),"\n")
    
    #reading complete analyzed genes for background
    x_back <- read.csv(file_back[[i]],header = T)
    x_back <- x_back[,1]
    
    
    #indexing background genes to subset and do not use the complete genome
    index_b <- ta_go_file$gene_id %in% x_back
    x_back_vect <- ta_go_file[index_b,]
    
    #creating a list of gene ids and their GO terms     
    gene_2_GO <- unstack(x_back_vect[,c(2,1)])
    
    #G2g <- inverseList(gene_2_GO) #convert from gene2GO to GO2gene
    #go2genes <- annFUN.gene2GO(whichOnto = "BP", gene2GO = gene_2_GO)
    
    ###############################
    ###############################
    ###############################
    
    # reading DEG files
    x <- read.csv(file[[i]],header = T)
    
    #splitting file in up and downregulated genes
    
    
    x_up <- x[which(x$status_name=="UP"),][,1]
    x_down <- x[which(x$status_name=="DOWN"),][,1]
    
    
    #joining in a list up and down to get results in a list
    genes_stat_list <- list(x_up,x_down)
    
    
    # running fisher tests by up and down respectively
    Fisher_results <- lapply(1:2,function(j){
      ###############################
      ###############################
      ###############################
      #j <- 1
      message("Running down and up regulated genes")
      message("Formatting to topGO format")
      #remover genes sin anotacion
      keep <- genes_stat_list[[j]] %in% x_back
      keep <- which(keep==TRUE)
      candidate_list <- genes_stat_list[[j]][keep]
      #converting gene list to factors 
      geneList <- factor(as.integer(x_back %in% candidate_list),levels = c(0,1))
      names(geneList) <- x_back
      
      ###############################
      ###############################
      ###############################
      message("Running topGO")
      #Create a topGO object
      GOdata <- new('topGOdata',
                    ontology='BP', 
                    allGenes = geneList, 
                    #annot = annFUN.GO2genes,#
                    annot = topGO::annFUN.gene2GO,
                    gene2GO = gene_2_GO)
      #new (get DAG level for each GO term)
      x_level <- buildLevels(GOdata@graph)
      x_level_l <- data.frame(GO=names(as.list(x_level$nodes2level)),
                              LEVEL=as.numeric(as.list(x_level$nodes2level)))
      

      #Run statistical tests
      resultFisher <- topGO::runTest(GOdata, 
                                     algorithm = "classic",
                                     statistic = "fisher")
      
      #returning all GO terms used
      allGO1 <- topGO::usedGO(GOdata)
      allGO1a <-allGO1[allGO1 %in% ta_go_file$go]
      
      ###############################
      ###############################
      ###############################
      message("Summarizing  topGO results")
      
      #summarizing results in one table
      all_res1 <- topGO::GenTable(GOdata, 
                                  classicFisher=resultFisher,
                                  orderBy='classicFisher',
                                  topNodes=length(allGO1))
      
      
      #getting scores to avoid mistakes in GenTable
      pValue.classic <- topGO::score(resultFisher)
      #creating a data.frame to join with all_res1
      pValue.classic <- data.frame(GO.ID=names(pValue.classic),p_val=pValue.classic)
      #joining accurate p values
      all_res1 <- dplyr::left_join(all_res1,pValue.classic,c("GO.ID"))
      #fitting to goa file
      ###############################
      message("Pre-filtering step: keeping only GO Terms used in the GOA input file")
      all_res1 <-all_res1[all_res1$GO.ID %in% allGO1a,]
      ###############################
      #adjust p value using fisher test
      all_res1$FDR <- stats::p.adjust(all_res1$p_val,method = "BH")
      #getting -log FDR for a Manhattan plot
      #all_res1$logpv <- -log10(all_res1$FDR)
      #creating categories significant and no significant
      all_res1$status_name <- NA
      all_res1$status_name[which(all_res1$FDR<pval)] <- "significant"
      all_res1$status_name[which(all_res1$FDR>=pval)] <- "no significant"
      #removing classicFisher column (this is replaced by p_val )
      all_res1$classicFisher <- NULL
      
      ###############################
      ###############################
      ###############################
      message("Subsetting only significant results")
      #subsetting only significant
      all_res1_s <- all_res1[which(all_res1$status_name=="significant"),]
      
      #Only apply filtering steps if there are enriched results!
      if(nrow(all_res1_s)>0){
      all_res1_s$DAG_level <- NA
      
      ###############################
      ###checking to filter GO:BP where children have better FDR values
      GO_IDs_initial <- data.frame(INITIAL = all_res1_s$GO.ID,
                                   SUGGESTED = NA)
      #CHILDREN =NA)
      
      ###############################
      message("First filtering step:Filtering results by children terms, evaluating each term.
              Please be patient")
      
      pb <-
        utils::txtProgressBar(min = 0,
                              max = nrow(GO_IDs_initial),
                              style = 3)
      
      SUGGESTED <- list()
      for(m in seq_len(nrow(GO_IDs_initial))){
        
        #message(m)
        utils::setTxtProgressBar(pb, m)
        
        #getting children terms per GO:ID enriched
        children_ids <- stack(as.list(GOBPCHILDREN[GO_IDs_initial$INITIAL[[m]]]))
        #children_ids <-  stack(GOBPCHILDREN_list[GO_IDs_initial$INITIAL[[m]]])
        
        #children terms FDR
        x_fdr_i <- all_res1[all_res1$GO.ID %in% children_ids$values,c("GO.ID","FDR")]
        #parental term FDR
        par_fdr <- data.frame(GO.ID = GO_IDs_initial$INITIAL[[m]],
                              FDR = all_res1_s$FDR[which(all_res1$GO.ID==GO_IDs_initial$INITIAL[[m]])])
        
        x_fdr_i <- rbind(x_fdr_i,par_fdr)
        #GO_IDs_initial$SUGGESTED[[m]] <- x_fdr_i$GO.ID[which.min(x_fdr_i$FDR)]
        SUGGESTED[[m]] <- x_fdr_i$GO.ID[which.min(x_fdr_i$FDR)]
        
        rm(children_ids,x_fdr_i,par_fdr)
      };rm(m)
      
      close(pb)
      
      GO_IDs_initial$SUGGESTED <- unlist(SUGGESTED)
      rm(SUGGESTED)
      
      #Getting unique filtered GO terms
      SUGGESTED_GO <- unique(GO_IDs_initial$SUGGESTED)
      #Filtering for all the output and obtain a filter table
      all_res1_s <- all_res1_s[all_res1_s$GO.ID %in% SUGGESTED_GO,]
      
      
      for(l in seq_len(nrow(all_res1_s))){
        all_res1_s$DAG_level[[l]] <- x_level_l[which(x_level_l$GO==all_res1_s$GO.ID[[l]]),2]
        all_res1_s
        };rm(l)
      
      #sorting significant by FDR
      all_res1_s <- all_res1_s[order(all_res1_s$FDR,decreasing = F),]
      
      } else {
        message("No significant enriched results!")
        all_res1_s <- all_res1_s
      }
      ###############################
      ###############################
      ###############################
      message("Returning results")
      #obtaining full table and filtered by pvalue in a list with two slots
      #full table and filtered
      FS <- list(full_table = all_res1, filtered = all_res1_s)
      return(FS)
    })
    message("   ")
    return(Fisher_results)

  })
  names(x) <- contrasts
  
  message("DONE!")
  message("#####")
  message("     ")
  return(x)
}



#x <- topGO_func_plaza_file(indir,folder,ta_go_file,custom_format=T,contrasts,pval)
################################################################################
################################################################################
################################################################################
#defining directories and files
#you can read the file directly from plaza in this case I am using a local file
#"https://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04_5/GO/go.tae.csv.gz"

# pval for cutoff
pval <- 0.05
#file where the results are available
indir <- "D:/TESIS_PHD/CHAPTER3/DEG"
#name of the folder where the DEG results are available

#contrasts provided for enrichment
################################################################################

param_list <- list(folder=c("COMPLETE"),
                   contrasts = list(c(
                     "TAluminum.AZU.10D-Control.AZU.10D",
                     "TAluminum.BGI.10D-Control.BGI.10D",
                     #"TAluminum.AZU_I.4HR-Control.AZU_I.4HR",
                     "TAluminum.AZU.4HR-Control.AZU.4HR"
                                      # "TAluminum.AZU_I-Control.AZU_I", 
                                      # "TAluminum.AZU-Control.AZU",
                                      # "TAluminum.BGI-Control.BGI"
                   )))

################################################################################
################################################################################
#run all in one line

summ_all_stud_list <- list()
  for(a in 1:length(param_list$folder)){
    
    
    #a <- 1
    #running the code
    out_dir <- "D:/TESIS_PHD/CHAPTER3/GO_DEG"
    ou_dir1 <- paste0(out_dir,"/",param_list$folder[[a]])#folder)
    if(!dir.exists(ou_dir1)){dir.create(ou_dir1)}
    
    ou_dir2 <- paste0(ou_dir1,"/",param_list$contrasts[[a]])#contrasts)
    for(i in 1:length(ou_dir2)){
      if(!dir.exists(ou_dir2[[i]])){dir.create(ou_dir2[[i]])}
    };rm(i)
    ################################################################################
    #run with plaza file
    ta_go_file <- "D:/TESIS_PHD/CHAPTER3/GO_PLAZA/go.osa.csv.gz" # plaza no filtered
    x <- topGO_func_plaza_file(indir,
                               folder=param_list$folder[[a]],
                               ta_go_file,
                               custom_format=F,
                               contrasts = param_list$contrasts[[a]],
                               pval)
    #run with Indee file
  
    
        # 
        # ta_go_file <- "/scratch/bis_klpoe/chsos/data/TA_GO_PLAZA/preprocessed_go_tae.tsv" # plaza filtered
        # x <- topGO_func_plaza_file(indir,
        #                        folder=param_list$folder[[a]],
        #                        ta_go_file,
        #                        custom_format=T,
        #                        contrasts = param_list$contrasts[[a]],
        #                        pval)
    #saving results
    
    summ_list <- list()
    for(i in 1:length(param_list$contrasts[[a]])){     #contrasts)){
      x2_up <-   x[[i]][[1]][[2]] #up filtered
      x2_down <- x[[i]][[2]][[2]] #down filtered
      x2_up_total <-   x[[i]][[1]][[1]] #up total
      x2_down_total <- x[[i]][[2]][[1]] #down total
      #top ten per filtered
      x2_up_top <- x2_up[1:10,]
      x2_down_top <- x2_down[1:10,]
      
      #saving tables filtered GO terms
      write.table(x2_up,paste0(ou_dir2[[i]],"/","topGO_UP_pval",pval,".tsv"),
                  quote=F,
                  na="",
                  row.names=F,
                  sep="\t")
      
      write.table(x2_down,paste0(ou_dir2[[i]],"/","topGO_DOWN_pval",pval,".tsv"),
                  quote=F,
                  na="",
                  row.names=F,
                  sep="\t")
      #saving total table
      
      write.table(x2_up_total,paste0(ou_dir2[[i]],"/","topGO_UP_total_pval",pval,".tsv"),
                  quote=F,
                  na="",
                  row.names=F,
                  sep="\t")
      
      write.table(x2_down_total,paste0(ou_dir2[[i]],"/","topGO_DOWN_total_pval",pval,".tsv"),
                  quote=F,
                  na="",
                  row.names=F,
                  sep="\t")
      
      #saving tables top ten of GO terms
      
      write.table(x2_up_top,paste0(ou_dir2[[i]],"/","topGO_UP_topten_pval_",pval,".tsv"),
                  quote=F,
                  na="",
                  row.names=F,
                  sep="\t")
      
      write.table(x2_down_top,paste0(ou_dir2[[i]],"/","topGO_DOWN_topten_pval_",pval,".tsv"),
                  quote=F,
                  na="",
                  row.names=F,
                  sep="\t")
      ########################################################################
      x_sum <- as.data.frame(matrix(ncol=5,nrow = 2))
      colnames(x_sum) <- c("status","total","filtered","study","contrast")
      x_sum[,1] <- c("DOWN","UP")
      x_sum[1,2] <- nrow(x2_down_total)
      x_sum[2,2] <- nrow(x2_up_total)
      x_sum[1,3] <- nrow(x2_down)
      x_sum[2,3] <- nrow(x2_up)
      x_sum$study <- param_list$folder[[a]]
      x_sum$contrast <- param_list$contrasts[[a]][[i]]
      
      summ_list[[i]]  <- x_sum
      
      rm(x2_up,x2_down,x2_down_total,x2_up_total,x_sum,x2_up_top,x2_down_top)
    };rm(i)
    
    #saving summary
    summ_list <- do.call(rbind,summ_list)
    summ_all_stud_list[[a]] <- summ_list
    rm(x)
    ################################################################################  
  };rm(a)

#obtaining total summary
summ_all_stud_list <- do.call(rbind,summ_all_stud_list)
#saving total summary
write.table(summ_all_stud_list,paste0(out_dir,"/","topGO_runs_summary_new",pval,".tsv"),
            quote=F,
            na="",
            row.names=F,
            sep="\t")
