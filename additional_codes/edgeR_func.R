
DEG_edgeR_func <- function(chunk,lfc, pval,out_dir,met_dir,levels,omit_samples){
  message(paste("chunk",chunk))
  library(ggplot2); library(tidyr); library(edgeR) ;
  library(tximport) ; library(readr); library(tidyr); 
  library(openxlsx);require(data.table);library(RColorBrewer);
  require(reshape2);require(ggpubr);library(ggvenn);#library(nVennR)

  #deifining Jaccard similarity function
  jaccard_similarity <- function(A, B) {
    intersection = length(intersect(A, B))
    union = length(A) + length(B) - intersection
    x <- data.frame(Adiff =length(unique(setdiff(A,B))),
                    Bdiff =length(unique(setdiff(B,A))),
                    intersection =intersection,
                    union=union,
                    jaccard=intersection/union)
    return (x)
  }
  
  
  data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",chunk)
  tx2gene <- read.table(paste0(data_dir,"/","salmon_tx2gene.tsv"))
  out_dir_1 <- paste0(out_dir,"/",chunk)
  if(!dir.exists(out_dir_1)){
    dir.create(paste0(out_dir,"/",chunk))
  }
  
  
  #tx2gene
  meta <- read.csv(paste0(met_dir,"/","run",chunk,".csv"), header = TRUE, sep=",")
  meta$trt <- sapply(strsplit(meta$sample,"_"),"[[",1)
  
  
  if(!is.null(omit_samples)){
    meta <- meta[!meta$sample %in% omit_samples,]
    #meta <- meta[which(meta$sample!=omit_samples),]
  }
  
  #meta$sample
  list1 = paste0(data_dir,"/", meta$sample,"/", "quant.sf")
  if(!is.null(omit_samples)){
    #list1 <-list1[which(list1!=paste0(data_dir,"/", omit_samples,"/", "quant.sf"))]
    list1 <- list1[list1!=paste0(data_dir,"/", omit_samples,"/", "quant.sf")]
    #list1 <-list1[which(list1!=paste0(data_dir,"/", omit_samples,"/", "quant.sf"))]
  }
  
  txi <- tximport(files=list1, type = "salmon", tx2gene = tx2gene)
  data <- as.data.frame(txi$counts)
  ab <-  as.data.frame(txi$abundance)
  df <- colSums(data)
  
  
  write.table(df, paste0(out_dir_1,"/","library_size.txt"), col.names = F, quote = F)
  colnames(data) <- meta$sample
  colnames(ab) <- meta$sample
  trt <- factor(meta$trt,levels = levels)
  
  # Save normalized (but not filtered) CPM
  y <- DGEList(counts = data,group = trt)
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y)
  write.table(cpm(y), paste0(out_dir_1,"/","CPM_normalized.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  
  #save tpm
  write.table(ab, paste0(out_dir_1,"/","abundance_drought_tpm_genelevel.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  
  # Filtering
  keep <- rowSums(cpm(data)>1)>=2 #it was 3 
  #keep <- filterByExpr(data)
  
  table(keep) # true false
  data_kept <- data[keep,] # 11515 genes were removed, 22132 were retainded
  #data_kept <- data[keep,]
  
  # Normalization and MDS plot
  d <- DGEList(data_kept,group =trt)
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  
  
  t <- plotMDS(d,plot = F)
  MDS_df <- data.frame(t$x, t$y)
  cols=brewer.pal(nrow(meta), "Set1")
  
  
  MDS<- ggplot(data = MDS_df) +
    geom_point(aes(x = t.x, y =t.y,  fill = meta$trt),size = 3, shape = 21)  +
    theme_bw() +
    aes(label = meta$sample)+ #allpoints) +
    # geom_text(hjust=0, vjust=0) +
    theme(legend.box="horizontal") +
    labs(fill = "Isolate") +
    xlab("Leading logFC dim1") +
    ylab("Leading logFC dim2") +
    scale_fill_manual(values = cols) +
    # scale_colour_brewer(palette="BrBG") +
    scale_colour_fermenter() +
    ggtitle("rice_drought") +
    theme(panel.grid =element_blank())
  #MDS
  
  ggsave(paste0(out_dir_1,"/","MDS.pdf"), MDS,height = 8, width = 8)
  
  # DE analysis
  #trt <- meta$trt
  trt <- trt
  #design <- model.matrix(~ trt + 0)
  #design <- model.matrix(~ trt ) #preferred
  #design <- model.matrix(~0 + trt, data = data)
  #rownames(design) <- meta$sample
  
  
  #design
  design <- model.matrix( ~0 + trt )
  colnames( design ) <- levels( trt )
  contrasts <- makeContrasts( DRO_CTRL = DROUGHT - CONTROL ,
                              # and so on for the other contrasts
                              levels=design)
  
  
  d <- estimateDisp(y = d, design = design,robust = T)
  pdf(paste0(out_dir_1,"/","plotBCV.pdf"),height = 8, width = 8)
  plotBCV(d)
  dev.off()
  
  #exact tests no glm (can be skipped!)
  # et <- exactTest(d)
  # t_et <- topTags(et)
  # View(t_et$table)
  
  
  # Fit NB model through replicates of a treatment
  fit <- glmQLFit(y = d, design = design,contrast=contrasts[,"DRO_CTRL"])
  
  pdf(paste0(out_dir_1,"/","plotQLDisp.pdf"),height = 8, width = 8)
  plotQLDisp(fit)
  dev.off()
  
  #testing glmQFTest
  qlf.2vs1 <- glmQLFTest(glmfit = fit,contrast = contrasts[,"DRO_CTRL"])#,contrast = c(0,1)) #teest
  qlf.2vs1$table$FDR <- p.adjust(qlf.2vs1$table$PValue,method = "BH")
  is.updown_QL <- decideTestsDGE(qlf.2vs1, p.value = pval, adjust.method = "BH")
  summary(is.updown_QL)
  qlf.2vs1$table$is.updown <- as.numeric(is.updown_QL)
  #tapply(qlf.2vs1$table$is.updown,qlf.2vs1$table$is.updown,length)
  #t_et <- topTags(qlf.2vs1,n = 20) #see dif exact test
  #summary(decideTests(qlf.2vs1, p.value=0.05))
  
  
  # write.csv(is.updown_QL, file= paste0(out_dir_1,"/","glmQLFTest","_pval_", as.character(pval),".csv"),
  #           row.names = T,quote = F)
  

  
  
  lfc_range <- paste0(round(max(qlf.2vs1$table$logFC[which(qlf.2vs1$table$is.updown==-1)]),2),
                      "-",
                      round(min(qlf.2vs1$table$logFC[which(qlf.2vs1$table$is.updown==1)]),2)
                      )

  qlf <- qlf.2vs1$table
  qlf$gene <- row.names(qlf)
  qlf <- qlf[,c(colnames(qlf)[7],colnames(qlf)[1:6])]
  
  write.csv(qlf, file= paste0(out_dir_1,"/","glmQLFTest","_pval_", as.character(pval),"_","fulltable.csv"),
            row.names = F,quote = F)
  
  
  pdf(paste0(out_dir_1,"/","plotMD_glmQLFTest.pdf"),height = 8, width = 8)
  plotMD(qlf.2vs1)
  abline(h=c(-1, 1), col="blue")
  dev.off()
  
  
  q1 <- data.frame(count = tapply(qlf.2vs1$table$is.updown,qlf.2vs1$table$is.updown,length))
  ###multiple lfc approach
  
  message("APPLYNG LFC threshold approach")
  lfc_list <- list()
  lfc_list_genes <- list()
  
  # L <- as.matrix(t(design))
  # row.names(L) <- colnames(fit)
  # L[2,][which(L[2,]==0)] <- -1
  # L[1,] <- 0
  #using glmTreat
  
  for(i in 1:length(lfc)){
    message(lfc[[i]])
  lfc_i <- lfc[[i]]
  LR <- glmTreat(glmfit = fit, contrast = contrasts[,"DRO_CTRL"],lfc = lfc_i)
  LR$table$FDR <- p.adjust(LR$table$PValue, "BH")
  is.updown <- decideTestsDGE(LR, p.value = pval, adjust.method = "BH")
  LR$table$is.updown <- as.numeric(is.updown)
  #tapply(LR$table$is.updown,LR$table$is.updown,length)
  
  pdf(paste0(out_dir_1,"/","plotMD_glmTreat_lfc_",as.character(lfc_i),".pdf"),height = 8, width = 8)
  plotMD(LR, status=is.updown, col=c("red","blue"), legend="topright")
  abline(h=c(-1, 1), col="blue")
  dev.off()
  
  
  # write.csv(is.updown, file= paste0(out_dir_1,"/","lfc_", as.character(lfc),"_pval_", as.character(pval),".csv"),
  #           row.names = F,quote = T)
  
  LR_t <- LR$table
  LR_t$gene <- row.names(LR)
  LR_t <- LR_t[,c(colnames(LR_t)[7],colnames(LR_t)[1:6])]
  
  write.csv(LR_t, file= paste0(out_dir_1,"/","lfc_", as.character(lfc_i),"_pval_", as.character(pval),"_","fulltable.csv"),
            row.names = F,quote = F)
  
  #creating a summary table for each approach 
  
  
  
  q2 <- data.frame(count = tapply(LR$table$is.updown,LR$table$is.updown,length))
  
  
  lfc_list[[i]] <- q2
  lfc_list_genes[[i]] <-  LR_t[which(LR_t$is.updown!=0),]
  rm(LR,is.updown,LR_t)
  };rm(i)
  
  #verifying counts of decideTestsDGE
  # q1a <- data.frame(count = tapply(qlf.2vs1$table$is.updown[which(qlf.2vs1$table$FDR<0.05)],
  #                                  qlf.2vs1$table$is.updown[which(qlf.2vs1$table$FDR<0.05)],
  #                                  length))
  # 
  # q2a <- data.frame(count = tapply(LR$table$is.updown[which(LR$table$FDR<0.05)],
  #                                  LR$table$is.updown[which(LR$table$FDR<0.05)],
  #                                  length))
  # 
  
  
  
  summary <- as.data.frame(matrix(ncol = 4,nrow = length(lfc)+1))
  colnames(summary) <- c("approach","-1","0","1")
  summary[,1] <- c("glmQLFTest",paste0("glmTreat_lfc",lfc))
  summary[1,1] <- paste0(summary[1,1],"s","(",lfc_range,")")
  
  
  ###Venn diagrams
  #sorting summary files to do a new one
  lfc_list[[length(lfc_list)+1]] <- q1 # <- append(lfc_list, q1)
  lfc_list <- lfc_list[c((length(lfc)+1),1:length(lfc))]
  
  approach <- summary$approach
  approach <- sub(pattern = "glmTreat_",replacement = "",approach)
  
  #IF approach length >1 calculate intersects
  
  if(length(lfc_list)>1){
  message("Plotting venn diagrams...")
  #sorting genes for Venn diagram
  lfc_list_genes[[length(lfc_list_genes)+1]] <-   qlf[which(qlf$is.updown!=0),] # <- append(lfc_list, q1)
  lfc_list_genes <- lfc_list_genes[c((length(lfc_list_genes)),1:(length(lfc_list_genes))-1)]
  
  
  message("UP diagram")
  #UP
  lfc_list_genes_up <- lapply(1:length(lfc_list_genes),function(i){
    x <- lfc_list_genes[[i]]$gene[which(lfc_list_genes[[i]]$is.updown==1)]
   return(x)
  })
  names(lfc_list_genes_up) <- approach
  lfc_list_genes_up <- lapply(1:length(lfc_list_genes_up),function(i){
    if(length(lfc_list_genes_up[[i]])==0){
      x <- NULL
    } else {
      x <- lfc_list_genes_up[[i]]
      #names(x) <- names(lfc_list_genes_up)[[i]]
    }
    
  })
  names(lfc_list_genes_up) <- approach
  lfc_list_genes_up <- lfc_list_genes_up[!unlist(lapply(lfc_list_genes_up, is.null))]
  lfc_list_genes_up2 <- lfc_list_genes_up[c(summary[1,1],"lfc0.25","lfc0.5")]
  
  #cols=brewer.pal(length(lfc_list_genes_up), "Set1")
  #x_up <- process_region_data(Venn(lfc_list_genes_up))
  #x_up$count[which(x_up$count==0)] <- NA
  
  
  
#  pdf(paste0(out_dir_1,"/","gvenn_up.pdf"),height = 8, width = 8)
  # myNV <- plotVenn(sets = lfc_list_genes_up,sNames = names(lfc_list_genes_up),nCycles = 20000,showPlot = F)
  # showSVG(myNV, opacity=0.2,showLegend = T,fontScale = 1.5,
  #         outFile = paste0(out_dir_1,"/","gvenn_up.svg"),borderWidth = 0,systemShow = F,labelRegions = T)
 # dev.off()
  
  gvenn <- ggvenn(lfc_list_genes_up2,stroke_size = 0,
                  set_name_size = 2.5,text_size = 3)
  
  
  ggsave(paste0(out_dir_1,"/","gvenn_up.pdf"), gvenn,height = 8, width = 10)
  

  ##counting summary tables UP
  x_up <- list()  
  for(j in 1:length(lfc_list_genes_up)){
    x_up_j <- list()
    for(k in 2:length(lfc_list_genes_up)){
      if(k>j){
      x <- jaccard_similarity(A = lfc_list_genes_up[[j]],B=lfc_list_genes_up[[k]])
      x$A <- names(lfc_list_genes_up)[[j]]
      x$B <- names(lfc_list_genes_up)[[k]]
      x_up_j[[k]] <- x
      }
    };rm(k)
    x_up[[j]] <- do.call(rbind,x_up_j)
    rm(x_up_j)
  };rm(j)
  x_up <- do.call(rbind,x_up)
  x_up <- x_up[c(6,7,1:5)]
  x_up$status <- "UP"
  
  message("DOWN diagram")
  #DOWN
  lfc_list_genes_down <- lapply(1:length(lfc_list_genes),function(i){
    x <- lfc_list_genes[[i]]$gene[which(lfc_list_genes[[i]]$is.updown==-1)]
    return(x)
  })
  names(lfc_list_genes_down) <- approach
  lfc_list_genes_down <- lapply(1:length(lfc_list_genes_down),function(i){
    if(length(lfc_list_genes_down[[i]])==0){
      x <- NULL
    } else {
      x <- lfc_list_genes_down[[i]]
      #names(x) <- names(lfc_list_genes_up)[[i]]
    }
    
  })
  names(lfc_list_genes_down) <- approach
  
  lfc_list_genes_down <- lfc_list_genes_down[!unlist(lapply(lfc_list_genes_down, is.null))]
  lfc_list_genes_down2 <- lfc_list_genes_down[c(summary[1,1],"lfc0.25","lfc0.5")]
  
  gvenn <- ggvenn(lfc_list_genes_down2,stroke_size = 0,
                  set_name_size = 2.5,text_size = 3)
  ggsave(paste0(out_dir_1,"/","gvenn_down.pdf"), gvenn,height = 8, width = 10)
  
  ##counting summary tables DOWN
  x_down <- list()  
  for(j in 1:length(lfc_list_genes_up)){
    x_down_j <- list()
    for(k in 2:length(lfc_list_genes_up)){
      if(k>j){
        x <- jaccard_similarity(A = lfc_list_genes_up[[j]],B=lfc_list_genes_up[[k]])
        x$A <- names(lfc_list_genes_up)[[j]]
        x$B <- names(lfc_list_genes_up)[[k]]
        x_down_j[[k]] <- x
      }
    };rm(k)
    x_down[[j]] <- do.call(rbind,x_down_j)
    rm(x_down_j)
  };rm(j)
  x_down <- do.call(rbind,x_down)
  x_down <- x_down[c(6,7,1:5)]
  x_down$status <- "DOWN"
  
  
  #all
  
  message("!=0 diagram (DEG)")
  lfc_list_genes_all <- lapply(1:length(lfc_list_genes),function(i){
    if(nrow(lfc_list_genes[[i]])==0){
      x <- NULL
    } else {
      x <-lfc_list_genes[[i]]$gene
      #names(x) <- names(lfc_list_genes_up)[[i]]
    }
    return(x)
  })
  
  names(lfc_list_genes_all) <- approach
  
  lfc_list_genes_all <- lfc_list_genes_all[!unlist(lapply(lfc_list_genes_all, is.null))] 
  lfc_list_genes_all2 <- lfc_list_genes_all[c(summary[1,1],"lfc0.25","lfc0.5")]
  
  gvenn <- ggvenn(lfc_list_genes_all2,stroke_size = 0,
                  set_name_size = 2.5,text_size = 3)
  ggsave(paste0(out_dir_1,"/","gvenn_all.pdf"), gvenn,height = 8, width = 10)
  
  ##counting summary tables DOWN
  x_all <- list()  
  for(j in 1:length(lfc_list_genes_up)){
    x_all_j <- list()
    for(k in 2:length(lfc_list_genes_up)){
      if(k>j){
        x <- jaccard_similarity(A = lfc_list_genes_up[[j]],B=lfc_list_genes_up[[k]])
        x$A <- names(lfc_list_genes_up)[[j]]
        x$B <- names(lfc_list_genes_up)[[k]]
        x_all_j[[k]] <- x
      }
    };rm(k)
    x_all[[j]] <- do.call(rbind,x_all_j)
    rm(x_all_j)
  };rm(j)
  x_all <- do.call(rbind,x_all)
  x_all <- x_all[c(6,7,1:5)]
  x_all$status <- "BOTH"
  
  } else {
    message("Ommiting Venn diagrams, only one approach available")
  }
  
  x_summary_counts <- rbind(x_up,x_down)
  x_summary_counts <- rbind(x_summary_counts,x_all)
  
  write.csv(x_summary_counts, file= paste0(out_dir_1,"/","gene_interesects_summary.csv"),
            row.names = F,quote = F)
  
  ###summary files
  for(i in 1:nrow(summary)){
    message(i)
    # if(i==1){
    #   x = q1
    # } else {
    #   x = lfc_list
    # }
    # 
    if(length(as.character(lfc_list[i][[1]][which(row.names(lfc_list[i][[1]])=="-1"),]))>0){
      summary[i,2] <- as.character(lfc_list[i][[1]][which(row.names(lfc_list[i][[1]])=="-1"),]) 
    } else {
      summary[i,2] <- NA
    }
    
    if(length(as.character(lfc_list[i][[1]][which(row.names(lfc_list[i][[1]])=="0"),]))>0){
      summary[i,3] <- as.character(lfc_list[i][[1]][which(row.names(lfc_list[i][[1]])=="0"),]) 
    } else {
      summary[i,3] <- NA
    }
    
    if(length(as.character(lfc_list[i][[1]][which(row.names(lfc_list[i][[1]])=="1"),]))>0){
      summary[i,4] <- as.character(lfc_list[i][[1]][which(row.names(lfc_list[i][[1]])=="1"),]) 
    } else {
      summary[i,4] <- NA
    }
  };rm(i)
  
  
  write.csv(summary, file= paste0(out_dir_1,"/","summary.csv"),
            row.names = F,quote = F)
  
  
  xx_sum <- melt(summary,"approach")
  xx_sum$value <- as.numeric(xx_sum$value)
  for(i in 1:nrow(xx_sum)){
    if(is.na(xx_sum$value[[i]])){
      xx_sum$value[[i]] <- NA
    } else {
      xx_sum$value[[i]] <- log10(xx_sum$value[[i]])
    }
  };rm(i)

  xx_sum_graph <- ggbarplot(xx_sum, x = "approach", y = "value", color = "variable",fill="variable",
            add = "mean_se", palette = c("red","gray","blue"),ylab="log10(counts)",
          #  angle=-90,
            position = position_dodge())
  xx_sum_graph <- xx_sum_graph + rotate_x_text(90)  
  
  ggsave(paste0(out_dir_1,"/","summary_graph.pdf"), xx_sum_graph,height = 8, width = 8)
  
  message("DONE!")
  
}


lfc=c(0.25,0.5,1,1.5,2)
chunk <- 2#2#2 #1
pval = 0.05
out_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"
met_dir <- "/scratch/bis_klpoe/chsos/data/sample_files/DONE"
levels <- factor(c("CONTROL","DROUGHT"))
omit_samples <- NULL# "CONTROL_SRR5225300" #$#NULL normally activate in 2
x <- DEG_edgeR_func(chunk,lfc, pval,out_dir,met_dir,levels,omit_samples)
#https://f1000research.com/articles/5-1438/v2


