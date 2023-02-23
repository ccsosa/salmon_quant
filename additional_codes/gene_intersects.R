################################################################################
require(data.table);require(UpSetR)
################################################################################
#parameters
data_dir <- "D:/TESIS_PHD/CHAPTER3/DEG"
studies <- c("COMPLETE")
pval <- 0.05
################################################################################
#getting 
contrasts_i <- lapply(seq_len(length(studies)),function(i){
  message(paste("Getting results for ",studies[[i]]))
  stu_dir <- paste0(data_dir,"/",studies[[i]])
  csv_stu_dir <- paste0(stu_dir,"/","csv")
  contrasts <- fread(paste0(stu_dir,"/","contrasts.csv"))
  contrasts <- colnames(contrasts)
  contrasts <- contrasts[-1]
  
  contrasts_j <- lapply(seq_len(length(contrasts)),function(j){
    x  <- fread(paste0(csv_stu_dir,"/","glmQLFTest_",contrasts[[j]],"_pval_",pval,"_filtered.csv"),
                header=T)
    if(nrow(x)>0){
      #message(paste("Getting results for ",contrasts[[j]]))
      x$study <- studies[[i]]
      x$contrasts <- contrasts[[j]]
      return(x)
    }
  })
  
  contrasts_j <- do.call(rbind,contrasts_j)
  return(contrasts_j)
})

#joining genes and contrasts
contrasts_i <- do.call(rbind,contrasts_i)
contrasts_i$comb <- paste0(contrasts_i$study,"/",contrasts_i$contrasts)
#defining up and downregulating gene lists
contrasts_i_up <- contrasts_i[which(contrasts_i$status_name=="UP"),]
contrasts_i_down <- contrasts_i[which(contrasts_i$status_name=="DOWN"),]
#getting unique genes and contrasts
up_un <- unique(contrasts_i_up$comb)
down_un <- unique(contrasts_i_down$comb)
genes_up <- unique(contrasts_i_up$V1)
genes_down <- unique(contrasts_i_down$V1)

################################################################################
#Obtaining a matrix with 1 and 0 for upregulated genes
genes_up_m <- as.data.frame(matrix(nrow = length(genes_up),
                                   ncol = length(up_un)+1))


colnames(genes_up_m) <- c("genes",up_un)
genes_up_m$genes <- genes_up

pb <-
  utils::txtProgressBar(min = 0,
                        max = nrow(genes_up_m),
                        style = 3)

for(i in 1:nrow(genes_up_m)){
  utils::setTxtProgressBar(pb, i)
  x_gen <- contrasts_i_up[which(contrasts_i_up$V1==genes_up_m$genes[[i]]),]
  genes_up_m[i,-1] <- colnames(genes_up_m)[-1] %in% x_gen$comb*1
  rm(x_gen)
};rm(i)
close(pb)

################################################################################
#Obtaining a matrix with 1 and 0 for downregulated genes
genes_down_m <- as.data.frame(matrix(nrow = length(genes_down),
                                   ncol = length(down_un)+1))


colnames(genes_down_m) <- c("genes",down_un)
genes_down_m$genes <- genes_down

pb <-
  utils::txtProgressBar(min = 0,
                        max = nrow(genes_down_m),
                        style = 3)

for(i in 1:nrow(genes_down_m)){
  utils::setTxtProgressBar(pb, i)
  x_gen <- contrasts_i_down[which(contrasts_i_down$V1==genes_down_m$genes[[i]]),]
  genes_down_m[i,-1] <- colnames(genes_down_m)[-1] %in% x_gen$comb*1
  rm(x_gen)
};rm(i)
close(pb)

################################################################################
#upset upregulated
x_upstep_up <- upset(genes_up_m,number.angles = 45,point.size = 3.5, line.size = 2,
                     mainbar.y.label = "Upregulated genes intersection", sets.x.label = "Genes Per study/contrast",
                     order.by = "freq", keep.order = TRUE,
                     nintersects=NA,
                     set_size.numbers_size=T,
                     text.scale=1.5,
                     #nsets = length(up_un),
                     #show.numbers	=T ,
                     #scale.intersections = "log10",
                      sets = c(
                        "COMPLETE/TAluminum.AZU.10D-Control.AZU.10D",
                        "COMPLETE/TAluminum.BGI.10D-Control.BGI.10D",
                        #"COMPLETE/TAluminum.AZU_I.4HR-Control.AZU_I.4HR",
                        "COMPLETE/TAluminum.AZU.4HR-Control.AZU.4HR"
                        
                               #  "COMPLETE/TAluminum.AZU_I-Control.AZU_I", 
                               # "COMPLETE/TAluminum.AZU-Control.AZU",
                               # "COMPLETE/TAluminum.BGI-Control.BGI"
                                            ))#, empty.intersections = "off", cutoff = 1)#, 
#                     text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))
#x_upstep_up


pdf(paste0(data_dir,"/","UPSET_UP.pdf"),height = 8, width = 30)#15
x_upstep_up
dev.off()


################################################################################
#upset downregulated

#genes_down_m$`COMPLETE/TAluminum.AZU_I-Control.AZU_I` <- 0
x_upstep_down <- upset(genes_down_m,number.angles = 45,point.size = 3.5, line.size = 2,
                     mainbar.y.label = "Downregulated genes intersection", sets.x.label = "Genes Per study/contrast",
                     order.by = "freq", keep.order = TRUE,
                     nintersects=NA,
                     set_size.numbers_size=T,
                     text.scale=1.5,
                     #nsets = length(up_un),
                     #show.numbers	=T ,
                     #scale.intersections = "log10",
                     sets = c(
                       "COMPLETE/TAluminum.AZU.10D-Control.AZU.10D",
                       "COMPLETE/TAluminum.BGI.10D-Control.BGI.10D",
                       #"COMPLETE/TAluminum.AZU_I.4HR-Control.AZU_I.4HR",
                       "COMPLETE/TAluminum.AZU.4HR-Control.AZU.4HR"
                       
                       #  "COMPLETE/TAluminum.AZU_I-Control.AZU_I", 
                       # "COMPLETE/TAluminum.AZU-Control.AZU",
                       # "COMPLETE/TAluminum.BGI-Control.BGI"
                     ))#, empty.intersections = "off", cutoff = 1)#, 
#pdf(paste0(data_dir,"/","UPSET_DOWN.pdf"),height = 8, width = 15)
pdf(paste0(data_dir,"/","UPSET_DOWN.pdf"),height = 8, width = 30)
x_upstep_down
dev.off()
################################################################################
#x_upstep_down
#
filtered_groups <-              c(
  
  
  "COMPLETE/TAluminum.AZU.10D-Control.AZU.10D",
  "COMPLETE/TAluminum.BGI.10D-Control.BGI.10D",
  #"COMPLETE/TAluminum.AZU_I.4HR-Control.AZU_I.4HR",
  "COMPLETE/TAluminum.AZU.4HR-Control.AZU.4HR"
                              # "COMPLETE/TAluminum.AZU_I-Control.AZU_I", 
                              # "COMPLETE/TAluminum.AZU-Control.AZU",
                              # "COMPLETE/TAluminum.BGI-Control.BGI"
                     )
require(vegan)
genes_up_m2 <- as.data.frame(t(genes_up_m[,filtered_groups]))
colnames(genes_up_m2) <- row.names(genes_up_m)
row.names(genes_up_m2) <- filtered_groups#colnames(genes_up_m[,-1])
################################################################################
genes_down_m2 <- as.data.frame(t(genes_down_m[,filtered_groups]))
colnames(genes_down_m2) <- row.names(genes_down_m)
row.names(genes_down_m2) <- filtered_groups# colnames(genes_down_m[,-1])

pdf(paste0(data_dir,"/","CLUSTER_SIMILARITIES.pdf"),height = 8, width = 15)

par(mfrow = c(1, 2)) # Create a 2 x 2 plotting matrix
plot(hclust(vegdist(genes_up_m2,binary = T),method = "ward.D"),main="Upregulated genes",xlab="",sub="")
plot(hclust(vegdist(genes_down_m2,binary = T),method = "ward.D"),main="Downregulated genes",xlab="",sub="")
dev.off()
################################################################################
up_summary <- data.frame(comb=colnames(genes_up_m[,-1]),counts=colSums(genes_up_m[,-1]),status="Upregulated")
down_summary <- data.frame(comb=colnames(genes_down_m[,-1]),counts=colSums(genes_down_m[,-1]),status="Downregulated")
summary_counts <- rbind(up_summary,down_summary)
write.table(summary_counts,paste0(data_dir,"/","Genes_counts_summary.tsv"),
            sep = "\t",na = "",
            quote = F,row.names = F)

