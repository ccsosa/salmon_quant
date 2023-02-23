################################################################################
library(data.table);library(ggVennDiagram);library(ggplot2);require(ggpubr);
require(UpSetR);require(forcats);require(pheatmap);require(gtools);
require(DBI);require(GO.db);require(ggheatmap);require(AnnotationDbi)
#require(ComplexHeatmap)
################################################################################
#calculating Jaccard
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
################################################################################
data_dir <- "D:/TESIS_PHD/CHAPTER3/GO_DEG"
#data_dir <- "D:/TESIS_PHD/INTERNSHIP/DROUGHT/DEG/GO_DEG"
################################################################################
#Listing studies
studies <- list.dirs(data_dir,full.names = F,recursive = F)
#studies <- studies[-1]
stud_dirs <- paste0(data_dir,"/",studies)
################################################################################
#compiling all studies in one dataset
x_i <- lapply(seq_len(length(studies)),function(i){
#  i <-1
  sub_dirs <- list.dirs(stud_dirs[[i]],full.names = F,recursive = F)
  
  x_j <- lapply(seq_len(length(sub_dirs)),function(j){
   # j <- 1
    print(j)
    x_sub_j1 <- list.files(paste0(stud_dirs[[i]],"/",sub_dirs[[j]]),pattern = "topGO_UP_pval") 
    x_sub_j2 <- list.files(paste0(stud_dirs[[i]],"/",sub_dirs[[j]]),pattern = "topGO_DOWN_pval") 
    x_sub_j <- c(x_sub_j1,x_sub_j2);rm(x_sub_j1,x_sub_j2)
    x_k <- lapply(seq_len(length(x_sub_j)),function(k){
      x <- data.table::fread(paste0(stud_dirs[[i]],"/",
                                    sub_dirs[[j]],"/",x_sub_j[[k]]))
      #creating extra columns for upset plot
      x$study <- studies[[i]]
      x$comparison <- sub_dirs[[j]]
      x$genelist <- x_sub_j[[k]]
      if(nrow(x)>0){
      x <- x
      } else {
        x <- NULL
      }
      return(x)
    })
    x_k <- do.call(rbind,x_k)
    return(x_k)
  })
  x_j <- do.call(rbind,x_j)
  if(nrow(x_j)>0) {
    x_j <- x_j
  } else {
    x_j <- NULL
  }
  return(x_j)
})

x_i<-x_i[!sapply(x_i,is.null)]
x_i <- do.call(rbind,x_i)
x_i$comb <- paste0(x_i$study,"/",x_i$comparison)


#saving compiled file
write.table(x_i,paste0(data_dir,"/","GO_DEG_compile.tsv"),
            sep = "\t",na = "",
            quote = F,row.names = F)
################################################################################
#subsetting in up and down datasets
x_i_up <- x_i[which(x_i$genelist=="topGO_UP_pval0.05.tsv"),]
x_i_down <- x_i[which(x_i$genelist=="topGO_DOWN_pval0.05.tsv"),]
################################################################################
#calculating Jaccard coefficient
un_comb <- unique(x_i$comb)
comb <- as.data.frame(t(combn(un_comb,2)))
comb$JACCARD_UP <- NA
comb$JACCARD_DOWN <- NA
colnames(comb)[1:2] <- c("SOURCE","TARGET")

for(i in seq_len(nrow(comb))){
  comb$JACCARD_UP[[i]] <- jaccard(x_i_up$Term[which(x_i_up$comb==comb[i,1])],
                                  x_i_up$Term[which(x_i_up$comb==comb[i,2])])
  
  comb$JACCARD_DOWN[[i]] <- jaccard(x_i_down$Term[which(x_i_down$comb==comb[i,1])],
                                    x_i_down$Term[which(x_i_down$comb==comb[i,2])])
  
};rm(i)
################################################################################
#up_intersect (getting shared GO terms)
x_up_count <- tapply(x_i_up$Term,x_i_up$Term,length)
x_up_count <- data.frame(GO_TERM = names(x_up_count), COUNT= x_up_count)
rownames(x_up_count) <- NULL
x_up_count <- x_up_count[order(x_up_count$COUNT,decreasing = T),]
x_up_count <- x_up_count[which(x_up_count$COUNT>1),]
x_up_count$STATUS <- "Up"
################################################################################
#down_intersect (getting shared GO terms)
x_down_count <- tapply(x_i_down$Term,x_i_down$Term,length)
x_down_count <- data.frame(GO_TERM = names(x_down_count), COUNT= x_down_count)
rownames(x_down_count) <- NULL
x_down_count <- x_down_count[order(x_down_count$COUNT,decreasing = T),]
x_down_count <- x_down_count[which(x_down_count$COUNT>1),]
x_down_count$STATUS <- "Down"
################################################################################
#Obtaining up and down shared GO terms
int_up_down <- intersect(x_up_count$GO_TERM,x_down_count$GO_TERM)

x_count_summ <- rbind(x_up_count,x_down_count)
colnames(x_count_summ)
################################################################################
#saving intersected GO terms
write.table(x_count_summ,paste0(data_dir,"/","GO_DEG_INT_SUMMARY.tsv"),
            sep = "\t",na = "",
            quote = F,row.names = F)
################################################################################
################################################################################
#getting unique genes and contrasts for upset
up_un <- unique(x_i_up$comb)
down_un <- unique(x_i_down$comb)
genes_up <- unique(x_i_up$GO.ID)
genes_down <- unique(x_i_down$GO.ID)
################################################################################
#top ten of GO terms per study/contrast
x_i_up_top_l <- lapply(seq_len(length(up_un)),function(i){
  x_i_up_top <- x_i_up[which(x_i_up$comb==up_un[[i]]),]
  x_i_up_top <- x_i_up_top[order(x_i_up_top$FDR,decreasing = T),]
  x_i_up_top <- x_i_up_top[1:10,]
  x_i_up_top <- x_i_up_top[which(!is.na(x_i_up_top$GO.ID)),]
  if(nrow(x_i_up_top)>0){
    return(x_i_up_top)
  }
})
x_i_up_top_l <- do.call(rbind,x_i_up_top_l)


write.table(x_i_up_top_l,paste0(data_dir,"/","top_ten_UP.tsv"),
            sep = "\t",na = "",
            quote = F,row.names = F)

# x_i_up_top_l_graph <- ggplot(data=x_i_up_top_l, aes(x=fct_reorder(Term,FDR,.desc = T), y=-log10(FDR), fill=comb)) +
#   geom_bar(stat="identity")+
#   xlab("")+
#   ylab("-Log10(p-value)")+
#   ggtitle("Top ten enriched GO terms per study/contrast by FDR values (Upregulated genes)") +
#   theme(
#     legend.title=element_blank(),
#   #legend.position = "none",
#   axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ggsave(paste0(data_dir,"/","top_ten_UP.pdf"), x_i_up_top_l_graph,height = 8, width = 10)


##############################################################
x_i_down_top_l <- lapply(seq_len(length(down_un)),function(i){
  x_i_down_top <- x_i_down[which(x_i_down$comb==down_un[[i]]),]
  x_i_down_top <- x_i_down_top[order(x_i_down_top$FDR,decreasing = T),]
  x_i_down_top <- x_i_down_top[1:10,]
  x_i_down_top <- x_i_down_top[which(!is.na(x_i_down_top$GO.ID)),]
  if(nrow(x_i_down_top)>0){
    return(x_i_down_top)
  }
})
x_i_down_top_l <- do.call(rbind,x_i_down_top_l)

write.table(x_i_down_top_l,paste0(data_dir,"/","top_ten_DOWN.tsv"),
            sep = "\t",na = "",
            quote = F,row.names = F)


# x_i_down_top_l_graph <- ggplot(data=x_i_down_top_l, aes(x=fct_reorder(Term,logpv,.desc = T), y=logpv, fill=comb)) +
#   geom_bar(stat="identity")+
#   xlab("")+
#   ylab("-Log10(p-value)")+
#   ggtitle("Top ten enriched GO terms per study/contrast by FDR values (Downregulated genes)") +
#   theme(
#     legend.title=element_blank(),
#     #legend.position = "none",
#     axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ggsave(paste0(data_dir,"/","top_ten_DOWN.pdf"), x_i_down_top_l_graph,height = 8, width = 10)
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
  x_gen <- x_i_up[which(x_i_up$GO.ID==genes_up_m$genes[[i]]),]
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
  x_gen <- x_i_down[which(x_i_down$GO.ID==genes_down_m$genes[[i]]),]
  genes_down_m[i,-1] <- colnames(genes_down_m)[-1] %in% x_gen$comb*1
  rm(x_gen)
};rm(i)
close(pb)

################################################################################


#upset upregulated
x_upstep_up <- upset(genes_up_m,number.angles = 45,point.size = 3.5, line.size = 2,
                     mainbar.y.label = "Upregulated genes GO terms intersection", sets.x.label = "GO terms enriched Per study/contrast",
                     order.by = "freq", keep.order = TRUE,
                     nintersects=NA,
                     set_size.numbers_size=T,
                     text.scale=1.5,)#,
#nsets = length(up_un),
#show.numbers	=T ,
#scale.intersections = "log10",
# sets = c("SRP356530/LEAVES_D-LEAVES_C", 
#          "SRP098756/CROWN_D-CROWN_C",
#          "SRP098756/ROOTS_D-CROWN_C", #delete
#          "SRP098756/ROOTS_D-ROOTS_C"
#          #"SRP098756/ROOTS_D-LEAVES_D",
#"SRP257474/CK0H-ABA6H",
#"SRP257474/DT6H-ABA6H"

#))#, empty.intersections = "off", cutoff = 1)#, 
#                     text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))
#x_upstep_up


pdf(paste0(data_dir,"/","UPSET_UP.pdf"),height = 8, width = 12)
x_upstep_up
dev.off()

write.table(genes_up_m,paste0(data_dir,"/","GO_presence_up.tsv"),
            sep = "\t",na = "",
            quote = F,row.names = F)

write.table(genes_down_m,paste0(data_dir,"/","GO_presence_down.tsv"),
            sep = "\t",na = "",
            quote = F,row.names = F)

################################################################################
#upset downregulated
x_upstep_down <- upset(genes_down_m,number.angles = 45,point.size = 3.5, line.size = 2,
                       mainbar.y.label = "Downregulated genes GO terms intersection", 
                       sets.x.label = "GO terms enriched Per study/contrast",
                       order.by = "freq", keep.order = TRUE,
                       nintersects=NA,
                       set_size.numbers_size=T,
                       text.scale=1.5,)#,
#nsets = length(up_un),
#show.numbers	=T ,
#scale.intersections = "log10",
# sets = c("SRP356530/LEAVES_D-LEAVES_C", 
#          "SRP098756/CROWN_D-CROWN_C",
#          #"SRP098756/ROOTS_D-CROWN_C",#
#          "SRP098756/ROOTS_D-ROOTS_C"#,
#          #"SRP098756/ROOTS_D-LEAVES_D",
#          #"SRP257474/CK0H-ABA6H",
#          #"SRP257474/DT6H-ABA6H"
# ))
pdf(paste0(data_dir,"/","UPSET_DOWN.pdf"),height = 8, width = 12)
x_upstep_down
dev.off()
################################################################################
#loadin
x_i_top <- rbind(x_i_up_top_l,x_i_down_top_l)
################################################################################
#get GO terms groups and groups for upregulated
GO_terms_top <- unique(x_i_top$GO.ID)
groups_top <- unique(x_i_top$comb)
################################################################################
#count matrix up
matrix_top_go_up <- as.data.frame(matrix(ncol=length(groups_top),
                                         nrow = length(GO_terms_top)))

colnames(matrix_top_go_up) <- groups_top; row.names(matrix_top_go_up) <- GO_terms_top
group_cont <- strsplit(groups_top,"/")

group_cont_m_up <- data.frame(STUDY=unlist(lapply(group_cont, `[[`, 1)),
                              CONTRAST=unlist(lapply(group_cont, `[[`, 2))
)

#filling out table with raw data
for(i in seq_len(length(groups_top))){
  x_g_i <- fread(paste0(data_dir,"/",
                        group_cont_m_up$STUDY[[i]],
                        "/",
                        group_cont_m_up$CONTRAST[[i]],
                        "/",
                        "topGO_UP_total_pval0.05.tsv"
  ),header = T)
  
  x_g_i <- x_g_i[x_g_i$GO.ID %in% row.names(matrix_top_go_up),]
  # x_g_i <- x_g_i[gtools::mixedorder(rownames(x_g_i)),]
  
  for(j in seq_len(nrow(matrix_top_go_up))){
    ind_j <- which(x_g_i$GO.ID==row.names(matrix_top_go_up)[[j]])
    if(length(ind_j)>0){
      matrix_top_go_up[j,i] <- x_g_i$FDR[ind_j]
    } else {
      matrix_top_go_up[j,i] <-  NA
    }
  };rm(j)
};rm(i)
################################################################################
#count matrix down
matrix_top_go_down <- as.data.frame(matrix(ncol=length(groups_top),
                                           nrow = length(GO_terms_top)))

colnames(matrix_top_go_down) <- groups_top; row.names(matrix_top_go_down) <- GO_terms_top
group_cont <- strsplit(groups_top,"/")

group_cont_m_down <- data.frame(STUDY=unlist(lapply(group_cont, `[[`, 1)),
                                CONTRAST=unlist(lapply(group_cont, `[[`, 2))
)

#filling out table with raw data
for(i in seq_len(length(groups_top))){
  #message(i)
  x_g_i <- fread(paste0(data_dir,"/",
                        group_cont_m_down$STUDY[[i]],
                        "/",
                        group_cont_m_down$CONTRAST[[i]],
                        "/",
                        "topGO_DOWN_total_pval0.05.tsv"
  ),header = T)
  
  x_g_i <- x_g_i[x_g_i$GO.ID %in% row.names(matrix_top_go_down),]
  # x_g_i <- x_g_i[gtools::mixedorder(rownames(x_g_i)),]
  for(j in seq_len(nrow(matrix_top_go_down))){
    ind_j <- which(x_g_i$GO.ID==row.names(matrix_top_go_down)[[j]])
    if(length(ind_j)>0){
      matrix_top_go_down[j,i] <- x_g_i$FDR[ind_j]
    } else {
      matrix_top_go_down[j,i] <-  NA
    }
  };rm(j)
};rm(i)
################################################################################
colnames(matrix_top_go_down) <- paste0(colnames(matrix_top_go_down),"")#," (DOWN)")
colnames(matrix_top_go_up) <- paste0(colnames(matrix_top_go_up)," ")#," (UP)")
matrix_top_go_final <- cbind(matrix_top_go_down,matrix_top_go_up)

ann_col <- data.frame(genelist=c(rep("DOWNREGULATED",3),#4),
                                 rep("UPREGULATED",3)))#4)))
row.names(ann_col) <- colnames(matrix_top_go_final)


matrix_top_go_final_2 <- matrix_top_go_final
#log10
#matrix_top_go_final_2 <- -log10(matrix_top_go_final_2)
#round
matrix_top_go_final_2 <- round(matrix_top_go_final_2,3)
#pvalue
matrix_top_go_final_2[matrix_top_go_final_2 >0.05] <- " "
#log10
#matrix_top_go_final_2[matrix_top_go_final_2 < -log10(0.05)] <- " "
#remove NA
matrix_top_go_final_2[is.na(matrix_top_go_final_2)] <- " "



GO_names <- stack(AnnotationDbi::Term(unlist(row.names(matrix_top_go_final))))
index_GO_n <-is.na(GO_names$ind)*1:nrow(GO_names)
index_GO_n <- index_GO_n[which(index_GO_n!=0)]

for(i in seq_len(length(index_GO_n))){
  GO_names$values[index_GO_n[i]] <-  row.names(matrix_top_go_final)[index_GO_n[i]]
};rm(i)

matrix_top_go_final2 <- matrix_top_go_final
row.names(matrix_top_go_final2) <- unlist(GO_names$values)
matrix_top_go_final2[is.na(matrix_top_go_final2)] <- 1
################################################################################
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
################################################################################
#pheatmap option
#log10
#matrix_top_go_final3 <- -log10(matrix_top_go_final)
                               
#matrix_top_go_final[which(matrix_top_go_final>0.05),] <-NA
#log10
#p <- pheatmap::pheatmap(matrix_top_go_final3,cluster_rows = F,
#pvalue
p <- pheatmap::pheatmap(matrix_top_go_final,cluster_rows = F,
         fontsize = 18,
         annotation_legend = T,#F,
         labels_row = unlist(GO_names$values),
         annotation_col = ann_col,
         annotation_colors =  list(genelist=c(DOWNREGULATED="blue",#"#006EC9",
                                              UPREGULATED="red")),#"#FF00FF")),
         drop_levels=F,
         display_numbers = matrix_top_go_final_2,#T,#TRU1E,
         number_color = "white",
         fontsize_row = 20,#12
         fontsize_col =  16,#12
         border_color="gray98",
         fontsize_number=14,#10
         cellwidth = 80, #48
         angle_col = 45,
         cellheight = 20, #16
          #pvalue
           legend_breaks = c(0,0.2,0.4,0.6,0.8,0.98,
                             max(matrix_top_go_final,na.rm = T)),
           legend_labels = c("0","0.2","0.4","0.6","0.8","1","p-value\n"), 
          #log10
         #legend_breaks = c(0,15,30,45,60,75,90,93),
         #legend_labels = c("0","15","30","45","60","75","90","-log10(p-value)\n"),
        #color
         color = hcl.colors(8, "Greens"), #Oranges
         
         #color =rev(hcl.colors(6, "BluYl")),#hcl.colors(10, "BluYl"),
         legend = T,
         annotation_names_col = F,
         na_col = "gray",
         number_format= "%.2f",
         silent = F

          )

# ################################################################################
save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP.pdf"),
                  width = 36,
                  height = 36 )

# ################################################################################
# p <- pheatmap::pheatmap(t(matrix_top_go_final),
#                         #cluster_rows = F,
#                         fontsize = 18,
#                         annotation_legend = T,#F,
#                         labels_col = unlist(GO_names$values),
#                         annotation_row = ann_col,
#                         cluster_rows =T,
#                         cluster_cols=F,
#                         annotation_colors =  list(genelist=c(DOWNREGULATED="blue",#"#006EC9",
#                                                              UPREGULATED="red")),#"#FF00FF")),
#                         drop_levels=F,
#                         display_numbers = t(matrix_top_go_final_2),#T,#TRU1E,
#                         number_color = "white",
#                         fontsize_row = 28,#12
#                         fontsize_col =  20,#12
#                         border_color="gray98",
#                         fontsize_number=16,#10
#                         cellwidth = 44, #48
#                         angle_col = 270,
#                         cellheight = 30, #16
#                         #pvalue
#                         legend_breaks = c(0,0.2,0.4,0.6,0.8,0.98,
#                                           max(matrix_top_go_final,na.rm = T)),
#                         legend_labels = c("0","0.2","0.4","0.6","0.8","1","p-value\n"), 
#                         #log10
#                         #legend_breaks = c(0,15,30,45,60,75,90,93),
#                         #legend_labels = c("0","15","30","45","60","75","90","-log10(p-value)\n"),
#                         #color
#                         color = hcl.colors(8, "Greens"), #Oranges
#                         
#                         #color =rev(hcl.colors(6, "BluYl")),#hcl.colors(10, "BluYl"),
#                         legend = T,
#                         annotation_names_col = T,
#                         na_col = "gray",
#                         number_format= "%.2f",
#                         silent = F
#                         
# )
# save_pheatmap_pdf(x = p,
#                   filename =paste0(data_dir,"/","HEATMAP2.pdf"),
#                   width = 58,
#                   height = 15 )

# pdf(paste0(data_dir,"/","HEATMAP.pdf"),



#     height = 30,#24,
#     width = 24)#20)
# p
# dev.off()
################################################################################
# genecol <- c("red","green")
# names(genecol) <- c("DOWNREGULATED","UPREGULATED")
# col <- list(genelist=genecol)
# col_metaData  <- data.frame(
#   genelist=c(rep("DOWNREGULATED",4),
#              rep("UPREGULATED",4)))
# row.names(col_metaData) <- colnames(matrix_top_go_final2)
# p <-ggheatmap(matrix_top_go_final2,
#               cluster_cols = T,
#               cluster_rows = F,
#               scale = "none",
#               legendName = "p-value",
#               #annotation_position_cols="bottom",
#               annotation_cols = col_metaData,
#               annotation_color = col,
#               annotation_rows = NULL,
#               dist_method="euclidean",
#               hclust_method="ward.D",
#               border = "gray"
#               
# )#%>%
# #p
# # ggheatmap_plotlist(p)
# 
# p <-  ggheatmap_theme(p,1,
#                   theme =list(
#                     theme(axis.text.x = element_text(angle = 90,face = "bold",size = 8),
#                           axis.text.y = element_text(colour = "black",face = "bold")),
#                     theme(legend.title = element_text(face = "bold")),
#                     theme(legend.title = element_text(face = "bold")),
#                     theme(legend.title = element_text(face = "bold")),
#                     theme(legend.title = element_text(face = "bold"))
#                   )) 
#  
# #p  
# ggsave(paste0(data_dir,"/","HEATMAP.PDF"), p,height = 15, width = 10)

