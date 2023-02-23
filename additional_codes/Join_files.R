require(dplyr);require(data.table);require(pheatmap);library(caret)
#Directories
#where are the salmon files
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
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
folder <- "D:/TESIS_PHD/CHAPTER3/DEG/COMPLETE"
folder_out <- "D:/TESIS_PHD/CHAPTER3/DEG"

vhRR_dir <- "D:/TESIS_PHD/CHAPTER3/vHHR/input"
x1 <- read.table(paste0(folder,"/","abundance_drought_tpm_genelevel.txt"),header = T)
x1$tx <- NULL
x1$gene_name <- NULL
print(paste(nrow(x1)," / ",1))
write.table(x1,paste0(vhRR_dir,"/atlas_al_total.tsv"),sep = "\t",na = "",row.names = T,col.names = F,quote = F)
write.table(x1,paste0(folder_out,"/atlas_al_total.tsv"),sep = "\t",na = "",row.names = T,col.names =T,quote = F)


samples <- read.table(paste0("D:/TESIS_PHD/CHAPTER3/","/SAMPLES.TSV"),header = T)

###########################
# require(data.table)
#contrasts
sets = c(
  # "COMPLETE/TAluminum.AZU_I-Control.AZU_I", 
  #          "COMPLETE/TAluminum.AZU-Control.AZU",
  #          "COMPLETE/TAluminum.BGI-Control.BGI"
  "COMPLETE/TAluminum.AZU.10D-Control.AZU.10D",
  "COMPLETE/TAluminum.BGI.10D-Control.BGI.10D",
  #"COMPLETE/TAluminum.AZU_I.4HR-Control.AZU_I.4HR",
  "COMPLETE/TAluminum.AZU.4HR-Control.AZU.4HR"
  )
pval <- 0.05
data_dir <- "D:/TESIS_PHD/CHAPTER3/DEG"        
sets_st_cont <- strsplit(sets,"/")

##loading all DEG for all contrasts
x_deg <- lapply(1:length(sets),function(i){
  message(i)
  x <- fread(paste0(data_dir,"/",sets_st_cont[[i]][1],"/csv/",
                    "glmQLFTest_",sets_st_cont[[i]][2],"_pval_",
                    pval,"_filtered.csv"))
  x$status <- NULL
  x$group <- sets[[i]]
  return(x)
})
x_deg <- do.call(rbind,x_deg)
#
#geeting all genes accross all studies
x_unique_genes <- as.data.frame(matrix(nrow=length(unique(x_deg$V1)),ncol=length(sets)+1))
x_unique_genes$V1 <- unique(x_deg$V1)
colnames(x_unique_genes) <- c("genes",sets)
row.names(x_unique_genes) <- unique(x_deg$V1)

pb <-
  utils::txtProgressBar(min = 0,
                        max = nrow(x_unique_genes),
                        style = 3)

#Presence of a gene in a set of contrasts
for(i in 1:nrow(x_unique_genes)){
  utils::setTxtProgressBar(pb, i)
  x_i <- x_deg[which(x_deg$V1==x_unique_genes$genes[[i]]),]
  x_unique_genes[i,-1] <- sets %in% x_i$group # transforming boolean in 1 and 0s
  x_unique_genes[i,-1] <- (x_unique_genes[i,-1]*1)
};rm(i)

close(pb)


#unique gene core 
x_unique_genes$total<- rowSums(x_unique_genes[,-1])
#subsetting
unique_genes <- unique(x_unique_genes$genes)
subset_pan <- x1[which(row.names(x1) %in% unique_genes),]
subset_pan <- subset_pan[rowSums(subset_pan[,]) > 0,]
write.table(subset_pan,paste0(vhRR_dir,"/atlas_al_pan.tsv"),sep = "\t",na = "",row.names = T,col.names = F,quote = F)
write.table(subset_pan,paste0(folder_out,"/atlas_al_pan.tsv"),sep = "\t",na = "",row.names = T,col.names = T,quote = F)

################################################################################
################################################################################
# process <- preProcess(subset_pan, method=c("range"))
# norm_scale <- predict(process, subset1)
#z-score
preproc1 <- preProcess(subset_pan, method=c("center", "scale"))
norm_scale <- predict(preproc1, subset_pan)
#quantile(stack(norm_scale)[,1])
# z-score
# for(i in 1:ncol(subset1)){
#   subset1[,i] <- (subset1[,i] - mean(subset1[,i],na.rm = T)) / sd(subset1[,i],na.rm = T)  
# };rm(i)

Breaks <- seq(floor(min(norm_scale)),ceiling(max(norm_scale)), by = 5);Breaks[length(Breaks)] <- max(norm_scale) #by=1
#Breaks <-  round(Breaks,1)
Breaks2 <- as.character(Breaks); Breaks2[length(Breaks2)] <- "Z-Score Normalization(TPM)\n"#""


ann_col <- data.frame(Variety.treatment=samples$CONDITION)
row.names(ann_col) <- samples$SAMPLE

p <- pheatmap::pheatmap(norm_scale,cluster_rows = T,
                        cluster_cols = T,
                        fontsize = 18,
                        annotation_legend = T,#F,
                        labels_row =NULL,
                        annotation_col = ann_col,
                         annotation_colors =  list(Variety.treatment=c(
                                                                Control.AZU.10D="gray",#"#006EC9",
                                                                TAluminum.AZU.10D="purple",
                                                                Control.BGI.10D="salmon",
                                                                TAluminum.BGI.10D="green",
                                                                TAluminum.AZU_I="brown",
                                                                Control.AZU.4HR="black",
                                                                TAluminum.AZU_I.4HR="yellow",
                                                                TAluminum.AZU.4HR ="orange",
                                                                Control.AZU_I.4HR ="lightgreen"
                           
                                                              #Control.AZU="gray",#"#006EC9",
                        #                                      TAluminum.AZU="purple",
                        #                                      Control.BGI="salmon",
                        #                                      TAluminum.BGI="green",
                        #                                      TAluminum.AZU_I="brown",
                        #                                      Control.AZU_I="black"
                                                               )),#"#FF00FF")),
                        drop_levels=F,
                        # display_numbers = matrix_top_go_final_2,#T,#TRU1E,
                        #number_color = "white",
                        fontsize_row = 0.00001,#12
                        fontsize_col =  12,#12
                        border_color="gray98",
                        clustering_distance_rows="correlation",
                        clustering_distance_cols="correlation",
                        fontsize_number=14,#10
                        cellwidth = 100, #48
                        angle_col = 45,
                        cellheight = 8, #16
                        #pvalue
                        legend_breaks = Breaks,
                        legend_labels = Breaks2,
                        #log10
                        #legend_breaks = c(0,15,30,45,60,75,90,93),
                        #legend_labels = c("0","15","30","45","60","75","90","-log10(p-value)\n"),
                        #color
                        #color = hcl.colors(16, "Greens"), #Oranges
                        color = hcl.colors(6, "BluYl"),
                        #color = hcl.colors(8, "BrBg"),
                        #color =rev(hcl.colors(6, "BluYl")),#hcl.colors(10, "BluYl"),
                        legend = T,
                        annotation_names_col = F,
                        na_col = "gray",
                        #number_format= "%.2f",
                        silent = T
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_TPM_NORM.pdf"),
                  width = 50,#20,
                  height = 610 )

#rm(norm_scale,preproc1)


########################################################################################
########################################################################################
########################################################################################
########################################################################################
#Preparing LFC heatmap

x_deg_total <- lapply(1:length(sets),function(i){
  message(i)
  x <- fread(paste0(data_dir,"/",sets_st_cont[[i]][1],"/csv/",
                    "glmQLFTest_",sets_st_cont[[i]][2],"_pval_",
                    pval,"_fulltable.csv"))
  x$status <- NULL
  x$group <- sets[[i]]
  return(x)
})
x_deg_total <- do.call(rbind,x_deg_total)


########################################################################################
#LFC plot pan

subset1a <- x_deg_total[x_deg_total$V1 %in% unique_genes,]
LFC_1a <- as.data.frame(matrix(ncol=length(sets)+1,nrow=length(unique(subset1a$V1))))
LFC_1a[,1] <- unique(subset1a$V1)
colnames(LFC_1a) <- c("gene",sets)
LFC_1a_p_val <- LFC_1a

for(i in 2:ncol(LFC_1a)){

  x_i <- subset1a[which(subset1a$group==sets[i-1]),]
  for(j in 1:nrow(LFC_1a)){
      LFC_1a[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_1a$gene[[j]]),2])
      LFC_1a_p_val[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_1a$gene[[j]]),6])
  };rm(j)
};rm(i)

for(i in 2:ncol(LFC_1a)){
  LFC_1a_p_val[,i] <- as.numeric(LFC_1a[,i])
  LFC_1a[,i] <- as.numeric(LFC_1a_p_val[,i])
  LFC_1a_p_val[,i][which(LFC_1a_p_val[,i]>=pval)] <- ""
  LFC_1a_p_val[,i][which(LFC_1a_p_val[,i]<pval)] <- "*"
  LFC_1a_p_val[,i][which(is.na(LFC_1a_p_val[,i]))] <- ""
  LFC_1a[,i][which(is.na(LFC_1a[,i]))] <- 0
};rm(i)




row.names(LFC_1a) <- LFC_1a$gene;
row.names(LFC_1a_p_val) <- LFC_1a$gene;
LFC_1a$gene <- NULL;LFC_1a_p_val$gene <- NULL
colnames(LFC_1a) <- gsub("COMPLETE/","",colnames(LFC_1a))

ann_col <- data.frame(Variety=c("Azucena.10D","BGI.10D","Azucena.4HR"))
                                #"Azucena_NIL.4HR","Azucena.4HR"))
row.names(ann_col) <- sets
row.names(ann_col) <- gsub("COMPLETE/","",row.names(ann_col))

p <- pheatmap::pheatmap(LFC_1a,cluster_rows = T,
                        cluster_cols = T,
                        fontsize = 18,
                        annotation_legend = T,#F,
                        labels_row =NULL,
                        annotation_col = ann_col,
                        annotation_colors = list(Variety=c(Azucena.10D="purple",
                                                          BGI.10D="brown",
                                                          #Azucena_NIL.4HR= "red",
                                                          Azucena.4HR="gray"
                                                          )),
                        drop_levels=F,
                        fontsize_row = 0.00001,#12
                        fontsize_col =  12,#12
                        border_color=NULL,
                        clustering_distance_rows="correlation",
                        clustering_distance_cols="correlation",
                        fontsize_number=14,#10
                        cellwidth = 150, #48
                        angle_col = 45,
                        cellheight = 8, #16
                        legend_breaks = c(-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,
                                          max(LFC_1a,na.rm = T)),
                        legend_labels = c("-12.5","-10","-7.5","-5","-2.5","0","2.5","5","7.5","10","log2 fold change\n"), 
                        color = hcl.colors(4, "BluYl"),
                        legend = T,
                        annotation_names_col = T,
                        na_col = "gray",
                        silent = T
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_LFC.pdf"),
                  width = 50,
                  height = 600)

################################################################################
################################################################################
#core genes
UP_GENES <- subset1a[which(subset1a$status_name=="UP"),]
up_gen_v <- tapply(UP_GENES$V1, UP_GENES$V1, length)
up_gen_v <- data.frame(genes=names(up_gen_v),count=as.numeric(up_gen_v),status="UP")
up_gen_v <- up_gen_v[which(up_gen_v$count==3),]#3
###
DOWN_GENES <- subset1a[which(subset1a$status_name=="DOWN"),]
down_gen_v <- tapply(DOWN_GENES$V1, DOWN_GENES$V1, length)
down_gen_v <- data.frame(genes=names(down_gen_v),count=as.numeric(down_gen_v),status="DOWN")
down_gen_v <- down_gen_v[which(down_gen_v$count==3),]#3
core_genes <- rbind(up_gen_v,down_gen_v)
################################################################################
subset_core <- x1[which(row.names(x1) %in% core_genes$genes),]
subset_core <- subset_core[rowSums(subset_core[,]) > 0,]
write.table(subset_core,paste0(vhRR_dir,"/atlas_al_core.tsv"),sep = "\t",na = "",row.names = T,col.names = F,quote = F)
write.table(subset_core,paste0(folder_out,"/atlas_al_core.tsv"),sep = "\t",na = "",row.names = T,col.names = T,quote = F)


################################################################################
preproc1 <- preProcess(subset_core, method=c("center", "scale"))
norm_scale <- predict(preproc1, subset_core)

# z-score
# for(i in 1:ncol(subset1)){
#   subset1[,i] <- (subset1[,i] - mean(subset1[,i],na.rm = T)) / sd(subset1[,i],na.rm = T)  
# };rm(i)

Breaks <- seq(floor(min(norm_scale)),ceiling(max(norm_scale)), by = 5);Breaks[length(Breaks)] <- max(norm_scale) #by=1
#Breaks <-  round(Breaks,1)
Breaks2 <- as.character(Breaks); Breaks2[length(Breaks2)] <- "Z-Score Normalization(TPM)\n"#""

ann_col <- data.frame(Variety.treatment=samples$CONDITION)
row.names(ann_col) <- samples$SAMPLE

p <- pheatmap::pheatmap(norm_scale,cluster_rows = T,
                        cluster_cols = T,
                        fontsize = 18,
                        annotation_legend = T,#F,
                        labels_row =NULL,
                        annotation_col = ann_col,
                        annotation_colors = list(Variety.treatment=c(
                          Control.AZU.10D="gray",#"#006EC9",
                          TAluminum.AZU.10D="purple",
                          Control.BGI.10D="salmon",
                          TAluminum.BGI.10D="green",
                          #TAluminum.AZU_I="brown",
                          Control.AZU.4HR="black",
                          #TAluminum.AZU_I.4HR="yellow",
                          TAluminum.AZU.4HR ="orange"
                          #Control.AZU_I.4HR ="lightgreen"
                            
                          #Control.AZU="gray",#"#006EC9",
                          #                                      TAluminum.AZU="purple",
                          #                                      Control.BGI="salmon",
                          #                                      TAluminum.BGI="green",
                          #                                      TAluminum.AZU_I="brown",
                          #                                      Control.AZU_I="black"
                        )),#"#FF00FF")),,
                        drop_levels=F,
                        # display_numbers = matrix_top_go_final_2,#T,#TRU1E,
                        #number_color = "white",
                        fontsize_row = 0.00001,#12
                        fontsize_col =  12,#12
                        border_color="gray98",
                        clustering_distance_rows="correlation",
                        clustering_distance_cols="correlation",
                        fontsize_number=14,#10
                        cellwidth = 50, #48
                        angle_col = 45,
                        cellheight = 8, #16
                        #pvalue
                        legend_breaks = Breaks,
                        legend_labels = Breaks2,
                        #log10
                        #legend_breaks = c(0,15,30,45,60,75,90,93),
                        #legend_labels = c("0","15","30","45","60","75","90","-log10(p-value)\n"),
                        #color
                        #color = hcl.colors(16, "Greens"), #Oranges
                        color = hcl.colors(4, "BluYl"),
                        #color = hcl.colors(8, "BrBg"),
                        #color =rev(hcl.colors(6, "BluYl")),#hcl.colors(10, "BluYl"),
                        legend = T,
                        annotation_names_col = F,
                        na_col = "gray",
                        #number_format= "%.2f",
                        silent = T
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_TPM_NORM_CORE.pdf"),
                  width = 50,#20,
                  height = 360 )

rm(norm_scale,preproc1)

########################################################################################
#LFC plot core

subset1a <- x_deg_total[x_deg_total$V1 %in% core_genes$genes,]
LFC_1a <- as.data.frame(matrix(ncol=length(sets)+1,nrow=length(unique(subset1a$V1))))
LFC_1a[,1] <- unique(subset1a$V1)
colnames(LFC_1a) <- c("gene",sets)
colnames(LFC_1a) <- gsub("COMPLETE/","",colnames(LFC_1a))
LFC_1a_p_val <- LFC_1a

for(i in 2:ncol(LFC_1a)){
  
  x_i <- subset1a[which(subset1a$group==sets[i-1]),]
  for(j in 1:nrow(LFC_1a)){
    LFC_1a[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_1a$gene[[j]]),2])
    LFC_1a_p_val[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_1a$gene[[j]]),6])
  };rm(j)
};rm(i)

for(i in 2:ncol(LFC_1a)){
  LFC_1a_p_val[,i] <- as.numeric(LFC_1a[,i])
  LFC_1a[,i] <- as.numeric(LFC_1a_p_val[,i])
  LFC_1a_p_val[,i][which(LFC_1a_p_val[,i]>=pval)] <- ""
  LFC_1a_p_val[,i][which(LFC_1a_p_val[,i]<pval)] <- "*"
  LFC_1a_p_val[,i][which(is.na(LFC_1a_p_val[,i]))] <- ""
  LFC_1a[,i][which(is.na(LFC_1a[,i]))] <- 0
};rm(i)




row.names(LFC_1a) <- LFC_1a$gene;
row.names(LFC_1a_p_val) <- LFC_1a$gene;
LFC_1a$gene <- NULL;LFC_1a_p_val$gene <- NULL


ann_col <- data.frame(Variety=c("Azucena.10D","BGI.10D","Azucena.4HR"))#"Azucena_NIL.4HR","Azucena.4HR"))
row.names(ann_col) <- sets
row.names(ann_col) <- gsub("COMPLETE/","",row.names(ann_col))

p <- pheatmap::pheatmap(LFC_1a,cluster_rows = T,
                        cluster_cols = T,
                        fontsize = 18,
                        annotation_legend = T,#F,
                        labels_row =NULL,
                        annotation_col = ann_col,
                        annotation_colors = list(Variety=c(Azucena.10D="purple",
                                                           BGI.10D="brown",
                                                           #Azucena_NIL.4HR= "red",
                                                           Azucena.4HR="gray"
                        )),
                        drop_levels=F,
                        fontsize_row = 0.00001,#12
                        fontsize_col =  12,#12
                        border_color=NULL,
                        clustering_distance_rows="correlation",
                        clustering_distance_cols="correlation",
                        fontsize_number=14,#10
                        cellwidth = 300, #48
                        angle_col = 45,
                        cellheight = 8, #16
                        legend_breaks = c(-10,-7.5,-5,-2.5,0,2.5,5,7.5,
                                          max(LFC_1a,na.rm = T)),
                        legend_labels = c("-10","-7.5","-5","-2.5","0","2.5","5","7.5","log2 fold change\n"), 
                        color = hcl.colors(4, "BluYl"),
                        #display_numbers = LFC_1a_p_val,
                        legend = T,
                        annotation_names_col = T,
                        na_col = "gray",
                        silent = T
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_CORE_LFC.pdf"),
                  width = 30,
                  height = 300)

################################################################################

