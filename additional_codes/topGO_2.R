

topGO_gen <- function(ta_go_file,custom_format,genes_id,pval,out_dir,output_name){
  
#calling libraries  
require(data.table);require(topGO);require(dplyr);
require(GO.db);require(parallel)

if(!file.exists(paste0(out_dir,"/",output_name,".tsv"))){
message("Reading GOA file")
#reading goa file downloaded from plaza
ta_go_file <- fread(ta_go_file)

if(custom_format==F){
  #changing "#gene_id" colname to "gene_id"
#  colnames(ta_go_file)[1] <- unlist(strsplit(colnames(ta_go_file)[1],"#"))[[2]]
  #subsetting to reduce size
  ta_go_file <- ta_go_file[,c(1,3)]#.(gene_id,go,description)]
  colnames(ta_go_file) <- c("gene_id","go")
} else {
#  colnames(ta_go_file)<- c("go","gene_id","na")
  #subsetting to reduce size
  ta_go_file <- ta_go_file[,c(2,1)]#.(gene_id,go)]  
  colnames(ta_go_file) <- c("gene_id","go")
}

#indexing background genes to subset and do not use the complete genome
index_b <- ta_go_file$gene_id %in% unique(ta_go_file$gene_id)
x_back_vect <- ta_go_file[index_b,]

#creating a list of gene ids and their GO terms     
gene_2_GO <- unstack(x_back_vect[,c(2,1)])

message("Formatting to topGO format")
#remover genes sin anotacion
keep <- genes_id %in% unique(ta_go_file$gene_id)
keep <- which(keep==TRUE)
candidate_list <- genes_id[keep]
#converting gene list to factors 
geneList <- factor(as.integer(unique(ta_go_file$gene_id) %in% candidate_list),levels = c(0,1))
names(geneList) <- unique(ta_go_file$gene_id)

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
allGO1a <- allGO1
#allGO1a <-allGO1[allGO1 %in% ta_go_file[,2]]


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

message("Subsetting only significant results")
#subsetting only significant
all_res1_s <- all_res1[which(all_res1$status_name=="significant"),]

#Only apply filtering steps if there are enriched results!

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

write.table(all_res1_s,paste0(out_dir,"/",output_name,".tsv"),row.names = F,na = "",sep = "\t")

} else {
  message("Already process, loading")
  all_res1_s <- fread(paste0(out_dir,"/",output_name,".tsv"),header=T)
}
return(all_res1_s)

}


#######################
#running for core
ta_go_file <- "D:/TESIS_PHD/CHAPTER3/GO_PLAZA/go.osa.csv.gz" # plaza filtered
custom_format =F #format
#getting gene ids
genes_id_file <- "D:/TESIS_PHD/CHAPTER3/DEG/atlas_al_core.tsv"
#genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_419.tsv"
#genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_52.tsv"
genes_id_file <- data.table::fread(genes_id_file)
genes_id <- genes_id_file$`#`

pval <- 0.05
out_dir <- "D:/TESIS_PHD/CHAPTER3/DEG" #outpur dir
output_name <- "GO_CORE" #outcome filename

x <- topGO_gen(ta_go_file,custom_format,genes_id,pval,out_dir,output_name)
#######################
#running for pan
ta_go_file <- "D:/TESIS_PHD/CHAPTER3/GO_PLAZA/go.osa.csv.gz" # plaza filtered
custom_format =F #format
#getting gene ids
genes_id_file <- "D:/TESIS_PHD/CHAPTER3/DEG/atlas_al_pan.tsv"
#genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_419.tsv"
#genes_id_file <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input/atlas_drought_52.tsv"
genes_id_file <- data.table::fread(genes_id_file)
genes_id <- genes_id_file$`#`

pval <- 0.05
out_dir <- "D:/TESIS_PHD/CHAPTER3/DEG" #outpur dir
output_name <- "GO_PAN" #outcome filename

x <- topGO_gen(ta_go_file,custom_format,genes_id,pval,out_dir,output_name)


