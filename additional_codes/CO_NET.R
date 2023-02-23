# Load EGAD and the data files 
library(EGAD);require(data.table);require(igraph)
#https://bioc.ism.ac.jp/packages/3.12/bioc/vignettes/EGAD/inst/doc/EGAD.html
data_dir <- "D:/TESIS_PHD/CHAPTER3/vHHR/input"
core_genes <- paste0(data_dir,"/","atlas_al_core.tsv")
core_files <- fread(core_genes,header = F)
row.names(core_files) <- core_files$V1;genes <- core_files$V1
core_files$V1 <- NULL
core_files <- core_files[rowSums(core_files)>0]
core_files2 <- core_files
core_files2 <- log2(core_files2)

core_files <- t(core_files); colnames(core_files) <- genes
core_files <- cor(core_files)
#core_files[abs(core_files)>=0.7] <- 1
row.names(core_files) <- genes;colnames(core_files) <- genes;

#diag(core_files) <- NA
#core_files[lower.tri(core_files)] <- NA
core_files[abs(core_files)<0.7] <- 0
core_files[abs(core_files)>=0.7] <- 1
#core_files[abs(core_files)<0.7] <- NA
# core_files <- core_files[rowSums(core_files,na.rm = T)>0,]
# core_files <- core_files[,colSums(core_files,na.rm = T)>0]



x2 <- igraph::graph_from_adjacency_matrix(core_files,mode = "undirected")
components <- igraph::clusters(x2, mode="weak")
biggest_cluster_id <- which.max(components$csize)

# ids
vert_ids <- V(x2)[components$membership == biggest_cluster_id]

# subgraph
x2 <- igraph::induced_subgraph(x2, vert_ids)


xx2 <- cluster_louvain(x2, weights = NULL, resolution = 1)

V(x2)$Louvain <- membership(xx2)
V(x2)$degree <- degree(x2)
V(x2)$betweenness <- betweenness(x2)

igraph::write.graph(x2,
                    file = paste0("D:/TESIS_PHD/CHAPTER3", "/", "COR.graphml"),
                    format = "graphml")


data.frame <- get.data.frame(x2, what= c("vertices") )
write.table(data.frame,
                    file = paste0("D:/TESIS_PHD/CHAPTER3", "/", "node_table.tsv"),na = "",sep = "\t",row.names = F)


data <- lapply(1:nrow(core_files),function(i){
  #i <- 1
  data_i <- data.frame(gene1 = colnames(core_files)[i],
                     gene2 = colnames(core_files)[which(core_files[1,]==1)])
  return(data_i)
})

data <- do.call(rbind,data)
data <- data[!duplicated(data[ , c("gene1", "gene2")]),]



coexp_net <- build_coexp_network(core_files2, genes, method="spearman", flag="rank")
colnames(coexp_net) <- genes;row.names(coexp_net) <- genes;

#load annotations
GO_files <- fread(paste0(data_dir,"/","goosa.tsv"),header = F)
GO_files$V3 <- NULL


# Store your annotation matrix
goterms <- unique(GO_files$V1)
annotations <- make_annotations(GO_files[,c(2,1)],genes,goterms)


hist <- plot_distribution(node_degree(coexp_net), b=100,
                          xlim=c(0,14), ylim=c(0,2), xlab="Node degree")
# 
# semantic_net <- build_semantic_similarity_network(annotations, genes)
# hist <- plot_distribution(node_degree(semantic_net), b=20, xlab="Node degree", bars=TRUE)

gba_auc_nv <- neighbor_voting(annotations, coexp_net, nFold=3, output="AUROC")
# or 
gba_pr_nv <- neighbor_voting(annotations, coexp_net, nFold=3, output="PR")
head(gba_auc_nv)
