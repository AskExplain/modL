library(mlbench)
library(rattle.data)
library(pdfCluster)
library(corpcor)
library(mvtnorm)
library(MASS)
library(rrcov)
library(DAAG)


setwd("~/Documents/main_files/GCP/explain-dpm/uploads/script_V5_GVN/")

num_cells <- 1000
dataset_type = "cuomo"
if (dataset_type == "cuomo"){
  load("./data/workspace/cuomo/processed_cuomo_stem_sc_V2.Rdata")
  clusters <- as.factor(cell_metadata_cols$day)
  gene_order <- order(apply(counts,1,var),decreasing=T)


  gene_IDS <- gene_order[1:10000]
  gene_names <- do.call('c',lapply(strsplit(row.names(counts),"_"),function(X){X[2]}))
  genes_interest <- c("SIX2","SIX1","CD31","MEIS1","HOXD11","WT1","PAX2","GATA3","CDH1","JAG1","UMOD","LTL","CUBN","NPHS1",
                      "POU5F1","SOX2","KLF4","MYC")
  gene_IDS <- unique(c(gene_IDS,which(gene_names%in%genes_interest)))


}

path="~/Documents/main_files/GCP/explain-dpm/uploads/script_V1_HDPM/hdpm-master/R/"
setwd(path)
main_dpm <- list.files(path = "./")
lapply(c(main_dpm),function(file_ID){
  source(paste(path,file_ID,sep=""))
})

Y=1
cluster_order <- order(unique(clusters))
y <- counts[gene_IDS,clusters == unique(clusters)[cluster_order[Y]]]
y <- y[,order(apply(y,2,var),decreasing = T)[1:num_cells]]
row.names(y) <- gene_names[gene_IDS]
y <- y[order(row.names(y),decreasing = F),]

library("cate")
cate_fa.em_model <- cate::fa.em(Y = as.matrix(y),r = 2)
hd_c_model <- HDclassif::hddc(data = as.matrix(y),K = 2)

row.names(y)[hd_c_model$class==1]


# y2 <- umap::umap(y,n_components=2)$layout
y30 <- umap::umap(y,n_components=30)$layout
set.seed(5)
main_net5 <- deepgmm(y=y30, k = c(2,2), W = 50,r = c(dim(y30)[2]-1,dim(y30)[2]-2))


main <- main_net5$main
print(table(main))


gene_aim <- do.call('rbind',lapply(sort(unique(main)),function(Y){
  data.frame(class=Y,genes=sort(row.names(y)[main==Y]))
}))
print(unique(main))

gene_aim_list <- lapply(c(sort(unique(main))),function(Y){
  data.frame(class=Y,genes=sort(row.names(y)[main==Y]))
})

names(gene_aim_list) <- sort(unique(main))

plot(y2,col=main,cex=0.3,pch=19)

plot(y2,col=as.factor(main),cex=0.3,pch=19)
par(mfcol=c(3,4))

lapply(c(1:length(unique(main))),function(X){
  plot(y2[main!=X,],col="black",main=X,cex=0.1,pch=19)
  points(y2[main==X,1],y2[main==X,2],col="red",main=X,cex=0.4,pch=19)
  
})

View(gene_aim)
View(gene_aim_list)

# save.image("../../../../N_data/cuomo_DPM_analysis_UMAP_GOOD.rds")





## enable hierarchical structures


main_dist <- as.matrix(dist((main_net5[["comb_main"]]+c(seq(2,56*2,2)))))
row.names(main_dist) <- colnames(main_dist) <- row.names(y)


main_dist[c("POU5F1","SOX2","MYC","NANOG","T","GATA6","EOMES","MIXL1","PIM1","PIM2","PIM3"),
          c("POU5F1","SOX2","MYC","NANOG","T","GATA6","EOMES","MIXL1","PIM1","PIM2","PIM3")]







main_dist <- as.matrix(dist((main_net4[["comb_main"]]+c(seq(2,54*2,2)))))
row.names(main_dist) <- colnames(main_dist) <- row.names(y)


main_dist[c("POU5F1","SOX2","MYC","NANOG","T","GATA6","EOMES","MIXL1","PIM1","PIM2","PIM3"),
          c("POU5F1","SOX2","MYC","NANOG","T","GATA6","EOMES","MIXL1","PIM1","PIM2","PIM3")]



View(main_net4)
View(main_net5)

