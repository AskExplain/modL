setwd("~/Documents/main_files/GCP/explain-dpm/uploads/script_Modular_Learning/CRAN/modL/")

## Step by step process to run modL

## Set up data
# X - Single cell reference data
# Y - NanoString spatial transcript data
# P - Pixels of NanoString spatial transcript data
# R - Pixels of NanoString data outside of spatial transcript regions

x <- read.csv("../../data/external/main/wu_counts.csv",header=T,row.names=1)
y <- read.delim("../../data/external/main/Kidney_Raw_TargetCountMatrix.txt",row.names = 1)
x <- x[,order(colnames(x),decreasing = T)]
p <- read.csv("../../data/workflow/TRAIN_pixels_color.csv",header=T,row.names=1)
r <- read.csv("../../data/workflow/TEST_pixels_color.csv",header=T,row.names=1)

gene_overlap_main <- row.names(x)[row.names(x) %in% row.names(y)]
gene_overlap <- unique(c("ACE","AGT","AGTR1","REN","KRT8","KRT18","PODXL","NPHS1","WT1","MAFB","GPX3","ATF1","RET",gene_overlap_main[gene_overlap_main %in% row.names(x)[order(apply(x,1,var),decreasing = T)[1:300]]]))

## Run corevec
# Inputs:
# P:  ROI images with Known genes
# R:  ROI images with Unknown genes
# Outputs:
# L: gene feature transformation
# Kappa: mapping from unknown spots to known

devtools::load_all("../../scripts/a_library/corevec/draft_v1/corevec/")
p.r_corevec <- corevec::corevec(y = t(as.matrix(p)),x = t(as.matrix(r)),covariate_x = rep(0,dim(r)[1]),covariate_y = rep(0,dim(p)[1]),k_dim = 30,min_iter = 30,max_iter = 30,log = F,factor = F,scale = F)

# x.y_corevec <- corevec::corevec(y = as.matrix(y[gene_overlap,]),x = as.matrix(x[gene_overlap,]),covariate_x = rep(0,dim(x)[2]),covariate_y = rep(0,dim(y)[2]),k_dim = 10,min_iter = 30,max_iter = 30)


## Run Hierarchical Deep Parameterised Mixtures
# Inputs:
# [L . P]: feature transformed ROI images with Known genes
# [L . R]: feature transformed ROI images with Unknown genes
# Outputs:
# Class Probabilities: probabilities for class labels
# Lambda: loadings from mixtures of factor analysis
# A:  projection transform to visualise latent variables in HDPM

library(mlbench)
library(rattle.data)
library(pdfCluster)
library(corpcor)
library(mvtnorm)
library(MASS)
library(rrcov)
library(DAAG)
path="~/Documents/main_files/GCP/explain-dpm/uploads/script_Modular_Learning/scripts/a_library/hdpm/draft_v1/hdpm-master/R/"
setwd(path)
main_dpm <- list.files(path = path)
lapply(c(main_dpm),function(file_ID){
  source(paste(path,file_ID,sep=""))
})
setwd("~/Documents/main_files/GCP/explain-dpm/uploads/script_Modular_Learning/CRAN/modL/")

LP <- t(as.matrix(p.r_corevec$statistics$alpha.L%*%p.r_corevec$statistics$y))
LR <- t(as.matrix(p.r_corevec$statistics$alpha.L%*%p.r_corevec$statistics$x))

best.bic = Inf
for (i in 1:50){
  set.seed(i)
  r.p.hdpm <- deepgmm(y=rbind(LP,LR), k = c(2,2), layers = 1,r = c(dim(LR)[2]-1,dim(LR)[2]-2), scale = T)

  if (Reduce('+',r.p.hdpm$bic)<best.bic){
    best.bic <- Reduce('+',r.p.hdpm$bic)
    best.r.p.hdpm <- r.p.hdpm
  }
}



## Run Three Stage Regression for transfer learning
# Inputs:
# X: reference cell dataset
# Y: observed ST dataset
# Class Probabilities: probabilities for class labels
# Outputs:
# Beta: parameter that maps reference cells to spots per class
# Delta: parameter that maps corevec features to reference cells
# Z: new ST data transferred from known to unknown spots

p.classProb <- best.r.p.hdpm$ps.y.main[1:dim(LP)[1],]
r.classProb <- best.r.p.hdpm$ps.y.main[-c(1:dim(LP)[1]),]

beta <- MASS::ginv(t(as.matrix(x[gene_overlap,]))%*%(as.matrix(x[gene_overlap,])))%*%t(as.matrix(x[gene_overlap,]))%*%as.matrix(y[gene_overlap,])



delta <- array(0,dim=c(dim(LP)[2],dim(beta)[1],dim(r.classProb)[2]))
for (cl in c(1:dim(p.classProb)[2])){
  print(cl)
  delta[,,cl] <- MASS::ginv(t(as.matrix(LP))%*%(as.matrix(LP)))%*%t(as.matrix(p.classProb[,cl]*LP))%*%t(beta)
}



## Run vecnet
# Inputs:
# Y: observed ST dataset
# Z: new ST data transferred from known to unknown spots
# Kappa: mapping from unknown spots to known spots
# Outputs:
# Alpha: gene network parameter from unknown spots to known spots
# Gamma:

source("../../scripts/a_library/vecnet/draft_v1/vecnet.R")

Z <- array(0,dim=c(dim(x[gene_overlap,])[1],dim(LR)[1],dim(r.classProb)[2]))

for (clr in c(1:dim(r.classProb)[2])){
  print(clr)
    Z[,,clr] <- as.matrix(x[gene_overlap,])%*%(t(as.matrix(delta[,,clr]))%*%t(r.classProb[,clr]*LR)) / sum(r.classProb[,clr])
}

alpha <- array(0,dim=c(length(gene_overlap),length(gene_overlap),dim(r.classProb)[2],dim(r.classProb)[2]))
for (clr1 in c(1:dim(r.classProb)[2])){
  for (clr2 in c(1:dim(r.classProb)[2])){
    print(c(clr1,clr2))

    Z.r1.r2_vecnet <- vecnet(y = Z[,rowSums(r.classProb[,clr2]>r.classProb[,-clr2])==dim(r.classProb)[2]-1,clr2], x = Z[,rowSums(r.classProb[,clr1]>r.classProb[,-clr1])==dim(r.classProb)[2]-1,clr1], covariate_y = rep(0,dim(Z[,rowSums(r.classProb[,clr2]>r.classProb[,-clr2])==dim(r.classProb)[2]-1,clr2])[2]), covariate_x = rep(0,dim(Z[,rowSums(r.classProb[,clr1]>r.classProb[,-clr1])==dim(r.classProb)[2]-1,clr1])[2]), min_iter = 30,max_iter = 30, log = F,factor = F)
    main_alpha <- Z.r1.r2_vecnet$statistics$alpha

    colnames(main_alpha) <- row.names(main_alpha) <- gene_overlap
    main_alpha <- (main_alpha)[order(row.names(main_alpha),decreasing = F),order(row.names(main_alpha),decreasing = F)]

    alpha[,,clr1,clr2] <- main_alpha
  }
}

for (clr1 in c(1:dim(r.classProb)[2])){
  for (clr2 in c(1:dim(r.classProb)[2])){
    print(c(clr1,clr2))

    test_alpha <- alpha[,,clr1,clr2]
    colnames(test_alpha) <- row.names(test_alpha) <- gene_overlap
    test_alpha <- (test_alpha)[order(row.names(test_alpha),decreasing = F),order(row.names(test_alpha),decreasing = F)]

    print(c(max(test_alpha["GPX3",]),max(test_alpha[,"MAFB"]),test_alpha["GPX3","MAFB"],test_alpha["MAFB","GPX3"]))
    print(c(max(test_alpha["ATF1",]),max(test_alpha["RET",]),test_alpha["ATF1","RET"],test_alpha["RET","ATF1"]))
  }
}




## Visualisation
all_alpha <- data.frame(rep(row.names(test_alpha),times=length(gene_overlap)),rep(row.names(test_alpha),each=length(gene_overlap)),c(test_alpha))


image(apply(t(apply(matrix(((apply(as.matrix(r.classProb),1,function(X){which(X==max(X))[1]}))),nrow=40,ncol=40,byrow = T),2,rev)),2,function(X){as.integer(as.factor(X))}))

par(mfrow=c(2,2))
image(t(apply(matrix(scale(Z[which(gene_overlap%in%"SPP1"),,2]),nrow=40,ncol=40,byrow = T),2,rev)))
image(t(apply(matrix(scale(Z[which(gene_overlap%in%"GOLGA4"),,1]),nrow=40,ncol=40,byrow = T),2,rev)))
image(t(apply(matrix(scale(Z[which(gene_overlap%in%"CRIM1"),,1]),nrow=40,ncol=40,byrow = T),2,rev)))
image(t(apply(matrix(scale(Z[which(gene_overlap%in%"IGFBP4"),,1]),nrow=40,ncol=40,byrow = T),2,rev)))


activation.atlas <- array(0,dim=c(dim(Z)[2],dim(p.r_corevec$statistics$alpha.L)[2],dim(r.classProb)[2]))

for (clr in c(1:dim(r.classProb)[2])){
  activation.atlas[,,clr] <- t(Z[,,clr])%*%(as.matrix(x[gene_overlap,])%*%t(as.matrix(delta[,,clr])))%*%MASS::ginv(t((as.matrix(x[gene_overlap,])%*%t(as.matrix(delta[,,clr]))))%*%((as.matrix(x[gene_overlap,])%*%t(as.matrix(delta[,,clr])))))%*%(p.r_corevec$statistics$alpha.L)%*%MASS::ginv(t(p.r_corevec$statistics$alpha.L)%*%(p.r_corevec$statistics$alpha.L))
}

visualise_activation.atlas <- array(activation.atlas[,4],dim=c(1600,3,32,32))

library(raster)

for (tile_id in 40){
  print(tile_id)
  plot_tile <- (visualise_activation.atlas[tile_id,,,])
  plot_tile <- aperm((plot_tile+abs(min(plot_tile)))/max(plot_tile+abs(min(plot_tile))),perm = c(2,3,1))
  plot(raster::as.raster(plot_tile))
}










# save data and make a dashboard in django-plotly-dash
# gene.names <- gene_overlap
#
# Z1 <- Z[,,1]
# Z2 <- Z[,,2]
# Z3 <- Z[,,3]
# Z4 <- Z[,,4]
#
# row.names(Z1) <- row.names(Z2) <- row.names(Z3) <- row.names(Z4) <- gene.names
#
# write.csv(Z1,"../../data/workflow/saved_proper/Z1.csv",quote=F,row.names=T)
# write.csv(Z2,"../../data/workflow/saved_proper/Z2.csv",quote=F,row.names=T)
# write.csv(Z3,"../../data/workflow/saved_proper/Z3.csv",quote=F,row.names=T)
# write.csv(Z4,"../../data/workflow/saved_proper/Z4.csv",quote=F,row.names=T)
#
# rcP <- r.classProb
# write.csv(rcP,"../../data/workflow/saved_proper/rcP.csv",quote=F,row.names=T)
#
# aa1 <- activation.atlas[,,1]
# aa2 <- activation.atlas[,,2]
# aa3 <- activation.atlas[,,3]
# aa4 <- activation.atlas[,,4]
#
# write.csv(aa1,"../../data/workflow/saved_proper/aa1.csv",quote=F,row.names=T)
# write.csv(aa2,"../../data/workflow/saved_proper/aa2.csv",quote=F,row.names=T)
# write.csv(aa3,"../../data/workflow/saved_proper/aa3.csv",quote=F,row.names=T)
# write.csv(aa4,"../../data/workflow/saved_proper/aa4.csv",quote=F,row.names=T)
#
#
#
#
# delta1 <- delta[,,1]
# delta2 <- delta[,,2]
# delta3 <- delta[,,3]
# delta4 <- delta[,,4]
#
# write.csv(delta1,"../../data/workflow/saved_proper/delta1.csv",quote=F,row.names=T)
# write.csv(delta2,"../../data/workflow/saved_proper/delta2.csv",quote=F,row.names=T)
# write.csv(delta3,"../../data/workflow/saved_proper/delta3.csv",quote=F,row.names=T)
# write.csv(delta4,"../../data/workflow/saved_proper/delta4.csv",quote=F,row.names=T)
#
#
#
# write.csv(x[gene_overlap,],"../../data/workflow/saved_proper/x.csv",quote=F,row.names=T)
# write.csv(p.r_corevec$statistics$alpha.L,"../../data/workflow/saved_proper/alpha_L.csv",quote=F,row.names=T)
#
#
# for (clr in c(1:4)){
#   main_transform = (as.matrix(x[gene_overlap,])%*%t(as.matrix(delta[,,clr])))%*%MASS::ginv(t((as.matrix(x[gene_overlap,])%*%t(as.matrix(delta[,,clr]))))%*%((as.matrix(x[gene_overlap,])%*%t(as.matrix(delta[,,clr])))))%*%(p.r_corevec$statistics$alpha.L)%*%MASS::ginv(t(p.r_corevec$statistics$alpha.L)%*%(p.r_corevec$statistics$alpha.L))
#   write.csv(main_transform,paste("../../data/workflow/saved_proper/main_transform",clr,".csv",sep=""),quote=F,row.names=T)
#
# }
#
#
#
#
#
# colnames(x) <- gsub("\\.","",gsub(pattern = "[[:digit:]]",replacement = "",colnames(x)))
#
# proper_x <- do.call('cbind',lapply(unique(colnames(x)),function(X){
#   return(
#     (x[gene_overlap,colnames(x) %in% X])[,order(apply(x[gene_overlap,colnames(x) %in% X],2,var),decreasing = T)[1:30]]
#   )
# }))
#
# write.csv(proper_x,"../../data/workflow/saved_proper/proper_x.csv",quote=F,row.names=T)












