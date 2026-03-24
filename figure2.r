####################################################
####------------------figure 2------------------####
####################################################

###code in R###
library(pcutils)
library(tidyverse)
library(ggplot2)
library(viridis)
library(reshape2)
library(colorspace)
library(corrplot)
library(paletteer)

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratObject)

library(ROGUE)
library(reshape2)

###--PAGA connectivities and Jaccard index for qM MuSCs subclusters--###
od<-c("Quiescent-1","Quiescent-2","Activated-1","Activated-2","Activated-3","Transitional","Myoblast-1","Myoblast-2","Myocyte")
data<-as.matrix(read.table("velocity/qM_musc_paga_connectivities_afterbbknn_newcluster.txt",row.names = 1,header=T,sep='\t'))
colnames(data)<-rownames(data)
a<-data[od,]
b<-a[,od]
data<-b
pdf("figure2e_paga.pdf")
corrplot(data,method="pie",type="lower",is.cor=FALSE,cl.ratio=0.2,
         diag=TRUE,col.lim=c(0,1),col=COL1('OrRd', 200),tl.col="black")
dev.off()

dataj<-as.matrix(read.csv("jaccard-qM-MUSC-newlabel.csv",row.names = 1))
colnames(dataj)<-rownames(dataj)
a<-dataj[od,]
b<-a[,od]
dataj<-b
matrix_data <- apply(dataj, c(1, 2), function(x) ifelse(x == 1, 0.5, x))
pdf("figure2e_jaccard.pdf")
corrplot(dataj,method="color",type="upper",is.cor=FALSE,
         col.lim=c(0,1),colorRampPalette(c(rep("#2171b5",1),"white",rep("firebrick",1)))(100),tl.col="black")
dev.off()

###--Proportions of subclusters in LLC and PF muscle groups--###
color_dict <- c(
  "Activated-1"= "#E4A031",
  "Activated-2" = "#F6E2C1",
  "Activated-3" = "#D75B4E",
  "Transitional" = "#A1D4A2",
  "Myoblast-1" = "#7C4D77",
  "Myoblast-2" = "#B6B3D6",
  "Myocyte" = "#26445E",
  "Quiescent-1" = "#DBE4FB",
  "Quiescent-2" = "#1f77b4"
)

data <- read.table('paired_group_name_perct_t.xls', header = TRUE, sep = "\t",row.names=NULL)
rownames(data)<-data$labels
data<-data[,-1]
pdf("figure2f.pdf",width = 7,height=7)
stackplot(data,flow = T,bar_params = list(width = 0.8, position = "stack"),topN = 10,number = T)+
  scale_fill_manual(values = color_dict)+ 
  theme_bw()+ 
  ggtitle("")+ xlab("")+ 
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_blank())
dev.off()      

###--Distribution of log2-transformed fold change in cell abundance--###
parse_h5ad <- function(adata){
  require(reticulate)
  ad <- import("anndata", convert = FALSE)
  ada <- ad$read_h5ad(adata)
  meta <- py_to_r(ada$obs)
  genemeta <- py_to_r(ada$var)
  emb <- py_to_r(ada$obsm['X_umap'])
  pca <- py_to_r(ada$obsm['X_pca'])
  if(class(ada$raw$X)[1] == "scipy.sparse._csr.csr_matrix"){ # or scipy.sparse.csr.csr_matrix
    exp <- t(py_to_r(ada$raw$X$toarray()))
  }else{
    exp <- t(py_to_r(ada$raw$X))
  }
  rownames(exp) <- rownames(py_to_r(ada$raw$var))
  colnames(exp) <- rownames(meta)
  rownames(emb) <- rownames(meta)
  rownames(pca) <- rownames(meta)
  return(
    list(
      metadata = meta,
      genemeta = genemeta,
      expression = exp,
      embedding = emb,
      pca = pca
    )
  )
}
h5ad <- parse_h5ad("qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding
pca <- h5ad$pca

data <- SingleCellExperiment(assays = list(logcounts = expr),
                            colData = meta,
                            reducedDims=SimpleList(PCA=pca,UMAP=emb)) #pca,umap calculated by scanpy

##1Create a Milo object
data_milo <-Milo(data)
##2Construct KNN graph
data_milo <- buildGraph(data_milo, k = 30, d = 30)
##3Defining representative neighbourhoods on the KNN graph
data_milo <- makeNhoods(data_milo, prop = 0.1, k = 30, d=30, refined = TRUE)
##4Counting cells in neighbourhoods
data_milo <- countCells(data_milo, meta.data = as.data.frame(colData(data_milo)), sample="sample")
##5Defining experimental design
data_design <- data.frame(colData(data_milo))[,c("sample", "paired_group_name", "batch")]
#Convert batch info from integer to factor
data_design$batch <- as.factor(data_design$batch) 
data_design <- distinct(data_design)
rownames(data_design) <- data_design$sample
data_design
##6Computing neighbourhood connectivity
data_milo <- calcNhoodDistance(data_milo, d=30)
##7 Do the DA test
da_res <- testNhoods(data_milo, design = ~ paired_group_name , design.df=data_design)
da_res %>% arrange(SpatialFDR) %>% head() 

da_res <- annotateNhoods(data_milo, da_res, coldata_col = "labels")
da_res$labels <- ifelse(da_res$labels_fraction < 0.6, "Mixed", da_res$labels)
da_res$labels <- factor(da_res$labels, levels = rev(c("Quiescent-1","Quiescent-2","Activated-1","Activated-2","Activated-3",
                                                           "Transitional","Myoblast-1","Myoblast-2","Myocyte","Mixed")))
pdf("figure2g.pdf", width = 6, height = 8)
plotDAbeeswarm(da_res, group.by = "labels")
dev.off()

###--ROGUE index for qM MuSCs subclusters--###
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)
rogue.res<-rogue(expr, labels = meta$labels, samples=meta$donor,platform = "UMI", span = 0.6)
res1<-melt(rogue.res)
res1$value<-as.numeric(res1$value)
res1$variable<-factor(res1$variable,levels = c("Quiescent-1","Quiescent-2","Activated-1","Activated-2","Activated-3",
                                               "Transitional","Myoblast-1","Myoblast-2","Myocyte"))
                                                
pdf("figureS2b.pdf",width = 8,height=7)
ggplot(res1,aes(x=variable,y=value,group=variable,color=variable))+
  geom_boxplot()+geom_jitter(width=0.1)+
  theme_bw()+
  scale_color_manual(values = color_dict)+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
  ylab("ROGUE INDEX")+xlab("")+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1,colour = "black",size = 15),
        axis.text.y = element_text(angle=0,hjust = 1,colour = "black",size = 15),
        axis.title.y = element_text(colour = "black",size = 15),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(linewidth =1, colour = "black"))
dev.off()
