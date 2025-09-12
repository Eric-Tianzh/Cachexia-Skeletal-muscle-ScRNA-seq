####################################################
####------------------figure 5------------------####
####################################################

###code in R###
###--The differential activity of 14 signaling pathways between the LLC and PF groups--###
library(Seurat)
library(progeny)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)
options(stringsAsFactors=FALSE)
library(reticulate)

library(openxlsx)
library(readr)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(circlize)
library(ComplexHeatmap)

parse_h5ad <- function(adata){
  require(reticulate)
  ad <- import("anndata", convert = FALSE)
  ada <- ad$read_h5ad(adata)
  meta <- py_to_r(ada$obs)
  genemeta <- py_to_r(ada$var)
  emb <- py_to_r(ada$obsm['X_umap'])
  if(class(ada$raw$X)[1] == "scipy.sparse._csr.csr_matrix"){ # or scipy.sparse.csr.csr_matrix
    exp <- t(py_to_r(ada$raw$X$toarray()))
  }else{
    exp <- t(py_to_r(ada$raw$X))
  }
  rownames(exp) <- rownames(py_to_r(ada$raw$var))
  colnames(exp) <- rownames(meta)
  
  rownames(emb) <- rownames(meta)
  return(
    list(
      metadata = meta,
      genemeta = genemeta,
      expression = exp,
      embedding = emb
    )
  )
}

h5ad <- parse_h5ad("../qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding

data <- CreateSeuratObject(counts = expr, meta.data = meta)
Idents(data)<-data@meta.data$labels

CellsClusters <- data.frame(Cell = names(Idents(data)), 
    CellType = as.character(Idents(data)),
    stringsAsFactors = FALSE)

Cellsgroup <- data.frame(Cell = rownames(data@meta.data), 
    CellGroup = as.character(data@meta.data$paired_group_name),
    stringsAsFactors = FALSE)

##compute the Progeny activity scores and add them to our Seurat object
data <- progeny(data, scale=FALSE, organism="Mouse", top=500, perm=1, return_assay = TRUE)
data <- Seurat::ScaleData(data, assay = "progeny") 

progeny_scores_df <- 
    as.data.frame(t(GetAssayData(data, slot = "scale.data", 
        assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 

progeny_scores_df <- inner_join(progeny_scores_df,Cellsgroup)
  <- inner_join(progeny_scores_df, CellsClusters)
write.csv(progeny_scores_df, "progeny_scores_df.csv", row.names = FALSE)

group="Transitioanl"
result<-read.csv("progeny_scores_df.csv")
result$cellcluster<-paste(result$CellType,result$CellGroup,sep="_")
resultsub<-result[grep(group,result$CellType),]

p<-ggplot(resultsub, aes(
  x = Pathway,             
  y = Activity,
  fill = CellGroup)) +
  geom_boxplot(aes(fill=CellGroup),staplewidth=0.5,outliers = TRUE,outlier.colour = "gray",outlier.alpha = 0.6)+
  geom_hline(yintercept = 0,color="#333333",linetype="dashed")+
  scale_fill_manual(values=c("#8D2F25","#4A5F7E"))+
  stat_compare_means(aes(group = CellGroup),method = "t.test",label = "p.signif")+
  labs(x="",y=group,title = '')+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text(angle=45,vjust =0.6,hjust=0.5),
    panel.grid.major=element_line(colour=NA),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(colour="#4292c6", fill="#4292c6"),
    strip.text.y=element_text(colour="white",face = "bold"))+  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 
pdf("figure5a.pdf",width = 6,height = 8)
print(p)
dev.off()

##Hypoxia related functional genes network##

net<-read.xlsx("model_mouse_full.xlsx")
netsub<-net[which(net$pathway=="Hypoxia"),]
gene<-c("Nfil3","Lox","Nr3c1","Smad3","Itgb1","Pdgfrb","Serpine1","Col1a1","Suco","Ang")
netfilter<-netsub[which(netsub$gene %in% gene),]
netplot<-data.frame("from"=netfilter$pathway,"to"=netfilter$gene,
                    "weight"=netfilter$weight)
netplot$weight<-as.numeric(netplot$weight)
recvinfo<-data.frame("name"=c(netfilter$gene,"Hypoxia"),"pvalue"=c(netfilter$p.value,0))
network <- graph_from_data_frame(d = netplot, directed = F,vertices =recvinfo)

col_fun <- colorRamp2(
  breaks = c(0, 0.05),  
  colors = c('#bd3d3f', 'azure')  
)
my_color <- col_fun(recvinfo$pvalue)

pdf(file = "figure5f.pdf", width = 5, height = 5)
par(bg = "white", mar = c(0, 0, 0, 0))
plot.igraph(network, 
     layout = layout_with_gem,
     vertex.color = my_color,  
     vertex.label.cex = 1.2,  
     vertex.label.color = "black",  
     vertex.frame.color = "transparent",  
     edge.width = E(network)$weight, 
     edge.color="#26445E",
     edge.curved = 0.2
)
lgd <- ComplexHeatmap::Legend(
  col_fun = col_fun,
  at = c(0, 0.01, 0.02, 0.03, 0.04, 0.05),
  title = "pvalue",
  direction = "vertical"
)
draw(lgd, x = unit(0.9, "npc"), y = unit(0.2, "npc")) 
dev.off()

###--Differentially expressed genes between LLC and PF in Transitional state--###
##volcano plot of DEGs in Transitional state##
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
diffresult<-read.table("diffresult_bbknn_newcluster_LLC_PF_raw_ALL.csv",sep = ",",header=T)
gene<-diffresult[which(diffresult$cluster=="Transitional"),]
geneshow<-c("Csmd1","Galnt15","Ctsl","Il6ra","Pvr","Gadd45g","Lox","Prelp","Igfbp5","Crip1","Lgals1")

pdf("figure5b.pdf",width=5,height=5.5)
EnhancedVolcano(gene,
                lab = gene$gene,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = geneshow,
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 1.2,
                labSize = 3,
                boxedLabels=F,
                drawConnectors=T,
                typeConnectors="open",
                arrowheads=F,
                legendPosition = "",
                xlim = c(- 1.5,1.5),
                ylim = c(0, 60),
                axisLabSize = 10,
                titleLabSize = 12,
                col = c("black","black","black","darkred"),
                colAlpha=0.8,
                title = '',
                subtitle = '')
dev.off()

##GSEA analysis of DEGs in Transitional state##
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)
library(ggsci)
library(enrichplot)

hallmark <- read.gmt("mh.all.v2024.1.Mm.entrez.gmt")
gene<-diffresult[which(diffresult$cluster=="Transitional"),]
gene.df <- bitr(gene$gene,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Mm.eg.db)
colnames(gene.df)<-c("gene","ENTREZID")
temp<-inner_join(gene,gene.df,by = "gene")
temp<-temp[order(temp$avg_log2FC,decreasing = T),]
  
genelist<-temp$avg_log2FC
names(genelist)<-temp$ENTREZID
  
gsea.re1<- GSEA(genelist,TERM2GENE = hallmark,pvalueCutoff = 1,pAdjustMethod = 'BH') 
result_GSEA<-as.data.frame(gsea.re1)
result_GSEA<-result_GSEA[order(result_GSEA$NES,decreasing = T),]
result_gsea<-result_GSEA

pdf("figure5c.pdf",width = 5,height = 5)
gseaplot2(gsea.re1,
          geneSetID = rownames(result_gsea)[c(8,41,48,50)],
          title = "",
          color = c("HALLMARK_HYPOXIA"="darkred",
                    "HALLMARK_MYOGENESIS"="#26445E",
                    "HALLMARK_FATTY_ACID_METABOLISM"="#4C7780",
                    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"="#474769"),
          base_size = 14,
          rel_heights = c(1, 0.2, 0.2),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)
dev.off()

###--The activity of 43 cytokines in Transitional state as predicted by Cytosig--###
library(ggplot2)
library(openxlsx)
library(readr)
library(reshape2)
library(stringr)

matrix<-read.csv("celltype_cyto_qm_LLCPF.csv")
rownames(matrix)<-matrix[,1]
matrix<-as.matrix(matrix[,-1])
result<-matrix
result<-melt(result)
colnames(result)<-c("celltype","cytokine","Activity")

group="Transitional_LLC"
resultsub<-result[result$celltype==group,]
resultsub<-resultsub[order(resultsub$Activity,decreasing = T),]
resultsub$cytokine<-str_to_title(tolower(resultsub$cytokine))

pdf("figure5d.pdf",width=8,height=5)
ggplot(resultsub, aes(
  x = factor(cytokine,levels = unique(cytokine)),          
  y = Activity,
  fill = Activity)) +
  geom_bar(stat = 'identity')+
  scale_fill_gradient2(low="#4A5F7E",mid="white", high="#8D2F25")+
  labs(x="",y=group)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust=1),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+                      
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 
dev.off()

###--Transitional-specific regulons--###
library(SCopeLoomR)
library(SCENIC)
library(AUCell)

h5ad <- parse_h5ad("qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding
inputDir='qM/MUSC/pySCENIC/newcluster'
scenicLoomPath=file.path(inputDir,'out_SCENIC.loom')
loom <- open_loom(scenicLoomPath) 

exprMat <- get_dgem(loom)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
save(regulons,"regulons.RData")
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC') #TFs x cells
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)
sub_regulonAUC <- regulonAUC[,match(rownames(meta),colnames(regulonAUC))]
selectedResolution <- "labels" 
cellsPerGroup <- split(rownames(meta),meta[,selectedResolution]) 
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 

#RAS
# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells) rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
#RAS Regulon activity score
write.csv(regulonActivity_byGroup,"RAS.csv")
write.csv(regulonActivity_byGroup_Scaled,"RAS_scaled.csv")

#RSS Regulon Specificity Score
cell=meta[colnames(sub_regulonAUC),selectedResolution]
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cell) 
rss=na.omit(rss) 
save(rss,file = "rss.RData")
write.csv(rss,"RSS.csv")

rssType <- sort(rss[, "Transitional"], decreasing = TRUE)
rssrank <- data.frame(
  regulon = names(rssType), rank = seq_along(rssType),
  rss = rssType
)
rasgroup<-ras[,c("X","Transitional")]
colnames(rasgroup)<-c("regulon","RAS")
rasgroup<-rasgroup[order(rasgroup$RAS, decreasing = TRUE),]
rssras<-merge(rssrank,rasgroup,by.x="regulon",by.y="regulon")
rssras<-rssras[order(rssras$rss,rssras$RAS,decreasing = T),]
write.csv(rssras,"transitional_rssras_mt.csv")

t<-rssras$regulon[1:5]
rssras$regulon[-which(rssras$regulon=="Nr3c1(+)")] <- NA
rssras$regulon[1:5]<-t
rssras$color <- FALSE
rssras$color[which(rssras$regulon=="Nfil3(+)")] <- TRUE
rssras$color[which(rssras$regulon=="Nr3c1(+)")] <- TRUE

pdf("figure5e.pdf",width = 5,height = 5)
ggplot(rssras, aes(x = RAS, y = rss,color=color)) +
  geom_point(size = 2, stroke = 0, shape = 16) +
  scale_color_manual(values = c("grey30","darkred")) +
  ggrepel::geom_text_repel(aes(label = regulon),box.padding = 0.25,size = 3.5) +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.position = "none",
    axis.line = element_line(linetype = 1, color = "black", size = 0.3),
    axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
  ) +
  xlab("Regulon Activity Score") +
  ylab("Regulon Specificity Score") 
dev.off()

