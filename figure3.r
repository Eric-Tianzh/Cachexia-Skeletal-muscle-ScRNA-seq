####################################################
####------------------figure 3------------------####
####################################################

###code in R###
###--MuSCs trajectory analysis--###
##Monocle3
library(monocle3)
library(ggplot2)
library(viridis)
library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)
library(ggridges)

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

setwd("/mnt/tzhdata/project/cachexia/scRNA/mouse-muscle/qM/MUSC/monocle3")
h5ad <- parse_h5ad("../qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding
gene_annotation <- data.frame(gene_short_name = rownames(expr))
rownames(gene_annotation) <- rownames(expr)

cds <- new_cell_data_set(expr, cell_metadata = meta, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,preprocess_method = "PCA")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- emb
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
saveRDS(cds, "bbknn/monocle3_qM_newcluster.rds")


# a helper function to identify the root principal points:
###!!!!!!!!!!!!!!!!!!!!!!!!!!change the time bin and cluster label
get_earliest_principal_node <- function(cds, time_bin="Quiescent-1"){ #start node
  cell_ids <- which(colData(cds)[, "labels"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
root_pr_nodes=get_earliest_principal_node(cds)
cds1<-order_cells(cds,root_pr_nodes=c(root_pr_nodes))#qm bbknn newcluster

p1<-plot_cells(cds1, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_roots = FALSE, label_branch_points = FALSE,
           cell_size = 1.5,trajectory_graph_color = "white",trajectory_graph_segment_size = 1)+
  theme(axis.text.x = element_text(angle=0,vjust = 1,colour = "black",size = 15),
        axis.text.y = element_text(angle=0,hjust = 1,colour = "black",size = 15),
        axis.title.x = element_text(colour = "black",size = 15),
        axis.title.y = element_text(colour = "black",size = 15))+
  scale_color_viridis(discrete=FALSE, name="pseudotime")
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 14, name = "Spectral")))
ggsave(p1, file='figure3b-left.pdf', width=10, height=10)

#pseudotime ridge plot
ps_tim <- data.frame(pseudotime=pseudotime(cds1),cell=names(pseudotime(cds1)))
meta<-data.frame(pData(cds))
metacluster<-data.frame(cluster=meta$labels,cell=rownames(meta))
psdata<-merge(metacluster,ps_tim)
psdata<-psdata[!is.infinite(psdata$pseudotime),]
psdata$cluster<-factor(psdata$cluster,levels = c("Quiescent-1","Quiescent-2","Activated-1","Activated-2","Activated-3",
                                               "Transitional","Myoblast-1","Myoblast-2","Myocyte"))

p2<-ggplot(psdata, aes(x = pseudotime, y = cluster, fill = cluster)) +
  geom_density_ridges(scale = 2) +
  theme_ridges(font_size = 16, grid = TRUE,center_axis_labels=TRUE) +
  #scale_fill_brewer(palette = 4)
  scale_fill_viridis(discrete=TRUE,name="cluster")+
  theme(legend.position = "None")
ggsave(p2,file='figure3b-right.pdf', width=10, height=10)

gene<-c('Pax7','Myf5','Myod1')
pg<-plot_genes_in_pseudotime(cds1[gene,], color_cells_by="labels", 
                         min_expr=0.5, ncol = 2)+scale_color_manual(values = color_dict)
ggsave(pg, file='figureS3c.pdf', width=8, height=8)

###--MuSCs trajectory analysis quiescence score (Q-score) and the myogenic score (M-score)--###
library(reticulate)
library(AUCell)
library(ggpubr)
h5ad <- parse_h5ad("../../qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding

gene<-list("Q"=c("Pax3","Pax7","Myf5"),"M"=c("Myod1","Fosl1","Srxn1"))
cells_rankings <- AUCell_buildRankings(expr, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(gene, cells_rankings)
selectedThresholds <- rowMeans(getAUC(cells_AUC))

colnames(emb)<-c("umap1","umap2")
cells_AUC_data <- as.data.frame(t(cells_AUC@assays@data$AUC))
cells_AUC_data<-cells_AUC_data[rownames(emb),]  
data <- cbind(emb,cells_AUC_data)
data$R<-(data$M+1)/(data$Q+1)
data$labels<-meta$labels
write.csv(data,"qm_musc_emb_QMSCORE.csv",row.names = T)

pdf("figure3c_left.pdf", width = 5, height = 5) 
ggplot(data, aes(umap1, umap2))  +
  geom_point(aes(colour  = Q),size=0.5) +
  scale_color_gradient(low="#f7ebe7",high="#8D2F25")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("figure3c_middle.pdf", width = 5, height = 5) 
ggplot(data, aes(umap1, umap2))  +
  geom_point(aes(colour  = M),size=0.5) +
  scale_color_gradient(low="#d1dbe4",high="#4A5F7E")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("figure3c_right.pdf", width = 5, height = 5) 
ggplot(data, aes(umap1, umap2))  +
  geom_point(aes(colour  = R),size=0.5) +
  scale_color_gradient2(low="#8D2F25",mid="white",high="#4A5F7E",midpoint = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("figureS3d.pdf", width = 5, height = 5) 
ggplot(data, aes(
  x = labels,           
  y = R, 
  color=labels)) +
  geom_boxplot(outliers = FALSE)+
  scale_fill_manual(values = color_dict)+scale_color_manual(values = color_dict)+
  theme_bw()+
  ylab("")+xlab("")+
  theme(axis.text.x = element_text(angle=90,vjust = 1,hjust = 1,colour = "black",size = 15),
        axis.text.y = element_text(angle=0,hjust = 1,colour = "black",size = 15),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(linewidth =1, colour = "black"))
dev.off()

###--MuSCs differentiation potency scores analysis by CytoTRACE2--###
library(CytoTRACE2)
library(dplyr)
library(reticulate)
library(ggpubr)
setwd("/mnt/tzhdata/project/cachexia/scRNA/mouse-muscle/qM/MUSC")
h5ad<-parse_h5ad("qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding
annotation <- data.frame(phenotype = meta$labels) %>% #qM-MUSC newcluster
  set_rownames(., rownames(meta))

results = cytotrace2(expr,species = "mouse",ncores = 30,seed=1234)
saveRDS(results, "./cytotrace/bbknn/newcluster/results.rds")

plots <- plotData(cytotrace2_result = results,expr,
                  annotation = annotation, 
                  is_seurat = FALSE)

results<-readRDS('results.rds')
group<-data.frame(meta$labels)
group$cell<-rownames(meta)
results$cell<-rownames(results)
cytoresult<-merge(group,results,by.x='cell',by.y='cell')
cytoresult$meta.labels<-factor(cytoresult$meta.labels,levels = c("Quiescent-1","Quiescent-2","Activated-1","Activated-2","Activated-3",
                                               "Transitional","Myoblast-1","Myoblast-2","Myocyte"))

p1<-ggboxplot(cytoresult, x="meta.labels", y="CytoTRACE2_Score", width = 0.8, 
          color = "black",#轮廓颜色
          fill="meta.labels",#填充
          palette = color_dict,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.3, #误差条大小
          size=0.8, #箱型图边线的粗细
          outlier.shape=NA #不显示outlier
          )+
  geom_jitter(aes(fill=meta.labels,color=meta.labels),width=0.03,alpha = 0.5,shape=21,size=1)+
  scale_fill_manual(values = color_dict)+
  scale_color_manual(values = color_dict)+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
  ylab("Potency score")+xlab("")+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1,colour = "black",size = 15),
        axis.text.y = element_text(angle=0,hjust = 1,colour = "black",size = 15),
        axis.title.y = element_text(colour = "black",size = 15),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(linewidth =1, colour = "black"))

my_comparisons <- list(c("Myoblast-2","Myocyte"),c("Myoblast-1","Myoblast-2"),c("Transitional","Myoblast-1"),
                       c("Activated-3","Transitional"),c("Activated-2","Activated-3"),
                       c("Activated-1","Activated-2"),c("Quiescent-2","Activated-1"),c("Quiescent-1","Quiescent-2")
                       )
pdf("figure3d.pdf",width=12,height = 10)
p1+stat_compare_means(comparisons = my_comparisons,label="p.signif",step.increase=0.1,
                      method = "wilcox.test")
dev.off()

###--GO biological process pathways--###
library(GOplot)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)

##Transitional
data<-read.csv("../qm_musc_res1.0_labels_marker.csv")
datagroup<-data[which(data$group=="Transitional"),]
gene<-datagroup[,c(3,5)]
colnames(gene)<-c("ID","logFC")

enrich<-read.csv("00GO-FCfilter-up/clusterTransitional_go.csv")
enrich<-data.frame(enrich$ONTOLOGY,enrich$ID,enrich$Description,enrich$p.adjust,enrich$geneID)
colnames(enrich)<-c("category","ID","term","adj_pval","genes")
enrich$genes<-gsub("/",",",enrich$genes)

circ <- circle_dat(enrich,gene)
library(stringr)
circ$genes<-str_to_title(tolower(circ$genes))
#specific term
term <- c("regulation of collagen biosynthetic process",
          "leukocyte activation involved in immune response",
          "response to hypoxia",
          "response to lipopolysaccharide")
chord <- data.frame(chord_dat(circ,gene,term))

chord<-chord[order(chord[,'logFC'],decreasing = T),]
chord1<-t(chord %>% dplyr::select(-logFC))
genefc<-circ[which(circ$term %in% term),c("term","genes","logFC")]
genefc<-genefc[order(genefc$logFC,decreasing = T),]
col_fun <- colorRamp2( #inter
  breaks = c(-1, 0, 1),  # foldchange范围
  colors = c('cornflowerblue', 'azure', 'brown1')  # 对应颜色
)
receiver_colors <- col_fun(genefc$logFC)
names(receiver_colors) <- genefc$genes
grid.col = c("regulation of collagen biosynthetic process" = "#26445E", 
             "leukocyte activation involved in immune response" = "#C76B60", 
             "response to hypoxia" = "#257D8B",
             "response to lipopolysaccharide"="#474769",receiver_colors)

pdf("figure3e-right.pdf",width = 6,height = 6)
circos.par(start.degree = 180, clock.wise = TRUE)
chordDiagram(chord1, grid.col = grid.col,
             #transparency=0.2,
             transparency=0.4,
             directional = 1, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             diffHeight = 0.05,
             link.target.prop=FALSE,
            # scale=T,
             #annotationTrack = c("name","grid"),
             annotationTrack = c("grid"),
             annotationTrackHeight = mm_h(c(2, 2)),
             preAllocateTracks = 1
            )
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 
              CELL_META$ylim[1], 
              CELL_META$sector.index,
              facing = "clockwise", 
              niceFacing = TRUE, 
              adj = c(0, 0.3), 
              cex = 0.8
  )

}, bg.border = NA)
lgd <- ComplexHeatmap::Legend(
  col_fun = col_fun,
  title = "log2FC",
  direction = "vertical"
)
draw(lgd, x = unit(0.9, "npc"), y = unit(0.1, "npc")) 
dev.off()
circos.clear()

data<-read.csv("../qm_musc_res1.0_labels_marker.csv")
datagroup<-data[which(data$group=="Activated-3"),]
gene<-datagroup[,c(3,5)]
colnames(gene)<-c("ID","logFC")

enrich<-read.csv("00GO-FCfilter-up/clusterActivated-3_go.csv")
enrich<-data.frame(enrich$ONTOLOGY,enrich$ID,enrich$Description,enrich$p.adjust,enrich$geneID)
colnames(enrich)<-c("category","ID","term","adj_pval","genes")
enrich$genes<-gsub("/",",",enrich$genes)

circ <- circle_dat(enrich,gene)
library(stringr)
circ$genes<-str_to_title(tolower(circ$genes))
#specific term
term <- c("muscle organ development","muscle cell proliferation",
          "muscle system process")
chord <- data.frame(chord_dat(circ,gene,term))

##activated-3
chord<-chord[order(chord[,'logFC'],decreasing = T),]
chord1<-t(chord %>% dplyr::select(-logFC))
genefc<-circ[which(circ$term %in% term),c("term","genes","logFC")]
genefc<-genefc[order(genefc$logFC,decreasing = T),]
col_fun <- colorRamp2( #act
  breaks = c(-2, 0, 2),  # foldchange范围
  colors = c('cornflowerblue', 'azure', 'brown1')  # 对应颜色
)
receiver_colors <- col_fun(genefc$logFC)
names(receiver_colors) <- genefc$genes
grid.col = c("muscle organ development" = "#3E90BF", 
             "muscle cell proliferation" = "#4C7780", 
             "muscle system process" = "#474769",receiver_colors)

pdf("figure3e-left.pdf",width = 6,height = 6)
circos.par(start.degree = 180, clock.wise = TRUE)
chordDiagram(chord1, grid.col = grid.col,
             #transparency=0.2,
             transparency=0.4,
             directional = 1, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             diffHeight = 0.05,
             link.target.prop=FALSE,
            # scale=T,
             #annotationTrack = c("name","grid"),
             annotationTrack = c("grid"),
             annotationTrackHeight = mm_h(c(2, 2)),
             preAllocateTracks = 1
            )
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 
              CELL_META$ylim[1], 
              CELL_META$sector.index,
              facing = "clockwise", 
              niceFacing = TRUE, 
              adj = c(0, 0.3), 
              cex = 0.8
  )

}, bg.border = NA)
lgd <- ComplexHeatmap::Legend(
  col_fun = col_fun,
  title = "log2FC",
  direction = "vertical"
)
draw(lgd, x = unit(0.9, "npc"), y = unit(0.1, "npc")) 
dev.off()
circos.clear()

###--Gene average expression heatmap--###
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ComplexHeatmap)
h5ad <- parse_h5ad("qM-MUSC.after_batch.bbknn.newcluster.marker.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding
data <- CreateSeuratObject(counts = expr, meta.data = meta)
Idents(data)<-data@meta.data$labels
data@meta.data$group<-paste(data@meta.data$labels,data@meta.data$paired_group_name,sep = "_")

fibrogene<-c('Itgb2','Tgfb1','Lox','Itgb1','Itgav','Smad2','Smad3','Pdgfrb','Ccn2','Prelp','Serpine1','Col1a1','Col1a2')
myogenesis<-c('Igf1','Igfbp7','Myf5','Fhl1','Tnnt1','Sparc','Jsrp1')
inflam<-c('Hif1a','Il6ra','Cxcl2','Bcl3','Pvr','Ccl7','Ahr','Osmr','Ier5')

#计算平均表达量
gene_cell_exp <- AverageExpression(data,
                                   features = c(fibrogene,myogenesis,inflam),
                                   group.by = 'group',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[,colnames(gene_cell_exp) %in% c("Transitional_LLC","Activated-3_LLC","Transitional_PF","Activated-3_PF")]


#complexheatmap作图
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = color_dict))#颜色设置

df<-data.frame(class=rep('Fibrosis',length(fibrogene)))
df<-rbind(df,data.frame(class=rep('Myogenesis',length(myogenesis))))
df<-rbind(df,data.frame(class=rep('Inflammation',length(inflam))))
rownames(df)<-c(fibrogene,myogenesis,inflam)
df$class<-factor(df$class,levels=c('Fibrosis','Myogenesis','Inflammation'))

right_anno = HeatmapAnnotation(df = df,
                             border = F,
                             simple_anno_size=unit(1,'mm'),annotation_width=unit(0.1, 'mm'),
                             show_annotation_name = F,
                             gp = gpar(col = NA),which = "row",
                             col = list(class=c('Fibrosis'="#719AAC",'Myogenesis'="#E29135",'Inflammation'="#72B063")))


#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
ht<-Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = T,
        show_column_names = T,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("#2171b5","#d8d8d8","darkred"))(100),
        border = '#d8d8d8',
        rect_gp = gpar(col = "#d8d8d8", lwd = 0.2),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        row_split = df$class,
        row_title = NULL,
        row_gap = unit(1, "mm"),
        right_annotation = right_anno)

pdf("figure3f.pdf",width=3,height=7)
draw(ht)
dev.off()

