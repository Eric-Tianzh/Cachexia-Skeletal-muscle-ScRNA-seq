####################################################
####------------------figure 4------------------####
####################################################

###code in R###
###--Community probability between tumor and MUSCs subclusters--###
###CellChat
library(CellChat)
library(Seurat)
options(stringsAsFactors = FALSE)

library(openxlsx)
library(viridis)
library(ggplot2)
library(RColorBrewer)

expr<-read.csv("qM-MUSC_after_bbknn_new_cluster_tumor/expr.csv",header=T,row.names=1)
meta<-read.csv("qM-MUSC_after_bbknn_new_cluster_tumor/meta.csv",header=T,row.names=1)
colnames(expr) <- gsub("^X", "", colnames(expr))   
colnames(expr) <- gsub("\\.", "-", colnames(expr))
expr<-as.matrix(expr)

cellchat <- createCellChat(object = expr,meta = meta,group.by = "labels")
CellChatDB <- CellChatDB.mouse 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- updateCellChat(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)    
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
head(df.net)
write.csv(df.net, "cellchat_bbknn/newcluster/cell-cell_communications.qM_MUSC_tumor_all.csv")
save(cellchat,file="cellchat_bbknn/newcluster/qM_MUSC_tumor_all_cellchat.RData")
#visualize the cell-cell communication network
cellchatdata<-read.xlsx("cell-cell_communications.qM_MUSC_tumor_all.xlsx",rowNames = TRUE)
cellchatsecret<-cellchatdata[cellchatdata$annotation=="Secreted Signaling",]
epithsecret<-cellchatsecret[which(cellchatsecret$source == "Epithelial" & cellchatsecret$target != "Epithelial"),]
epithsecret<-epithsecret[order(epithsecret$prob,decreasing = T),]
epithsecret$interaction_name_2<-factor(epithsecret$interaction_name_2,rev(unique(epithsecret$interaction_name_2)))
epithsecret$source<-"Tumor"
epithsecret$LR<-paste(epithsecret$source,epithsecret$target,sep="->")
write.xlsx(epithsecret,"epithsecret.xlsx")

display.brewer.all()
my_palette <- rev(colorRampPalette(brewer.pal(10,"Spectral"))(50))
pdf("figure4a.pdf",width=6,height=6)
ggplot(epithsecret,aes(x=LR,y=interaction_name_2,color=prob))+
  geom_point(shape=19,size=4)+
  scale_color_gradient(low="#F6E2C1", high="#CC5200")+
  theme_bw()+
  labs(x="",y="",,color = "Community Prob")+
  theme(axis.text.x = element_text(angle=90,vjust =1,hjust=0.2),
        panel.grid.minor = element_blank())
dev.off()

##split data in t1 and t2 time points
data <- CreateSeuratObject(counts = expr, meta.data = meta)

t1_cells <- rownames(data@meta.data)[which(data@meta.data$time_group == "t1")]
t1_expression <- GetAssayData(data, slot = "counts")[, t1_cells]
t1_metadata <- data@meta.data[t1_cells, ]
t2_cells <- rownames(data@meta.data)[which(data@meta.data$time_group == "t2")]
t2_expression <- GetAssayData(data, slot = "counts")[, t2_cells]
t2_metadata <- data@meta.data[t2_cells, ]

cellchatt1 <- createCellChat(object = t1_expression,meta = t1_metadata,group.by = "labels")
cellchatt2 <- createCellChat(object = t2_expression,meta = t2_metadata,group.by = "labels")

CellChatDB <- CellChatDB.mouse 
cellchatt1@DB <- CellChatDB
cellchatt2@DB <- CellChatDB

cellchatt1 <- subsetData(cellchatt1) 
cellchatt1 <- identifyOverExpressedGenes(cellchatt1)
cellchatt1 <- identifyOverExpressedInteractions(cellchatt1)
cellchatt2 <- subsetData(cellchatt2) 
cellchatt2 <- identifyOverExpressedGenes(cellchatt2)
cellchatt2 <- identifyOverExpressedInteractions(cellchatt2)

cellchatt1 <- computeCommunProb(cellchatt1, raw.use = TRUE)    
cellchatt1 <- filterCommunication(cellchatt1, min.cells = 10)
df.net <- subsetCommunication(cellchatt1)
write.csv(df.net, "cellchat_bbknn/newcluster/cell-cell_communications_qmmusc_tumor_t1.all.csv")
save(cellchatt1,file="cellchat_bbknn/newcluster/qM_MUSC_tumor_cellchat_t1.RData")
cellchatt2 <- computeCommunProb(cellchatt2, raw.use = TRUE)    
cellchatt2 <- filterCommunication(cellchatt2, min.cells = 10)
df.net <- subsetCommunication(cellchatt2)
write.csv(df.net, "cellchat_bbknn/newcluster/cell-cell_communications_qmmusc_tumor_t2.all.csv")
save(cellchatt2,file="cellchat_bbknn/newcluster/qM_MUSC_tumor_cellchat_t2.RData")

cellchatt1 <- aggregateNet(cellchatt1)
groupSize1 <- as.numeric(table(cellchatt1@idents))
cellchatt2 <- aggregateNet(cellchatt2)
groupSize2 <- as.numeric(table(cellchatt2@idents))

pdf("figure4b.pdf")
netVisual_aggregate(cellchatt1, signaling = "SPP1",weight.scale = T,label.edge= F,edge.color="#26445E",layout ="circle")
netVisual_aggregate(cellchatt2, signaling = "SPP1",weight.scale = T,label.edge= F,edge.label.color="#26445E",layout ="circle")
dev.off()

###--Ligand-receptor pair between tumor and Transitional state--###
###NicheNet
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggpubr)

library(UpSetR)

##1.Read in NicheNet’s networks##
organism <- "mouse"
if(organism == "human"){
  lr_network <- readRDS("/mnt/tzhdata/data/NicheNet_network/lr_network_human_21122021.rds")
  ligand_target_matrix <- readRDS("/mnt/tzhdata/data/NicheNet_network/ligand_target_matrix_nsga2r_final.rds")
  weighted_networks <- readRDS("/mnt/tzhdata/data/NicheNet_network/weighted_networks_nsga2r_final.rds")
} else if(organism == "mouse"){
  lr_network <- readRDS("/mnt/tzhdata/data/NicheNet_network/lr_network_mouse_21122021.rds")
  ligand_target_matrix <- readRDS("/mnt/tzhdata/data/NicheNet_network/ligand_target_matrix_nsga2r_final_mouse.rds")
  weighted_networks <- readRDS("/mnt/tzhdata/data/NicheNet_network/weighted_networks_nsga2r_final_mouse.rds")
}
##2.Read in the expression data of interacting cells##
h5ad <- parse_h5ad("TUMOR_merge_qmmusc_llcandpf.after_batch.epith_aneuploid.group.h5ad")
expr <- h5ad$expression
meta <- h5ad$metadata
emb <- h5ad$embedding
data <- CreateSeuratObject(counts = expr, meta.data = meta)
Idents(data)<-data@meta.data$labels

##3.Perform the NicheNet analysis-step by step##
#3.1.Define a set of potential ligands for both the sender-agnostic and sender-focused approach#
receiver = "Transitional"
expressed_genes_receiver <- get_expressed_genes(receiver, data, pct = 0.05)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver) 
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique() 

sender_celltypes <- c("Epithelial")
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, data, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique() 
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 


#3.2.Define the gene set of interest#
condition_oi <-  "LLC"
condition_reference <- "PF"
seurat_obj_receiver <- subset(data, idents = receiver)              

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "paired_group_name",
                                  min.pct = 0.05) %>% rownames_to_column("gene") 
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene) 
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

seurat_obj_epith <- subset(data, idents = "Epithelial")
DE_table_epith_all <-  FindMarkers(object = seurat_obj_epith,
                                  ident.1 = "t2", ident.2 = "t1",
                                  group.by = "time_group",
                                  logfc.threshold =0,min.pct = 0.05) %>% rownames_to_column("gene")
write.csv(DE_table_epith_all, "geneset_oi_epith_alldiff_in_t2_t1.csv", row.names = FALSE)

#3.3.Define the background genes#
#(All expressed genes in the receiver cell population (that are also in the ligand-target matrix) is defined as the ‘background set’)
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#3.4.Perform NicheNet ligand activity analysis#
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
#Ligands are ranked based on the area under the precision-recall curve (AUPR) 
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
write.csv(ligand_activities, "ligand_activities_epith2inter.csv", row.names = FALSE)

cellchatdata<-read.xlsx("cell-cell_communications.qM_MUSC_tumor_all.xlsx",rowNames = TRUE)
dfgene<-read.csv("geneset_oi_epith_alldiff_in_t2_t1.csv")
nichenetdata<-read.csv("ligand_activities_epith2inter.csv")

nichenetligand<-nichenetdata$test_ligand
cellchatligand<-unique(epithsecret$ligand[epithsecret$target=="Transitional"])
tumordiffgene<-dfgene$gene[dfgene$group2 != "none"]
upset_list <- list(tumordiffgene,nichenetligand,cellchatligand)
names(upset_list)<-c("tumor_diffgene","NicheNet_ligand","CellChat_ligand")
pdf("figure4c.pdf",height=5,width=6)
upset(fromList(upset_list), 
      nsets = 100,   
      nintersects = 40, 
      order.by = "freq", 
      keep.order = F, 
      mb.ratio = c(0.7,0.3),  
      text.scale = 2,
      point.size = 4,line.size = 1,sets.x.label = ""
)
dev.off()


##DEGs in tumor violin plot##
Idents(seurat_obj_epith)<-seurat_obj_epith@meta.data$time_group
sub_seurat_obj_epith <- subset(seurat_obj_epith, idents = c("t1","t2"))
feature<-c("Spp1","Hbegf","Ereg","Adm","Wnt7b")
my_comparisons <- list(c("t1","t2"))
plots_violins <-VlnPlot(sub_seurat_obj_epith, features = feature,group.by="time_group",add.noise=FALSE,pt.size=0,cols=c("#26445E","#E4A031"))

for(i in 1:length(plots_violins)) {
  plots_violins[[i]] <- plots_violins[[i]] +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = 'black',face = "bold", size = 12),
        axis.text.y = element_text(color = 'black', face = "bold",size=12),
        axis.title.y = element_text(color = 'black', face = "bold", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),
        panel.spacing = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic"),
        legend.position = 'none')+
        ylim(0, max(FetchData(seurat_obj_epith, feature[i]), na.rm = TRUE) * 1.4)+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",hide.ns = F,label="p.signif",bracket.size=0.8,tip.length=0,size=6)
}
pdf("figure4d.pdf",height = 8,width = 7)
plots_violins
dev.off()

##NicheNet ligand-target and ligand-receptor heatmaps##
gene<-c("Spp1","Hbegf","Ereg","Adm","Wnt7b")
genetarget<-gene %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = genetarget,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25) 

order_ligands <- intersect(gene, colnames(active_ligand_target_links)) %>% rev()
order_targets <- genetarget$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "#E4A031", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "#E4A031")
pdf("figure4e_down.pdf",width=5,height = 3)
p_ligand_target
dev.off()

genereceptor<-get_weighted_ligand_receptor_links(
  gene, expressed_receptors,
  lr_network, weighted_networks$lr_sig)
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  genereceptor,
  gene,
  order_hclust = "receptors")
p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "#26445E", legend_title = "Prior interaction potential")+
  scale_fill_gradient2(low = "whitesmoke",  high = "#26445E")
pdf("figure4e_top.pdf",width=5,height = 3)
p_ligand_receptor
dev.off()