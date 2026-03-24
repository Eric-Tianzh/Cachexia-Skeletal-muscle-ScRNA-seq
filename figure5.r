####################################################
####------------------figure 5------------------####
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
pdf("figure5a.pdf",width=6,height=6)
ggplot(epithsecret,aes(x=LR,y=interaction_name_2,color=prob))+
  geom_point(shape=19,size=4)+
  scale_color_gradient(low="#F6E2C1", high="#CC5200")+
  theme_bw()+
  labs(x="",y="",,color = "Community Prob")+
  theme(axis.text.x = element_text(angle=90,vjust =1,hjust=0.2),
        panel.grid.minor = element_blank())
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
pdf("figure5b.pdf",height = 8,width = 7)
plots_violins
dev.off()
